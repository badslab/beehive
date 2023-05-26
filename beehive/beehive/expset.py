import logging
from functools import lru_cache, partial
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import polars as pl
import yaml

import beehive.exceptions as bex
from beehive import util
from beehive.util import find_prq, get_geneset_db

WARNED_NO_GENE_COL = False

lg = logging.getLogger(__name__)
lg.setLevel(logging.INFO)


diskcache = partial(
    util.diskcache, where=util.get_datadir("cache"), refresh=True)


# NOTE:
# datasets now depend on the view
# if all views are active at the same time, DATASETS should not be a global
# variable.
# @lru_cache(1)

@lru_cache(16)
def get_views(dsid: str):
    datadir = util.get_datadir("h5ad")
    yamlfile = datadir / f"{dsid}.yaml"
    with open(yamlfile, "r") as F:
        y = yaml.load(F, Loader=yaml.SafeLoader)
    return y["use_in_view"]


@lru_cache(32)
def get_datasets(has_de: bool = False,
                 view_name: Optional[str] = None):
    """Return a dict with all dataset."""

    datasets: Dict[str, Dict] = {}
    datadir = util.get_datadir("h5ad")
    lg.debug(f"checking h5ad files in {datadir}")

    for yamlfile in datadir.glob("*.yaml"):
        use = True
        basename = yamlfile.name
        basename = basename.replace(".yaml", "")
        lg.debug(f"Considering dataset {yamlfile}")
        with open(yamlfile, "r") as F:
            y = yaml.load(F, Loader=yaml.SafeLoader)
            if (view_name is not None) \
                    and (view_name not in y["use_in_view"]):
                use = False
            authors = y["author"].split(",")
            authors = [x.strip() for x in authors]
            y["first_author"] = authors[0]
            if len(authors) < 3:
                y["short_author"] = y["author"]
            else:
                y["short_author"] = authors[0] + " .. " + authors[-1]

            if "short_title" not in y:
                if len(y["title"]) >= 65:
                    y["short_title"] = y["title"][:57] + "..."
                else:
                    y["short_title"] = y["title"]
        if use is True:
            lg.debug("Using dataset {yamlfile}f for {view_name}")
            datasets[basename] = y

    if has_de:
        # return only datasets with diffexp data
        DSDE = {a: b for (a, b) in datasets.items()
                if len(b.get("diffexp", {})) > 0}
        # lg.debug(
        #    f"expset datadir is {datadir}, found {len(DSDE)} "
        #    f"(out of {len(DATASETS)}) sets with DE data"
        # )
        return DSDE
    else:
        lg.debug(f"expset datadir is {datadir}, found {len(datasets)} sets")
        for a, b in datasets.items():
            lg.debug(f"  - {a}")
        return datasets


def get_dataset_siblings(dsid: str, view_name: str) -> dict:
    """Return datasets with the same study_id."""
    dsets = get_datasets(view_name=view_name)
    dset = dsets[dsid]
    study_id = dset["study"]
    siblings = {}
    for d, dd in dsets.items():
        if dd["study"] == study_id:
            siblings[d] = dd
    return siblings


# @lru_cache(32)
def get_dataset(dsid, view_name):
    """Return metadata on a single dataset."""
    rv = get_datasets(view_name=view_name)[dsid]
    lg.info(f"Returning dataset {dsid}")
    return rv


def get_facet_options(dsid: str,
                      only_categorical: bool = False,
                      include_skip=False, view_name: str = "") \
        -> List[Tuple[str, str]]:
    """
    Return obs columns that a dataset can be facetted on.
    These have to be categorical

    param dataset - dictionary as loaded from the .yaml files
    """
    dataset = get_dataset(dsid, view_name)
    rv = []
    for key, data in dataset.get('obs_meta', {}).items():
        use = False
        if only_categorical:
            if data.get('dtype') == 'categorical':
                use = True
            if data.get('dtype') == 'skip' and include_skip:
                use = True
        else:
            if data.get('dtype') == 'categorical':
                use = True
            elif data.get('dtype') == 'numerical' \
                    and data.get('facet_on', True):
                use = True

        if use:
            name = data.get('name', key)
            rv.append((key, name))
    return list(sorted(rv))


def get_legend_of_obs(dsid, meta):
    datadir = util.get_datadir("h5ad")
    legend = ""
    for yamlfile in datadir.glob("*.yaml"):
        basename = yamlfile.name.replace(".yaml", "")
        if dsid == basename:
            with open(yamlfile, "r") as F:
                y = yaml.load(F, Loader=yaml.SafeLoader)
                try:
                    legend = y["obs_meta"][meta]["legend"]
                except:  # noqa: E722
                    continue
    return legend


def get_facet_options_numerical(dsid: str, view_name: str = "") \
        -> List[Tuple[str, str]]:
    """
    Return obs columns that a dataset can be facetted on.
    These have to be categorical

    param dataset - dictionary as loaded from the .yaml files
    """
    dataset = get_dataset(dsid, view_name)
    rv = []
    for key, data in dataset.get('obs_meta', {}).items():
        use = False
        if data.get('dtype') == 'numerical':
            use = True
        if use:
            name = data.get('name', key)
            rv.append((key, name))
    return list(sorted(rv))


def get_gene_meta_multi_aggregate(
        dsid: str, gene: str,
        metafields: List[str],
        view_name: str = "",
        nobins: int = 8,
        remove_NONE: bool = True):

    """Aggregate."""
    data = get_gene_meta_multi(
        dsid=dsid, gene=gene,
        metafields=metafields,
        view_name=view_name, nobins=nobins,
        remove_NONE=remove_NONE)

    # functions instead of lambda's yield better names in the
    # pandas aggregate..

    def q01(x): return np.quantile(x, 0.25)
    def q25(x): return np.quantile(x, 0.25)
    def q75(x): return np.quantile(x, 0.25)
    def q99(x): return np.quantile(x, 0.9)

    aggdata = data.groupby(metafields)['_gene_'].agg(
        [np.mean, np.median, max, min, np.std,
         q01, q25, q75, q99
        ])
    return aggdata


def get_gene_meta_multi(
        dsid: str, gene: str,
        metafields: List[str],
        view_name: str = "",
        nobins: int = 8,
        remove_NONE: bool = True):
    """Get gene and metadata for multiple fields."""

    genedata = get_gene(dsid, gene)
    if genedata is None:
        raise bex.GeneNotFoundException

    rvr = {"_gene_": genedata[gene]}

    for i, meta in enumerate(metafields):
        if meta == '--':
            continue
        m = get_meta(dsid, meta, nobins=nobins,
                     view_name=view_name, raw=True).to_pandas()
        rvr[meta] = m[meta]

    rv = pd.DataFrame(rvr)

    # TODO: This is probably not the fastest method
    if remove_NONE:
        rv = rv[rv.apply(lambda x: 'NONE' not in list(x), axis=1)]

    return rv

def get_gene_meta_three_facets(
        dsid: str, gene: str, meta1: str, meta2: str, meta3: str = "mouse.id",
        nobins: int = 8, view_name: str = "",mean_option: int = 0):
    """"
    function to return the aggregate of 3 metadata options
    first two metadatas will be displayed in a boxplot doubled
    and the third will be specific to every possible combination
    we end up with means of the third metadata
    and means of (metadata1 * metadata2) groups.
    """

    genedata = get_gene(dsid, gene)
    if genedata is None:
        raise bex.GeneNotFoundException

    # TODO: Rewrite this whole function? or better document

    # gdseries = genedata.to_pandas()[gene]
    # metanew = get_meta_multi(dsid, [meta1, meta2, meta3], nobins=nobins,
    #                         view_name=view_name, raw=True)
    # metanew.columns = 'meta_' + metanew.columns
    # metanew['gene_' + gene] = gdseries

    # print(metanew.head(3).T)

    emptyDF = False

    metadata = get_meta(dsid, meta1, nobins=nobins,
                        view_name=view_name, raw=True)  # facet1
    metadata = metadata.select((pl.col(meta1)).cast(str))

    if genedata.columns == metadata.columns:
        metadata.columns = [metadata.columns[0] + "_category"]
        new_meta = metadata.columns[0]
    else:
        new_meta = meta1

    groupby_columns1 = [new_meta]
    groupby_columns2 = [new_meta]

    if meta2 == "--":
        metadata2 = pl.DataFrame([])
        new_meta2 = meta2
        emptyDF = True

    else:
        metadata2 = get_meta(dsid, meta2, nobins=nobins,
                             view_name=view_name, raw=True)  # facet2
        metadata2 = metadata2.select((pl.col(meta2)).cast(str))

        if genedata.columns == metadata2.columns:
            metadata2.columns = [metadata2.columns[0] + "_category"]
        new_meta2 = metadata2.columns[0]

    if meta3 == "--":
        metadata3 = pl.DataFrame([])
        new_meta3 = meta3
    else:
        metadata3 = get_meta(dsid, meta3, nobins=nobins,
                             view_name=view_name, raw=True)
        metadata3 = metadata3.select((pl.col(meta3)).cast(str))
        if genedata.columns == metadata3.columns:
            metadata3.columns = [metadata3.columns[0] + "_category"]
        new_meta3 = metadata3.columns[0]

    if new_meta2 != "--":
        groupby_columns1 = groupby_columns1 + [new_meta2]
        groupby_columns2 = groupby_columns2 + [new_meta2]
    if new_meta3 != "--":
        groupby_columns2 = groupby_columns2 + [new_meta3]

    if genedata is None or metadata is None \
            or metadata2 is None \
            or metadata3 is None:
        return None


    rv_combined = (
        pl.concat([genedata, metadata, metadata2, metadata3], how="horizontal")
        .groupby(groupby_columns1)
        .agg(
            [
                pl.count(),
                pl.mean(gene).alias("mean"),
                pl.std(gene).alias("std"),
                pl.median(gene).alias("median"),
                pl.quantile(gene, 0.25).alias("q25"),
                pl.quantile(gene, 0.75).alias("q75"),
                pl.quantile(gene, 0.01).alias("q01"),
                pl.quantile(gene, 0.99).alias("q99"),
                pl.min(gene).alias("min"),
                pl.max(gene).alias("max"),
            ]
        )
    )
    
    if meta3 == "--":
        rv_mouseid = pl.DataFrame([])
    else:
        if mean_option == 0:
            rv_mouseid = (
                pl.concat([genedata, metadata, metadata2, metadata3], how="horizontal")
                .groupby(groupby_columns2)
                .agg(
                    [
                        pl.count().alias("count_" + new_meta3),  # count_mouseid
                        pl.mean(gene).alias("mean_" + new_meta3),  # mean_mouseid
                    ]
                )
            )
        else:
            rv_mouseid = (
                pl.concat([genedata, metadata, metadata2, metadata3], how="horizontal")
            .groupby(groupby_columns2)
            .agg(
                [
                    pl.count().alias("count_" + new_meta3),  # count_mouseid
                    pl.median(gene).alias("median_" + new_meta3),  # mean_mouseid
                ]
            ) 
            )  

    gc1 = groupby_columns1[0]
    go1 = get_order_of_obs(dsid, meta1)

    if not(go1): #no ordering available
                #check if they are digits.
        list_of_gc1_labels = list(filter(lambda x: x != "NONE",list(metadata.unique())[0]))
        if all(elem.isdigit() for elem in list_of_gc1_labels):
                #sort it numerically
            list_of_gc1_labels = sorted(list_of_gc1_labels, key=int)
        else:
            #sort it alphanumerically
            list_of_gc1_labels = sorted(list_of_gc1_labels)

        go1 = {value: index for index, value in enumerate(list_of_gc1_labels)}

    go1x = max(go1.values()) + 2 if go1 else 2

    rv_combined = rv_combined.with_columns([
        pl.col(gc1).apply(lambda x: go1.get(x, go1x)).alias('order1')])

    if len(groupby_columns1) == 1:
        rv_combined = rv_combined.rename({"order1": "order"})
    else:
        gc2 = groupby_columns1[1]
        go2 = get_order_of_obs(dsid, gc2)
        if not(go2): #no ordering available
                #check if they are digits.
            list_of_gc2_labels = list(filter(lambda x: x != "NONE",list(metadata2.unique())[0]))
            if all(elem.isdigit() for elem in list_of_gc2_labels):
                    #sort it numerically
                list_of_gc2_labels = sorted(list_of_gc2_labels, key=int)
            else:
                #sort it alphanumerically
                list_of_gc2_labels = sorted(list_of_gc2_labels)

            go2 = {value: index for index, value in enumerate(list_of_gc2_labels)}


        go2x = max(go2.values()) + 2 if go2 else 2

        rv_combined = rv_combined.with_columns([
            pl.col(gc2).apply(lambda x: go2.get(x, go2x)).alias('order2')])
        rv_combined = rv_combined.with_columns([
            ((rv_combined['order2'] * go2x) + rv_combined['order1'])
            .alias('order')])
        rv_combined = rv_combined.drop(["order1", "order2"])

    # switch to dfs
    df_combined = rv_combined.to_pandas()
    if meta3 != "--":
        df_mouseid = rv_mouseid.to_pandas()

    # calculate other measures
    df_combined['perc'] = 100 * df_combined['count'] / \
        df_combined['count'].sum()
    df_combined['_segment_top'] = df_combined['q99']
    df_combined['_bar_top'] = df_combined['q75']
    df_combined['_bar_median'] = df_combined['median']
    df_combined['_bar_bottom'] = df_combined['q25']
    df_combined['_segment_bottom'] = df_combined['q01']

    # manipulation to get the ["X"] columns as the factors
    if not (emptyDF):
        df_combined["x"] = df_combined[[
            new_meta2, new_meta]].apply(tuple, axis=1)
        df_combined.sort_values(
            by="x", key=lambda col: col.str[0], inplace=True)
        df_combined.sort_values(
            by="x", key=lambda col: col.str[1], inplace=True)
        
        if meta3 != "--":
            df_mouseid["x"] = df_mouseid[[
                new_meta2, new_meta]].apply(tuple, axis=1)
            df_mouseid.sort_values(
                by="x", key=lambda col: col.str[0], inplace=True)
            df_mouseid.sort_values(
                by="x", key=lambda col: col.str[1], inplace=True)

    # merge on the just created ["X"] column
        if meta3 != "--":
            final_rv = pd.merge(df_combined, df_mouseid, on="x")
        else:
            final_rv = df_combined
        # sort..
        final_rv.sort_values(by="x", key=lambda col: col.str[0], inplace=True)
        final_rv.sort_values(by="x", key=lambda col: col.str[1], inplace=True)

        final_rv = final_rv.rename(columns={"x": "cat_value"})

        # sort the column for colors to be displayed..
        if meta3 != "--":
            final_rv.sort_values(by=f'{new_meta}_y', inplace=True)
        else:
            final_rv.sort_values(by=f'{new_meta}', inplace=True)
    else:
        if meta3 != "--":
            final_rv = pd.merge(df_combined, df_mouseid, on=new_meta)
        else:
            final_rv = df_combined
        final_rv = final_rv.rename(columns={new_meta: "cat_value"})

    return final_rv


# might be better to merge these two functions into one:
# get_colors_of_obs()
# get_order_of_obs()
def get_colors_of_obs(dsid: str, meta: str,special = False,view_name = ""):
    final_dict = {}
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    backup_palette = plt.cm.Dark2
    # TODO need to incorporate it into yaml... cant just have it here, special for one dataset/one facet.
    #backup_palette2 = ["#ce4d19","#e4561c","#fa5e1f","#ff7e33","#ff931f”,”#FFA029","#ffad33","#ffbf60","#ffcf87","#ffd79a","#ffdead","#FFE8C2"]
    
    ii = 0

    datadir = util.get_datadir("h5ad")
    for yamlfile in datadir.glob("*.yaml"):
        basename = yamlfile.name.replace(".yaml", "")
        if dsid == basename:
            with open(yamlfile, "r") as F:
                y = yaml.load(F, Loader=yaml.SafeLoader)
                #if there are no values in the yaml, means that there is no coloring..
                if y["obs_meta"][meta].get("values"):
                    for key, data in y["obs_meta"][meta]["values"].items():
                        name = key
                        color = data.get("color")
                        # default is grey
                        if color is None:
                            final_dict[name] = mpl.colors.to_hex(
                                backup_palette(ii))
                            ii += 1
                        else:
                            final_dict[name] = color
                else:
                    metadata = get_meta(dsid, meta, view_name=view_name, raw=True)
                    metadata = metadata.select((pl.col(meta)).cast(str))
                    list_of_gc1_labels = list(filter(lambda x: x != "NONE",list(metadata.unique())[0]))
                    for key in list_of_gc1_labels:
                        final_dict[key] = mpl.colors.to_hex(backup_palette(ii))
                        ii+=1

    return final_dict


def get_order_of_obs(dsid: str, meta: str):
    final_dict: Dict[str, Any] = {}
    datadir = util.get_datadir("h5ad")
    if meta == "--": ##only one level to show in x axis..
        return final_dict
    for yamlfile in datadir.glob("*.yaml"):
        basename = yamlfile.name.replace(".yaml", "")
        if dsid == basename:
            with open(yamlfile, "r") as F:
                y = yaml.load(F, Loader=yaml.SafeLoader)
                if "category" in meta:
                    meta = meta.split("_category")[0]
                if y["obs_meta"][meta].get("values"):
                    for key, data in y["obs_meta"][meta]["values"].items():
                        name = key
                        order = data.get("order")
                        # there is no order in the yaml..
                        
                        if not (order):
                            return final_dict
                        final_dict[name] = order
    return final_dict


def get_gene_meta_agg(dsid: str, gene: str, meta: str, nobins: int = 8):
    """Return gene and observation."""

    genedata = get_gene(dsid, gene)
    metadata = get_meta(dsid, meta, nobins=nobins)

    if genedata is None:
        return None
    if metadata is None:
        return None

    # in new yamls, there is sometimes names of column same as that of gene
    # e.g. in man2m.. meta = APOE, gene = APOE
    # this ruins the concatination..
    if genedata.columns == metadata.columns:
        metadata.columns = [metadata.columns[0] + "_category"]

    new_meta = metadata.columns[0]

    rv = (
        pl.concat([genedata, metadata], how="horizontal")
        .groupby(new_meta)
        .agg(
            [
                pl.count(),
                pl.mean(gene).alias("mean"),
                pl.std(gene).alias("std"),
                pl.median(gene).alias("median"),
                pl.quantile(gene, 0.25).alias("q25"),
                pl.quantile(gene, 0.75).alias("q75"),
                pl.quantile(gene, 0.01).alias("q01"),
                pl.quantile(gene, 0.99).alias("q99"),
                pl.min(gene).alias("min"),
                pl.max(gene).alias("max"),
            ]
        )
    )

    rv = rv.to_pandas()
    rv = rv.rename(columns={new_meta: "cat_value"})
    try:
        rv2 = rv.copy()
        rv2['cat_value'] = rv['cat_value'].astype(float)
        rv2 = rv2.sort_values(by=['cat_value', 'mean'])
        rv = rv.loc[list(rv2.index)]

    except ValueError:
        rv = rv.sort_values(by=["cat_value", "mean"])

    return rv


def get_gene(dsid, gene):
    """Return expression values for this dataset."""

    try:
        rv = pl.read_parquet(find_prq(dsid, 'X'), columns=[gene])
    except pl.exceptions.SchemaError:
        return None
    return rv


def get_defields(dsid, view_name=None):
    ds = get_dataset(dsid, view_name)
    dex = ds.get("diffexp")
    return [(a, b.get('name', a))
            for (a, b) in dex.items()]


def get_gsea_columns(dsid, column=None):
    gsdb = get_geneset_db(dsid)
    sql = '''
        SELECT DISTINCT(column)
        FROM gsea '''
    rv = pd.read_sql(sql, gsdb)
    return list(rv['column'])


def get_gsea_data(
        dsid, defield=None,
        gsfilter="",
        sort_on='nes',
        return_no=50,):

    assert sort_on in ['fdr', 'nes', 'absnes']

    gsdb = get_geneset_db(dsid)

    sql = ''' SELECT
                gsea.geneset_hash,
                gsea.nes,
                gsea.fdr,
                gsea.lead_genes,
                abs(gsea.nes) as abs_nes,
                gsea.column as de_column,
                gs.title as geneset_title,
                gs.type as geneset_type,
                gs.direction as direction,
                gs.genes as genes,
                gr.organism as organism,
                gr.study_title as study_title,
                gr.study_hash,
                gr.study_author as study_author,
                gr.study_year as study_year
            FROM genesets as gs, groups as gr,
                    gsea as gsea
            WHERE gs.group_hash = gr.group_hash
                AND gsea.geneset_hash = gs.geneset_hash '''

    if defield != '-':
        sql += f""" AND gsea.column = "{defield}" """

    if gsfilter != "":
        sql += f""" AND gs.title LIKE "%{gsfilter}%" """

    if sort_on == 'fdr':
        sql = sql + ' ORDER BY gsea.fdr '
    elif sort_on == 'nes':
        sql = sql + ' ORDER BY gsea.nes DESC '
    elif sort_on == 'absnes':
        sql = sql + ' ORDER BY abs_nes DESC '

    sql += ' LIMIT 100 '
    rv = pd.read_sql(sql, gsdb)

    rv['genes'] = rv['genes'].str.split()
    rv['no_genes'] = rv['genes'].apply(len)
    rv['lead_genes'] = rv['lead_genes'].str.split()
    rv['no_lead_genes'] = rv['lead_genes'].apply(len)

    return rv


def get_dedata_simple(dsid, field,
                      ranktype: str ='lfc'):

    global WARNED_NO_GENE_COL
    cols = get_varfields(dsid)

    if 'gene' in cols:
        gcol = 'gene'
    else:
        gcol = cols[-1]
        if not WARNED_NO_GENE_COL:
            lg.debug(f'Could not find a gene column, using: {gcol}')
            WARNED_NO_GENE_COL = True

    if ranktype == 'slp':
        padj_col = field[:-3] + 'padj'
        slp_col = field[:-3] + 'slp'
        rv = pl.read_parquet(find_prq(dsid, 'var'),
                             [gcol, field, padj_col])
        rv = rv.to_pandas().rename(columns={gcol: 'gene'})
        rv[slp_col] = \
            (-1 * np.log10(rv[padj_col])).clip(0.001,100) \
            * rv[field]
        rv = rv[['gene', slp_col]]

    else:
        rv = pl.read_parquet(find_prq(dsid, 'var'), columns=[gcol, field])
        rv = rv.to_pandas().rename(columns={gcol: 'gene'})

    return rv


def get_dedata(dsid, categ, genes, view_name: str = ""):
    """Return diffexp data."""

    ds = get_dataset(dsid, view_name)
    dex = ds.get("diffexp")
    assert categ in dex

    if isinstance(genes, str):
        genes = [genes]

    rvpl = pl.read_parquet(find_prq(dsid, 'var'),
                           columns=["field"] + genes)
    rv = rvpl.to_pandas()
    rvx = rv["field"].str.split("__", expand=True)
    rvx.columns = ["categ", "cat_value", "measurement"]
    rv = pd.concat([rv, rvx], axis=1)
    rv = rv[rv["categ"] == categ].copy()
    del rv["categ"]
    del rv["field"]

    rv = rv.pivot(index="cat_value", columns="measurement", values=genes)
    return rv


## TODO: next two functions can  be merged into 1?


def get_dedata_new(dsid, categ):

    # to get 'gene' column

    # TODO: gene column MUST be called 'gene'
    ### SOME datasets have it called field
    ### Some datasets have it called gene
    ### some datasets have it at the last columns with __index__level..
    ## loss of ".index" when switched to prqs on old and new datasets.

    last_col = len(pl.read_parquet(find_prq(dsid, 'var')).columns)

    try:
        rv1 = pl.read_parquet(find_prq(dsid, 'var'), columns=["gene"])
    except:
        try:
            rv1 = pl.read_parquet(find_prq(dsid, 'var'), columns=["field"])
        except:
            rv1 = pl.read_parquet(find_prq(dsid, 'var'), columns=[last_col - 1])
    # get logfoldchange and padjusted for one category
    # (example injection__None)
    rv2 = pl.read_parquet(find_prq(dsid, 'var'),
                          columns=[categ + "__lfc"] + [categ + "__padj"])

    rv = pl.concat([rv2, rv1], how="horizontal")

    rv = rv.to_pandas()
    # in new format, padjusted and lfc are already on columns, genes on rows.
    # get genes has same index as that for the rows.
    return rv


def get_dedata_quadrant(dsid, categ1, categ2):

    # to get 'gene' column
    last_col = len(pl.read_parquet(find_prq(dsid, 'var')).columns)
    rv1 = pl.read_parquet(find_prq(dsid, 'var'), columns=["gene"])
    
    # get logfoldchange and padjusted for one category
    # (example injection__None)
    rv2 = pl.read_parquet(
        find_prq(dsid, 'var'),
        columns=[categ1 + "__lfc"] + [categ2 + "__lfc"] +
                [categ1 + "__padj"] + [categ2 + "__padj"])

    rv = pl.concat([rv2, rv1], how="horizontal")

    # avoid where chosen lfc's are the same ones..
    rv = rv.rename({categ1 + "__lfc": categ1 + "__lfc_1",
                    categ2 + "__lfc": categ2 + "__lfc_2",
                    categ1 + "__padj": categ1 + "__padj_1",
                    categ2 + "__padj": categ2 + "__padj_2"})
    rv = rv.to_pandas()
    # in new format, padjusted and lfc are already on columns, genes on rows.
    # get genes has same index as that for the rows.
    return rv1


def get_dedata_abundance(dsid, categ):

    # to get 'gene' column
    last_col = len(pl.read_parquet(find_prq(dsid, 'var')).columns)

    rv1 = pl.read_parquet(find_prq(dsid, 'var'), columns=[last_col - 1])

    rv2 = pl.read_parquet(
        find_prq(dsid, 'var'),
        columns=[categ + "__lfc"] + [categ + "__padj"] + [categ + "__lcpm"])

    rv = pl.concat([rv2, rv1], how="horizontal")
    rv = rv.to_pandas()

    return rv


def get_defaults(dsid, view_name: str = ""):
    final_dict = {}
    datadir = util.get_datadir("h5ad")
    for yamlfile in datadir.glob("*.yaml"):
        basename = yamlfile.name.replace(".yaml", "")
        if dsid == basename:
            with open(yamlfile, "r") as F:
                y = yaml.load(F, Loader=yaml.SafeLoader)
                # find which view_name:
                try:
                    for viewidx in range(len(y["defaults"])):
                        if view_name == list(y["defaults"][viewidx].keys())[0]:
                            final_dict = y["defaults"][viewidx]
                except:  # noqa: E722
                    return final_dict
    return final_dict


def get_meta_multi(dsid: str,
                   cols: List[str],
                   raw: bool = False,
                   nobins: int = 8,
                   view_name: str = ""):

    find_m = {}
    for c in cols:
        find_m[c] = get_meta(dsid, col=c, raw=raw, nobins=nobins,
                             view_name=view_name).to_series()
    rv = pd.DataFrame.from_dict(find_m)
    return rv


def get_meta(dsid, col, raw=False, nobins=8,
             view_name: str = ""):
    """Return one obs column."""
    rv = pl.read_parquet(find_prq(dsid, 'obs'), columns=[col])
    
    if raw:
        # just return whatever is in the db.
        return rv

    ds = get_dataset(dsid, view_name)

    dscol = ds["obs_meta"][col]

    if dscol["dtype"] == "categorical":
        rv[col] = rv[col].cast(str)

    elif dscol["dtype"] == "numerical":
        rvq = pd.qcut(rv.to_pandas()[col], nobins,
                      duplicates="drop", precision=2)

        rvcat = pd.DataFrame(
            dict(no=range(1, len(rvq.cat.categories) + 1),
                 q=rvq.cat.categories)
        ).set_index("q")

        rvcat["cic"] = rvq.value_counts()
        rvcat["cic"] = (100 * rvcat["cic"]) / rvcat["cic"].sum()

        rvcat = rvcat.reset_index()

        rvcat["name"] = rvcat.apply(
            lambda r: f"{r['no']:02d} {r['q']} - {r['cic']:.1f}%".format(**r),
            axis=1
        )

        rvq = rvq.cat.rename_categories(list(rvcat["name"]))

        rv[col] = rvq.astype(str)

    return rv


# @diskcache()
def get_genes(dsid):
    """Return a list fo genes for this datset."""
    X = pl.scan_parquet(find_prq(dsid, 'X'))
    return X.columns


# @diskcache()
def obslist(dsid):
    lg.warning("`obslist` is deprecated function, use get_obsfields")
    return get_obsfields(dsid)


def get_obsfields(dsid):
    """Return a list of obs columns for this datset."""
    X = pl.scan_parquet(find_prq(dsid, 'obs'))
    return X.columns


def get_varfields(dsid):
    """Return a list of var columns for this datset."""
    var = find_prq(dsid, 'var')
    
    X = pl.scan_parquet(find_prq(dsid, 'var'))
    return X.columns


def units_of_gene_expression(dsid):
    datadir = util.get_datadir("h5ad")
    for yamlfile in datadir.glob("*.yaml"):
        basename = yamlfile.name.replace(".yaml", "")
        if dsid == basename:
            with open(yamlfile, "r") as F:
                y = yaml.load(F, Loader=yaml.SafeLoader)
                # find which view_name:
                return y.get("unit_gene_expression", "")
