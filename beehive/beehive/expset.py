from ast import Pass
from functools import partial, lru_cache
from itertools import groupby
import logging
from typing import Dict, List, Tuple


import pandas as pd
import numpy as np
import polars as pl
import yaml

from beehive import util


lg = logging.getLogger(__name__)
lg.setLevel(logging.DEBUG)




diskcache = partial(
    util.diskcache, where=util.get_datadir("cache"), refresh=True)

#NOTE:
#datasets now depend on the view
#if all views are active at the same time, DATASETS should not be a global variable.
# @lru_cache(1)
def get_datasets(has_de: bool = False,view_name: str = ""):
    """Return a dict with all dataset."""
    # DATASETS: Dict[str, Dict] = {}
    DATASETS = dict()

    datadir = util.get_datadir("h5ad")
    #global DATASETS
    if len(DATASETS) == 0:
        for yamlfile in datadir.glob("*.yaml"):
            use = True
            basename = yamlfile.name
            basename = basename.replace(".yaml", "")
            with open(yamlfile, "r") as F:
                y = yaml.load(F, Loader=yaml.SafeLoader)
                if view_name not in y["use_in_view"]:
                    use = False
                authors = y["author"].split(",")
                authors = [x.strip() for x in authors]
                y["first_author"] = authors[0]
                if len(authors) < 3:
                    y["short_author"] = y["author"]
                else:
                    y["short_author"] = authors[0] + ".." + authors[-1]

                if "short_title" not in y:
                    if len(y["title"]) >= 65:
                        y["short_title"] = y["title"][:57] + "..."
                    else:
                        y["short_title"] = y["title"]
            if use == True:
                DATASETS[basename] = y

    if has_de:
        # return only datasets with diffexp data
        DSDE = {a: b for (a, b) in DATASETS.items()
                if len(b.get("diffexp", {})) > 0}
        lg.info(
            f"expset datadir is {datadir}, found {len(DSDE)} "
            f"(out of {len(DATASETS)}) sets with DE data"
        )
        return DSDE
    else:
        lg.info(f"expset datadir is {datadir}, found {len(DATASETS)} sets")
        return DATASETS


def get_dataset_siblings(dsid: str,view_name: str) -> dict:
    """Return datasets with the same study_id."""
    dsets = get_datasets(view_name = view_name)
    dset = dsets[dsid]
    study_id = dset["study"]
    siblings = {}
    for d, dd in dsets.items():
        if dd["study"] == study_id:
            siblings[d] = dd
    return siblings


# @lru_cache(32)
def get_dataset(dsid,view_name):
    """Return metadata on a single dataset."""
    rv = get_datasets(view_name = view_name)[dsid]
    lg.info(f"Returning dataset {dsid}")
    return rv


def get_facet_options(dsid: str,
                      only_categorical: bool = False, include_skip = False, view_name: str = "") -> List[Tuple[str, str]]:
    """
    Return obs columns that a dataset can be facetted on.
    These have to be categorical

    param dataset - dictionary as loaded from the .yaml files
    """
    dataset = get_dataset(dsid,view_name)
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

def get_legend_of_obs(dsid,meta):
    datadir = util.get_datadir("h5ad")
    legend = ""
    for yamlfile in datadir.glob("*.yaml"):
        basename = yamlfile.name.replace(".yaml", "")
        if dsid == basename:
            with open(yamlfile, "r") as F:
                y = yaml.load(F, Loader=yaml.SafeLoader)
                try:
                    legend = y["obs_meta"][meta]["legend"]
                except:
                    continue
    return legend

def get_facet_options_numerical(dsid: str, view_name: str = "") -> List[Tuple[str, str]]:
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

def get_gene_meta_three_facets(dsid: str, gene: str, meta1: str, meta2: str, meta3: str = "mouse.id",nobins:int = 8,view_name: str = ""):
    """"
    function to return the aggregate of 3 metadata options
    first two metadatas will be displayed in a boxplot doubled
    and the third will be specific to every possible combination
    we end up with means of the third metadata
    and means of (metadata1 * metadata2) groups.
    """
    emptyDF = False

    genedata = get_gene(dsid,gene)
    metadata = get_meta(dsid, meta1, nobins=nobins,view_name = view_name, raw = True) #facet1
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
        metadata2 = get_meta(dsid,meta2,nobins=nobins,view_name = view_name, raw = True) #facet2
        metadata2 = metadata2.select((pl.col(meta2)).cast(str))
        if genedata.columns == metadata2.columns:
            metadata2.columns = [metadata2.columns[0] + "_category"]
        new_meta2 = metadata2.columns[0]

    if meta3 == "--":
        metadata3 = pl.DataFrame([])
        new_meta3 = meta3
    else:
        metadata3 = get_meta(dsid,meta3,nobins=nobins,view_name = view_name, raw = True) #most likely mouse_id
        metadata3 = metadata3.select((pl.col(meta3)).cast(str))
        if genedata.columns == metadata3.columns:
            metadata3.columns = [metadata3.columns[0] + "_category"]
        new_meta3 = metadata3.columns[0]

    if new_meta2 != "--":
        groupby_columns1 = groupby_columns1 + [new_meta2]
        groupby_columns2 = groupby_columns2 + [new_meta2]
    if new_meta3 != "--":
        groupby_columns2 = groupby_columns2 + [new_meta3]


    if genedata is None or metadata is None or metadata2 is None or metadata3 is None:
        return None
    
    rv_combined = (
        pl.concat([genedata, metadata,metadata2,metadata3], how="horizontal")
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
    rv_mouseid = (
        pl.concat([genedata, metadata,metadata2,metadata3], how="horizontal")
        .groupby(groupby_columns2)
        .agg(
            [
                pl.count().alias("count_" + new_meta3), #count_mouseid
                pl.mean(gene).alias("mean_" + new_meta3), #mean_mouseid
            ]
        )
    )

    #switch to dfs
    rv_combined = rv_combined.to_pandas()
    rv_mouseid = rv_mouseid.to_pandas()

    #calculate other measures
    rv_combined['perc'] = 100 * rv_combined['count'] / rv_combined['count'].sum()
    rv_combined['_segment_top'] = rv_combined['q99']
    rv_combined['_bar_top'] = rv_combined['q75']
    rv_combined['_bar_median'] = rv_combined['median']
    rv_combined['_bar_bottom'] = rv_combined['q25']
    rv_combined['_segment_bottom'] = rv_combined['q01']

    #manipulation to get the ["X"] columns as the factors 
    if not(emptyDF):
        rv_combined["x"] = rv_combined[[new_meta2,new_meta]].apply(tuple, axis=1)
        rv_combined.sort_values(by = "x",key=lambda col: col.str[0],inplace = True)
        rv_combined.sort_values(by = "x",key=lambda col: col.str[1],inplace = True)

        rv_mouseid["x"] = rv_mouseid[[new_meta2,new_meta]].apply(tuple, axis=1)
        rv_mouseid.sort_values(by = "x",key=lambda col: col.str[0],inplace = True)
        rv_mouseid.sort_values(by = "x",key=lambda col: col.str[1],inplace = True)

    #merge on the just created ["X"] column
        final_rv = pd.merge(rv_combined,rv_mouseid,on="x")
        #sort..
        final_rv.sort_values(by = "x",key=lambda col: col.str[0],inplace = True)
        final_rv.sort_values(by = "x",key=lambda col: col.str[1],inplace = True)

        final_rv = final_rv.rename(columns={"x": "cat_value"})
        
        #sort the column for colors to be displayed..
        final_rv.sort_values(by = f'{new_meta}_y',inplace = True)
    else:
        final_rv = pd.merge(rv_combined,rv_mouseid,on=new_meta)
        final_rv = final_rv.rename(columns={new_meta: "cat_value"})

    return final_rv 


##might be better to merge these two functions into one:
#get_colors_of_obs()
#get_order_of_obs()
def get_colors_of_obs(dsid: str, meta: str):
    final_dict = {}
    datadir = util.get_datadir("h5ad")
    for yamlfile in datadir.glob("*.yaml"):
        basename = yamlfile.name.replace(".yaml", "")
        if dsid == basename:
            with open(yamlfile, "r") as F:
                y = yaml.load(F, Loader=yaml.SafeLoader)
                if y["obs_meta"][meta].get("values"):
                    for key,data in y["obs_meta"][meta]["values"].items():
                        name = key
                        color = data.get("color")
                        #default is grey
                        final_dict[name] = "grey" if color == None else color
    return final_dict

def get_order_of_obs(dsid: str, meta: str):
    final_dict = {}
    datadir = util.get_datadir("h5ad")
    for yamlfile in datadir.glob("*.yaml"):
        basename = yamlfile.name.replace(".yaml", "")
        if dsid == basename:
            with open(yamlfile, "r") as F:
                y = yaml.load(F, Loader=yaml.SafeLoader)
                if y["obs_meta"][meta].get("values"):
                    for key,data in y["obs_meta"][meta]["values"].items():
                        name = key
                        order = data.get("order")
                        #there is no order in the yaml..
                        if not(order):
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

    #in new yamls, there is sometimes names of column same as that of gene
    #e.g. in man2m.. meta = APOE, gene = APOE
    #this ruins the concatination..
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
    datadir = util.get_datadir("h5ad")
    try:
        rv = pl.read_parquet(datadir / f"{dsid}.X.prq", [gene])
    except pl.exceptions.SchemaError:
        return None
    return rv


def get_defields(dsid, view_name):
    ds = get_dataset(dsid, view_name)
    dex = ds.get("diffexp")
    return list(dex.keys())


def get_dedata(dsid, categ, genes, view_name: str = ""):
    """Return diffexp data."""
    ds = get_dataset(dsid, view_name)
    dex = ds.get("diffexp")
    assert categ in dex

    if isinstance(genes, str):
        genes = [genes]

    datadir = util.get_datadir("h5ad")
    rv = pl.read_parquet(datadir / f"{dsid}.var.prq", ["field"] + genes)
    rv = rv.to_pandas()
    rvx = rv["field"].str.split("__", expand=True)
    rvx.columns = ["categ", "cat_value", "measurement"]
    rv = pd.concat([rv, rvx], axis=1)
    rv = rv[rv["categ"] == categ].copy()
    del rv["categ"]
    del rv["field"]

    rv = rv.pivot(index="cat_value", columns="measurement", values=genes)
    return rv

# TODO next two functions can  be merged into 1?
def get_dedata_new(dsid, categ):
    datadir = util.get_datadir("h5ad")

    #to get 'gene' column
    last_col = len(pl.read_parquet(datadir / f"{dsid}.var.prq").columns)

    rv1 = pl.read_parquet(datadir / f"{dsid}.var.prq", [last_col - 1])

    # get logfoldchange and padjusted for one category
    # (example injection__None)
    rv2 = pl.read_parquet(
        datadir / f"{dsid}.var.prq", [categ + "__lfc"] + [categ + "__padj"])

    rv = pl.concat([rv2, rv1], how="horizontal")

    rv = rv.to_pandas()
    # in new format, padjusted and lfc are already on columns, genes on rows.
    # get genes has same index as that for the rows.
    return rv


def get_dedata_quadrant(dsid, categ1,categ2):
    datadir = util.get_datadir("h5ad")

    #to get 'gene' column
    last_col = len(pl.read_parquet(datadir / f"{dsid}.var.prq").columns)

    rv1 = pl.read_parquet(datadir / f"{dsid}.var.prq", [last_col - 1])

    # get logfoldchange and padjusted for one category
    # (example injection__None)
    rv2 = pl.read_parquet(
        datadir / f"{dsid}.var.prq", [categ1 + "__lfc"] + [categ2 + "__lfc"] + [categ1 + "__padj"] + [categ2 + "__padj"])

    rv = pl.concat([rv2, rv1], how="horizontal")

    #avoid where chosen lfc's are the same ones.. 
    rv = rv.rename({categ1 + "__lfc": categ1 + "__lfc_1",
                   categ2 + "__lfc": categ2 + "__lfc_2",
                   categ1 + "__padj": categ1 + "__padj_1",
                   categ2 + "__padj": categ2 + "__padj_2"})
    rv = rv.to_pandas()
    # in new format, padjusted and lfc are already on columns, genes on rows.
    # get genes has same index as that for the rows.
    return rv

def get_dedata_abundance(dsid,categ):
    datadir = util.get_datadir("h5ad")

    #to get 'gene' column
    last_col = len(pl.read_parquet(datadir / f"{dsid}.var.prq").columns)

    rv1 = pl.read_parquet(datadir / f"{dsid}.var.prq", [last_col - 1])

    rv2 = pl.read_parquet(
        datadir / f"{dsid}.var.prq", [categ + "__lfc"] + [categ + "__padj"] + [categ + "__lcpm"])
    
    rv = pl.concat([rv2, rv1], how="horizontal")
    rv = rv.to_pandas()

    return rv




def get_meta(dsid, col, raw=False, nobins=8, view_name: str = ""):
    """Return one obs column."""
    datadir = util.get_datadir("h5ad")
    rv = pl.read_parquet(datadir / f"{dsid}.obs.prq", [col])

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
    datadir = util.get_datadir("h5ad")
    lg.info("getting genes from " + str(datadir / f"{dsid}.X.prq"))
    X = pl.scan_parquet(datadir / f"{dsid}.X.prq")
    return X.columns


# @diskcache()
def obslist(dsid):
    lg.warning("`obslist` is deprecated function, use get_obsfields")
    return get_obsfields(dsid)


def get_obsfields(dsid):
    """Return a list of obs columns for this datset."""
    datadir = util.get_datadir("h5ad")
    X = pl.scan_parquet(datadir / f"{dsid}.obs.prq")
    return X.columns


def get_varfields(dsid):
    """Return a list of var columns for this datset."""
    datadir = util.get_datadir("h5ad")
    X = pl.scan_parquet(datadir / f"{dsid}.var.prq")
    return X.columns
