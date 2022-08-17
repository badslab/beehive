from ast import Pass
from functools import partial, lru_cache
import logging
from typing import Dict, List, Tuple


import pandas as pd
import numpy as np
import polars as pl
import yaml

from beehive import util


lg = logging.getLogger(__name__)
lg.setLevel(logging.DEBUG)


DATASETS: Dict[str, Dict] = {}

diskcache = partial(
    util.diskcache, where=util.get_datadir("cache"), refresh=True)


@lru_cache(1)
def get_datasets(has_de: bool = False):
    """Return a dict with all dataset."""
    datadir = util.get_datadir("h5ad")
    global DATASETS

    if len(DATASETS) == 0:
        for yamlfile in datadir.glob("*.yaml"):
            basename = yamlfile.name
            basename = basename.replace(".yaml", "")
            with open(yamlfile, "r") as F:
                y = yaml.load(F, Loader=yaml.SafeLoader)
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


def get_dataset_siblings(dsid: str) -> dict:
    """Return datasets with the same study_id."""
    dsets = get_datasets()
    dset = dsets[dsid]
    study_id = dset["study"]
    siblings = {}
    for d, dd in dsets.items():
        if dd["study"] == study_id:
            siblings[d] = dd
    return siblings


@lru_cache(32)
def get_dataset(dsid):
    """Return metadata on a single dataset."""
    rv = get_datasets()[dsid]
    lg.info(f"Returning dataset {dsid}")
    return rv


def get_facet_options(dsid: str,
                      only_categorical: bool = False, include_skip = False) -> List[Tuple[str, str]]:
    """
    Return obs columns that a dataset can be facetted on.
    These have to be categorical

    param dataset - dictionary as loaded from the .yaml files
    """
    dataset = get_dataset(dsid)
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


def get_facet_groups(dsid,facet):
    """
    Return metadata of obs column
    """
    dataset = get_dataset(dsid)
    result = {}
    for key, data in dataset.get('obs_meta', {}).items():
        if key == facet:
            for group in data.get("values",{}).items():
                name = str(group[1].get("name"))
                order = group[1].get("order") #will default to none if not existing
                color = group[1].get("color") #will default to none if not existing
                result[name] = {'order': order,'color': color}

    return result

def get_legend_of_obs(dsid,meta):
    datadir = util.get_datadir("h5ad")
    legend = ""
    for yamlfile in datadir.glob("*.yaml"):
        basename = yamlfile.name.replace(".yaml", "")
        if dsid == basename:
            with open(yamlfile, "r") as F:
                y = yaml.load(F, Loader=yaml.SafeLoader)
                legend = y["obs_meta"][meta]["legend"]
    return legend

def get_facet_options_numerical(dsid: str) -> List[Tuple[str, str]]:
    """
    Return obs columns that a dataset can be facetted on.
    These have to be categorical

    param dataset - dictionary as loaded from the .yaml files
    """
    dataset = get_dataset(dsid)
    rv = []
    for key, data in dataset.get('obs_meta', {}).items():
        use = False
        if data.get('dtype') == 'numerical':
            use = True
        if use:
            name = data.get('name', key)
            rv.append((key, name))
    return list(sorted(rv))

def get_gene_meta_three_facets(dsid: str, gene: str, meta1: str, meta2: str, meta3: str = "mouse.id",nobins:int = 8):
    genedata = get_gene(dsid,gene)
    metadata = get_meta(dsid, meta1, nobins=nobins, raw = True) #facet1
    metadata2 = get_meta(dsid,meta2,nobins=nobins, raw = True) #facet2
    metadata3 = get_meta(dsid,meta3,nobins=nobins,raw = True) #most likely mouse_id
    """"
    function to return the aggregate of 3 metadata options
    first two metadatas will be displayed in a boxplot doubled
    and the third will be specific to every possible combination
    we end up with means of the third metadata
    and means of (metadata1 * metadata2) groups.
    """

    if genedata is None or metadata is None or metadata2 is None or metadata3 is None:
        return None

    if genedata.columns == metadata.columns:
        metadata.columns = [metadata.columns[0] + "_category"]
    if genedata.columns == metadata2.columns:
        metadata2.columns = [metadata2.columns[0] + "_category"]
    if genedata.columns == metadata3.columns:
        metadata3.columns = [metadata3.columns[0] + "_category"]

    new_meta = metadata.columns[0]
    new_meta2 = metadata2.columns[0]
    new_meta3 = metadata3.columns[0]
    
    #group on both metadatas => excluding the third
    rv_combined = (
        pl.concat([genedata, metadata,metadata2,metadata3], how="horizontal")
        .groupby([new_meta,new_meta2])
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

    #group on all => get the means of the third.
    rv_mouseid = (
        pl.concat([genedata, metadata,metadata2,metadata3], how="horizontal")
        .groupby([new_meta,new_meta2,new_meta3])
        .agg(
            [
                pl.count().alias("count_" + meta3), #count_mouseid
                pl.mean(gene).alias("mean_" + meta3), #mean_mouseid
            ]
        )
    )

    #switch to dfs
    rv_combined = rv_combined.to_pandas()
    rv_mouseid = rv_mouseid.to_pandas()

    #need to exclude NONE here. (when yaml is updated...)
    #waiting on reply from nico
    #what to do when all intersected groups have NONE?
    #resulting dataframe will not include any data at all.
    # rv_combined = rv_combined[rv_combined[new_meta] != 'NONE']
    # rv_combined = rv_combined[rv_combined[new_meta2] != 'NONE']

    # rv_mouseid = rv_mouseid[rv_mouseid[new_meta3] != 'NONE']

    #calculate other measures
    rv_combined['perc'] = 100 * rv_combined['count'] / rv_combined['count'].sum()
    rv_combined['_segment_top'] = rv_combined['q99']
    rv_combined['_bar_top'] = rv_combined['q75']
    rv_combined['_bar_median'] = rv_combined['median']
    rv_combined['_bar_bottom'] = rv_combined['q25']
    rv_combined['_segment_bottom'] = rv_combined['q01']


    #manipulation to get the ["X"] columns as the factors 
    rv_combined["x"] = rv_combined[[new_meta2,new_meta]].apply(tuple, axis=1)
    rv_combined["x"] = sorted(rv_combined["x"],key=lambda tup: tup[0])
    rv_combined["x"] = sorted(rv_combined["x"],key=lambda tup: tup[1])

    rv_mouseid["x"] = rv_mouseid[[new_meta2,new_meta]].apply(tuple, axis=1)
    rv_mouseid["x"] = sorted(rv_mouseid["x"],key=lambda tup: tup[0])
    rv_mouseid["x"] = sorted(rv_mouseid["x"],key=lambda tup: tup[1])

    #merge on the just created ["X"] column
    final_rv = pd.merge(rv_combined,rv_mouseid,on="x")
    #sort..
    final_rv["x"] = sorted(final_rv["x"],key=lambda tup: tup[0])
    final_rv["x"] = sorted(final_rv["x"],key=lambda tup: tup[1])

    final_rv = final_rv.rename(columns={"x": "cat_value"})
    #sort the column for colors to be displayed..

    final_rv[new_meta+"_y"] = sorted(final_rv[new_meta+"_y"])

    #note, returning factors as we need it for the x_range of the plot.
    return final_rv #returns the different degrees of freedom = factors (meta1 x meta2) and the final rv 


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
                        final_dict[name] = "black" if color == None else color
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
                        final_dict[name] = 0 if order == None else order
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

    # Does the YAML have sort_order fields? If so use these
    # if False:  # if there are order fields:
    #     # TODO: Raghid
    #     pass
    # else:
    # Otherwise - do the following:
    try:
        # attempt a numerical sort on the cat_value
        # try to prevent alphanumerically sorted numeric
        # values 1 10 11 2 3

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

#NOT ADDED TO VIEWS YET
#Suggested format:
#Views
def get_usable_views(dsid):
    dataset = get_dataset(dsid)
    return


def get_defields(dsid):
    ds = get_dataset(dsid)
    dex = ds.get("diffexp")
    return list(dex.keys())


def get_dedata(dsid, categ, genes):
    """Return diffexp data."""
    ds = get_dataset(dsid)
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

    # TODO should be a stable name "gene" in all datasets
    # get last column,, no way to access it just by name
    # sometimes it's called gene, sometimes it's index__0__
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

    # TODO should be a stable name "gene" in all datasets
    # get last column,, no way to access it just by name
    # sometimes it's called gene, sometimes it's index__0__
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



def get_meta(dsid, col, raw=False, nobins=8):
    """Return one obs column."""
    datadir = util.get_datadir("h5ad")
    rv = pl.read_parquet(datadir / f"{dsid}.obs.prq", [col])

    if raw:
        # just return whatever is in the db.
        return rv

    ds = get_dataset(dsid)
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
