from functools import partial, lru_cache
import logging
from typing import Dict, List, Tuple


import pandas as pd
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
                      only_categorical: bool = False) -> List[Tuple[str, str]]:
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


def get_gene_meta_agg(dsid: str, gene: str, meta: str, nobins: int = 8):
    """Return gene and observation."""
    genedata = get_gene(dsid, gene)
    metadata = get_meta(dsid, meta, nobins=nobins)

    # def chksum(p: pl.DataFrame):
    #     import hashlib
    #     md5 = hashlib.md5()
    #     md5.update(p.to_csv().encode())
    #     return md5.hexdigest()

    if genedata is None:
        return None
    if metadata is None:
        return None

    rv = (
        pl.concat([genedata, metadata], how="horizontal")
        .groupby(meta)
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

    # print(dsid, gene, meta, nobins,
    #       chksum(genedata), chksum(metadata),
    #       chksum(rv))

    rv = rv.to_pandas()
    rv = rv.rename(columns={meta: "cat_value"})
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



def get_dedata_new(dsid,categ):
    datadir = util.get_datadir("h5ad")

    #get logfoldchange and padjusted for one category (example injection__None)
    rv = pl.read_parquet(datadir / f"{dsid}.var.prq",[categ+ "__lfc"] + [categ + "__padj"])

    rv = rv.to_pandas()
    #in new format, padjusted and lfc are already on columns, genes on rows.
    #get genes has same index as that for the rows.
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
