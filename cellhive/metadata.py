
import logging
from typing import Any, Dict, Optional

import anndata
import pandas as pd
import scanpy as sc
import scipy
from anndata import AnnData

from .util import query_pubmed

lg = logging.getLogger(__name__)


mandatory_metadata = dict([
    a.strip().split(None, 1) for a in
    """
    author          who's responsible?
    title           Title for this dataset
    organism        Organism this data dervies from
    study           short study identifier to group experiments
    experiment      short experiment identifiers unique for this dataset
    """.strip().split("\n")])


def prepare(adata: AnnData) -> None:
    """Check h5ad."""

    if 'cellhive' not in adata.uns:
        adata.uns['cellhive'] = {}

    for sub in ['layers', 'metadata', 'obs', 'obsm']:
        if sub not in adata.uns['cellhive']:
            adata.uns['cellhive'][sub] = {}

    for old in ['obs_keep', 'behave_obs_meta', 'obs_skip',
                'termite']:
        if old in adata.uns:
            del adata.uns[old]


def check(adata: AnnData) -> None:
    """
    Checks if metadata fields are defined in AnnData object.

    Parameters:
    - adata (AnnData): Annotated data object.

    Returns:
    - None

    Raises:
    - None
    """
    # Get the keys of the metadata in the AnnData object
    mdkeys = adata.uns['cellhive']['metadata'].keys()

    msg = []

    # Check if each mandatory metadata field is defined
    for k, v in mandatory_metadata.items():
        if k not in mdkeys:
            # Print a warning message if the field is not defined
            msg.append(f"Warning: Metadata field '{k}' not defined: {v}")

    adata.uns['cellhive']['metadata']['version'] = \
        str(adata.uns['cellhive']['metadata']['version'])

    if len(adata.obsm_keys()) == 0:
        msg.append("No UMAP/TSNE data (obsm is empty)")

    if len(msg) == 0:
        print("All seems fine")
    else:
        for m in msg:
            print(m)


def set1(adata: AnnData, key: str, value: Any) -> None:
    """
    Set one key/value in the cellhive data structure.

    Parameters:
    - adata (AnnData): Annotated data object.
    - key (str): Key to set in the metadata field.
    - value (Any): Value to assign to the specified key.

    Returns:
    - None

    Raises:
    - None
    """
    # Ensure the cellhive data structure is prepared
    prepare(adata)

    # Set the key/value in the metadata field of the cellhive data structure
    adata.uns['cellhive']['metadata'][key] = value


def md(adata, **kwargs):
    """
    Show or set metadata.

    If you provide pubmed_id, metadata is read from pubmed.
    """

    if 'pubmed' in kwargs:
        pubmed_id = kwargs['pubmed']
        pmd = query_pubmed(pubmed_id)
        set1(adata, 'author', pmd['author'])
        set1(adata, 'title', pmd['title'])
        set1(adata, 'year', pmd['year'])
        set1(adata, 'abstract', pmd['abstract'])
        set1(adata, 'doi', pmd['doi'])

    for k, v in kwargs.items():
        set1(adata, k, v)

    if 'cellhive' not in adata.uns:
        print("no annotations")
        return

    return adata.uns['cellhive']['metadata']


def get_layerdata(adata):
    """Create stats on the layers in this adata"""
    ldata = adata.uns['cellhive']['layers']

    def check_layer(name: str, X):

        if name not in ldata:
            ldata[name] = {}

        if isinstance(X, anndata._core.sparse_dataset.SparseDataset):
            X = X.value.todense()
        elif isinstance(X, scipy.sparse._csr.csr_matrix):
            X = X.todense()

        row: Dict[str, Any] = {}
        row['name'] = name
        row['dtype'] = X.dtype
        row['rows'] = X.shape[0]
        row['columns'] = X.shape[1]
        row['entries'] = entries = X.shape[0] * X.shape[1]
        row['min'] = X.min()
        row['max'] = X.max()

        zeros = len((X==0).nonzero()[0])
        perc_zeros = 100 * zeros / entries

        row['no_zeros'] = zeros
        row['% zeros'] = perc_zeros
        row['ignore?'] = ldata.get(name, {}).get('ignore', '-')
        row['description'] = ldata.get(name, {}).get('description', '-')

        chtype = ldata.get(name, {}).get('type', '-')
        if chtype == '-':
            looks_like_int = (row['max'] == int(row['max']))
            if looks_like_int and row['max'] > 500:
                # I think this may be raw counts
                lg.warning(f"Setting layer {name} to type count!")
                chtype = ldata[name]['type'] = 'count'
            elif not looks_like_int and row['max'] > 500:
                # I think this may be rpm counts
                lg.warning(f"Setting layer {name} to type: rpm!")
                chtype = ldata[name]['type'] = 'rpm'
            elif not looks_like_int and row['max'] < 50:
                # I think this may be rpm counts
                lg.warning(f"Setting layer {name} to type: logrpm!")
                chtype = ldata[name]['type'] = 'logrpm'
        row['data type'] = chtype
        return row

    layerdata = []
    layerdata.append(check_layer('X', adata.X))
    for layer in adata.layers.keys():
        layerdata.append(check_layer(layer, adata.layers[layer]))

    return pd.DataFrame(layerdata)


def layers(adata,
           name: str | None = None,
           ltype: str | None = None,
           ignore: bool | None = None,
           description: str | None = None) -> Optional[pd.DataFrame]:

    """Get or set layer data.

    ltype
    Usage:
       layers(): returns a dataframe
       layers(name='layername', [type='layertype', [load=True/False]]):
    """
    if name is not None:
        ldata = adata.uns['cellhive']['layers']
        if name not in ldata:
            ldata[name] = {}
        if ltype is not None:
            ltype = ltype.lower()
            if ltype not in ['count', 'rpm', 'logrpm', 'cell_abundance']:
                lg.warning(
                    f"Layer type for {name} is not count, rpm, logrpm, or cell_abundance!")
            ldata[name]['type'] = ltype
        if ignore is not None:
            ldata[name]['ignore'] = ignore
        if description is not None:
            ldata[name]['description'] = description

    df = get_layerdata(adata)
    return df.T


def obsm(adata: sc.AnnData,
         name: str | None = None,
         ignore: bool | None = None,
         description: str | None = None) -> pd.DataFrame | None:
    """
    Retrieve and process dimensionality reduction matrici data from `adata`.

    Args:
        adata (sc.AnnData): Anndata object containing the data.
        name (str | None, optional): Name of the `obsm` data to retrieve. Defaults to None.
        ignore (bool | None, optional): Boolean value indicating whether to ignore the `obsm` data. Defaults to None.

    Returns:
        pd.DataFrame: Transposed DataFrame containing the retrieved `obsm` data.

    Raises:
        AssertionError: If `name` is not in the keys of `adata.obsm`.

    """

    prepare(adata)
    ao = adata.uns['cellhive']['obsm']

    if name is not None:
        assert name in adata.obsm_keys()
        if name not in ao:
            ao[name] = {}
        if ignore:
            ao[name]['ignore'] = ignore
        if description is not None:
            ao[name]['description'] = description

    rv = []

    if len(adata.obsm_keys()) == 0:
        lg.warning("No dimred data in obsm?")
        return None

    for o in adata.obsm_keys():
        if o.startswith('_'):
            continue

        rv.append(dict(
            name=o,
            dim=adata.obsm[o].shape[1],
            ignore=ao.get(o, {}).get('ignore', '-'),
            description=ao.get(o, {}).get('description', '-'),
        ))

    return pd.DataFrame(rv).set_index('name').T


def obs_conv_int(adata, column):
    """Convert column to integers"""
    adata.obs[column] = adata.obs[column].astype(int)

def obs_conv_float(adata, column):
    """Convert column to float"""
    adata.obs[column] = adata.obs[column].astype(float)

def obs_conv_categorical(adata, column):
    """Convert column to categorical"""
    adata.obs[column] = adata.obs[column].astype(str)


OBS_CONV_FUNCTIONS = {
    'int': obs_conv_int,
    'float': obs_conv_float,
    'cat': obs_conv_categorical
    }

def get_obs_col_md(adata: AnnData, col: str, key: str,
                   default: Any = None) -> Any:
    """Get a metadata field on an obs column."""

    return adata.uns['cellhive']['obs']\
        .get(col, {})\
        .get(key, default)


def set_obs_col_md(adata: AnnData, col: str, key: str, val: Any) -> None:
    """Set a metadata field on an obs column."""
    if not key in adata.uns['cellhive']['obs']:
        adata.uns['cellhive']['obs'][col] = {}
    adata.uns['cellhive']['obs'][col][key] = val


def get_obs_column_metadata(adata: AnnData,
                            column: str,
                            max_categories: int = 30):
    """
    Determine what datatye a column is
    """

    # Guess the format from the raw data
    obs_col = adata.obs[column]
    no_uniq = len(obs_col.unique())
    no_example = 3
    example = ", ".join(map(str, obs_col.head(no_example)))

    forced = False

    if get_obs_col_md(adata, column, 'dtype'):
        odtype = get_obs_col_md(adata, column, 'dtype')
        forced = True
    else:
        odtype = str(obs_col.dtype)


    if 'int' in odtype and len(obs_col.unique()) < 20:
        if ('leiden' in column.lower()) or ('leuven' in column.lower()):
            #these are likely cluster names
            odtype = 'cat'
            forced = True
            set_obs_col_md(adata, column, 'dtype', 'cat')
            lg.warning("Expect {column} to be cluster IDs, forcing to categorical")


    print("??", odtype)
    if odtype in ['categorical', 'category', 'object', 'bool', 'str']:
        # guess categorical
        uniq = len(obs_col.unique())
        example = ", ".join(map(str, obs_col.value_counts().sort_values()[:no_example].index))
        dtype = 'cat'
        # unless....
        if not forced:
            if uniq > max_categories:
                # too many categories
                dtype = 'skip'
            elif uniq == 1:
                # not very interesting either...
                dtype = 'skip'

    elif 'int' in odtype:
        dtype = 'int'
    elif odtype.startswith('float'):
        dtype = 'float'
        example=", ".join(map(lambda x: f"{x:.3g}",
                             obs_col.head(no_example)))

    elif odtype == 'skip':
        dtype = 'skip'
    else:
        lg.warning(f"Unknown data type for {column} - {obs_col.dtype}")
        dtype = 'skip'

    return dict(name=column,
                dtype=dtype,
                no_uniq=no_uniq,
                example=example)


def prepare_obs(adata: AnnData) -> pd.DataFrame:

    md_obscol = {}
    td = adata.uns['cellhive']

    # regular obs columns
    lg.info("Processing obs table")

    for column in adata.obs.keys():
        if column.startswith('_'):
            continue
        md_obscol[column] \
            = get_obs_column_metadata(adata, column)

        md_obscol[column]['description'] = \
            get_obs_col_md(adata, column, 'description', '-')

    return md_obscol


def obs(adata: AnnData,
        column: str | None = None,
        dtype: str | None = None,
        ignore: bool | None = None,
        description: str | None = None,
        select: str = None) \
        -> pd.DataFrame | pd.Series | None:

    print('-====')
    prepare(adata)

    if column is not None:
        if dtype is not None:
            if dtype in ['str', 'categorical']:
                dtype = 'cat'

            if dtype not in ['int', 'float', 'cat']:
                lg.error("dtype must be one of: int, float or cat")
                return None
            else:
                lg.info(f"Convert obs column {column} to {dtype}")
                OBS_CONV_FUNCTIONS[dtype](adata, column)
                set_obs_col_md(adata, column, 'dtype', dtype)

        if description is not None:
            set_obs_col_md(adata, column, 'description', description)

    obs_data = pd.DataFrame(prepare_obs(adata)).T
    obs_data = obs_data.reset_index(drop=True)

    if column is not None:
        return obs_data[obs_data['name'] == column].T
    else:
        if select in ['str', 'categorical', 'cat']:
            obs_data = obs_data[obs_data['dtype'] == 'categorical']
        elif select == 'int':
            obs_data = obs_data[obs_data['dtype'] == 'int']
        elif select == 'float':
            obs_data = obs_data[obs_data['dtype'] == 'float']
        elif select in ['num', 'numerical']:
            obs_data = obs_data[obs_data['dtype'].isin(['int', 'float'])]

        return obs_data



# def todb(logrpm: bool = True,
#          skip_layers: bool = False,
#          skip_obs: bool = False,
#          skip_obsm: bool = False,
#          dim2load: int=2 ) -> None:

#     adata = globals.adata
#     check(adata)

#     # remove obs_names - enforce numbers to reduce database size
#     expname = adata.uns['cellhive']['metadata']['experiment']

#     exp_id = cellhive.db.autoincrementor('dataset_md', 'experiment', expname)
#     lg.info(f"Storing experiment {expname} with id {exp_id}")

#     if not skip_obs:
#         cellhive.h5ad.import_experiment_obs(exp_id, adata)

#     if not skip_obsm:
#         for (o, od) in adata.uns['cellhive']['obsm'].items():
#             if not od.get('load'):
#                 continue
#             lg.info(f"loading obsm: {o}")
#             # default - load first 2 dimensions
#             obsm = adata.obsm[o]
#             for i in range(min(dim2load, obsm.shape[1])):
#                 c = pd.DataFrame(
#                     pd.Series(adata.obsm[o][:,i],
#                               index=adata.obs_names))
#                 colname = f"{o}/{i:02d}"
#                 cellhive.h5ad.store_one_obs_col(
#                     exp_id=exp_id,
#                     colname=colname,
#                     col=c,
#                     original_name='-',
#                     dtype='dimred',
#                     dimred_name=o,
#                     dimred_dim=i)


#     layerdata = adata.uns['cellhive']['layers']

#     logrpm_in_adata = False
#     for lname, ldata in layerdata.items():
#         if ldata.get('type') == 'logrpm':
#             logrpm_in_adata = True

#     for lname, ldata in layerdata.items():
#         if not ldata.get('load'):
#             lg.info(f"Skipping load of layer {lname}")
#         ltype = ldata['type']
#         if not ltype:
#             lg.error(f"cannot load layer {lname} "
#                      + "with no type annotated")
#             exit(-1)

#         dataset_id = cellhive.h5ad.import_experiment_md(
#             adata, lname, ltype)

#         if not skip_layers:
#             cellhive.h5ad.import_counts(
#                 dataset_id, adata, lname)


#         if ltype == 'raw' and (not logrpm_in_adata) and logrpm:
#             dataset_id_2 = cellhive.h5ad.import_experiment_md(
#                 adata, 'auto_logrpm', 'logrpm')
#             if not skip_layers:
#                 cellhive.h5ad.import_counts(
#                     dataset_id_2, adata, lname, normalize='logrpm')

