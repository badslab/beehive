"""helpers for the gene expression app."""
import logging
import os
from pathlib import Path
import typer
import yaml
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as colors
import markdown

app = typer.Typer()

lg = logging.getLogger(__name__)

"""Small app for h5ad wrangling:"""
"""Creates .prq files needed for the visualization application"""
"""Updates the yaml file and extends on it"""
"""Required:"""
"""scanpy.h5ad file with "obs" and "var" """
"""'experiment.yaml' file with metadata information."""

def read_yaml_helper(outbase,new):
    lg.info(f"Filename for IO: {outbase}")

    dir_path = os.path.dirname(outbase)
    if new:
        yamlfile = os.path.join(dir_path, "experiment_new.yaml")
    else:
        yamlfile = os.path.join(dir_path, "experiment.yaml")

    if not(os.path.exists(yamlfile)):
        print("No yaml file found to print metadata. Try again.")
        return None
    else:
        with open(yamlfile) as F:
            yml = yaml.load(F, Loader=yaml.SafeLoader)
            return yml


@app.command("metadata")
def h5ad_uns(h5ad_file: Path = typer.Argument(..., exists=True),
            new: bool = typer.Option(False, "--new")):
    """
    Prints out the metadata yaml from the provided h5ad file.
    Will either look for an expermient.yaml or experiment_new.yaml
    """
    import sys
    outbase = h5ad_file
    yml = read_yaml_helper(outbase,new)
    if yml:
        yaml.dump(yml, sys.stdout)
    return


def value_typer(key, val):
    """Convert values to the correct datatype - if necessary"""
    if key in ['year']:
        return int(val)
    return val

@app.command("set")
def h5ad_set(h5ad_file: Path = typer.Argument(..., exists=True),
             key: str = typer.Argument(...),
             val: str = typer.Argument(...),
             new: bool = typer.Option(False, "--new")):
    """
    Sets a new key value in the old/new yaml metadata file.
    Requires a parameter new, which is either True or False.
    If True,  use --new. It will work on the new experiment_new.yaml. If False, skip the parameter. It will work on the old experiment.yaml file.\n
    Note: the experiment_new.yaml file is automatically created once 'beehive h5ad prepare' is ran. \n 
        Example: beehive h5ad set /dir/to/file/scanpy.h5ad keytoremove keytoadd valtoadd --new \n
    This will add {'keytoadd': 'valtoadd'} to the experiment_new.yaml
"""
    # check if the output yaml is there as well
    outbase = h5ad_file
    yml = read_yaml_helper(outbase,new)
    if not(yml):
        lg.warning(f"No experiment.yaml or experiment_new.yaml found! Try again.")
        return
    ##set key value
    yml[key] = val
    
    #overwrite
    dir_path = os.path.dirname(outbase)
    if new:
        yamlfile = os.path.join(dir_path, "experiment_new.yaml")
    else:
        yamlfile = os.path.join(dir_path, "experiment.yaml")

    with open(yamlfile, 'w') as F:
        yaml.dump(yml, F, Dumper=yaml.SafeDumper)
    return


@app.command("del")
def h5ad_del(h5ad_file: Path = typer.Argument(..., exists=True),
             key: str = typer.Argument(...),
             new: bool = typer.Option(False, "--new")):
    """
    Deletes an existing key value in the old/new yaml metadata file.
    Requires a parameter new, which is either True or False.
    If True, it will work on the new experiment_new.yaml. If False, it will work on the old experiment.yaml file.\n
    Note: the experiment_new.yaml file is automatically created once 'beehive h5ad prepare' is ran. \n
    Example: beehive h5ad del /dir/to/file/scanpy.h5ad keytoremove --new \n
    This will remove 'keytoremove' from the experiment_new.yaml
    """

    outbase = h5ad_file
    yml = read_yaml_helper(outbase,new)
    if not(yml):
        lg.warning(f"No experiment.yaml or experiment_new.yaml found! Try again.")
        return
    ##delete key value
    try:
        del yml[key]
    except:
        pass

    #overwrite
    dir_path = os.path.dirname(outbase)
    if new:
        yamlfile = os.path.join(dir_path, "experiment_new.yaml")
    else:
        yamlfile = os.path.join(dir_path, "experiment.yaml")

    with open(yamlfile, 'w') as F:
        yaml.dump(yml, F, Dumper=yaml.SafeDumper)
    return


def string_cleanup(variable_str):
    if variable_str[0] == "_":
        variable_str = variable_str[1:]
    if variable_str[-1] == "_":
        variable_str = variable_str[:-1]
    return variable_str.title().replace("_"," ")
    

@ app.command("prepare")
def h5ad_convert(h5ad_file: Path = typer.Argument(..., exists=True),
                 num_categories: int = 20):
    """
    Convert to polars/parquet dataframes. 
    Provide a h5ad file according to the following:
    beehive h5ad prepare /example_dir/example.h5ad 
    Note: In the example_dir, you should also have a yaml file called: experiment.yaml 
    The new yaml generated will be called experiment_new.yaml. It will also generate a simple markdown file to be viewed on the website. \n
    Option --num-categories will be used to dictate which obs_meta and/or var_meta categorical
    variables to include. If a categorical variable has <= num_categories, then it will be included.
    Otherwise, it will not. For example, a 'full_id' variable will have a huge number of variables, and 
    therefore will be excluded. For number of clusters, if 10 clusters are expected, then num_categories 
    should be higher than that to include cluster categorical variable in the yaml. The default number is set to 20
    """
    
    import pandas as pd
    import polars as pl
    import scanpy as sc

    outbase = h5ad_file
    lg.info(f"Filename for IO: {outbase}")

    adata = sc.read_h5ad(h5ad_file)
    dir_path = os.path.dirname(outbase)

    # check if the output yaml is there as well
    yamlfile = os.path.join(dir_path, "experiment.yaml")
    assert os.path.exists(yamlfile), """yaml.experiment file is needed to complete the preparation of .prq files"""
    with open(yamlfile) as F:
        yml = yaml.load(F, Loader=yaml.SafeLoader)

    ### 1. Cleanup yaml with basic metadata such as author title...###
    yml.pop("cat_covariates")
    yml.pop("cont_covariates")
    yml["meta"] = {}
    yml["diffexp"] = {}
    yml["dimred"] = []

    keys_permanent = ["author","name","year","access","full_author_list","group_id","n_var_genes","normalised",
                        "organism","published_in","target_no_clusters","url","doi","version"]
    keys_permament_types = ["string","string","int","string","string","string","int","bool",
                            "string","string","int","string","string","int"]
    for ind,key in enumerate(keys_permanent):
        if key not in yml.keys():
            yml[key] = f'TBD ({keys_permament_types[ind]})'

    #### 2. Manipulation to scanpy object to extract var obs ####
    dfx = adata.to_df()
    obs = adata.obs
    var = adata.var
    var.index = adata.var_names #gene names.
    var.index.name = 'gene'

    #### 3. Add obs variables to yaml. Type will be determined from the value of the columns. ####
    #### Categorical variables have a special case. An order will be assigned, a color will be assigned ###
    #### For Categorical variables, only those that have categories less than num_categories will be included. ####
    #### For categorical variables in obs, diff exp __padj __lfc ___cell_frac will be calculated as well. ####

    var_variables = []
    obs_variables = []
    var_dict = {}
    obs_dict = {}
    for k, v in obs.iteritems():
        k_type = 'numerical'
        obs_variables = obs_variables + [k]
        variable_name = string_cleanup(k)
        if k == "full_id" or k == "subject":
            continue
        # polars/parquest does not like categories
        if str(v.dtype) == 'category':
            print(f'Found categorical variable: {k}.')

            list_cats = obs[k].cat.categories.tolist()
            if len(list_cats) <= num_categories:
                obs[k] = v.astype('str')
                k_type = 'categorical'
                cat_values = {}
                #alpha_numeric sort:
                #use this list in yaml.
                list_cats = sorted(list_cats)
                #colors
                colormap = plt.cm.get_cmap('viridis')
                rgba_colors = [colormap(i) for i in np.linspace(0, 1, len(list_cats))]
                hex_colors = [colors.rgb2hex(c) for c in rgba_colors]

            #apply test
                print(f'applying gene expression ranking method using Wilcoxon Test.')
                sc.tl.rank_genes_groups(adata,k,method="wilcoxon",pts = True)
            ##use this list to make var columns.
                cat_list_ranked_genes_groups = list(adata.uns["rank_genes_groups"]["logfoldchanges"].dtype.names)

                for ind,cat in enumerate(list_cats):
                    color = hex_colors[ind]
                    cat_name = string_cleanup(cat)
                    cat_values[cat] = {"name":cat_name, "color":color,"order":ind}
            

                for ind,cat in enumerate(cat_list_ranked_genes_groups):
                    ##Generate k__groupname__padj, k__groupname__lfc, k__groupname__pt
                    
                    cat_lfc = [x[ind] for x in adata.uns["rank_genes_groups"]["logfoldchanges"]]
                    cat_padj = [x[ind] for x in adata.uns["rank_genes_groups"]["pvals_adj"]]
                    cat_cell_frac = adata.uns["rank_genes_groups"]["pts"][cat].to_list()
                    gene_names = [x[ind] for x in adata.uns["rank_genes_groups"]["names"]]
                    data={f'{k}_{cat}__lfc': cat_lfc, f'{k}_{cat}__padj' : cat_padj, f'{k}_{cat}__cell_frac': cat_cell_frac}
                    df = pd.DataFrame(data, index = gene_names)
                    var = pd.merge(var,df,left_index=True, right_index=True)
                    print(f'Generated {k}_{cat}__lfc, {k}_{cat}__padj, {k}_{cat}__cell_frac for group {cat} against all others.')
                    print("********************")

                obs_dict[k] = {"name": variable_name, "dtype": k_type, "values": cat_values}
            else:
                print(f'Skipping categorical variable {k} since its categories ({len(list_cats)}) exceed {num_categories}.')
        else:
            obs_dict[k] = {"name": variable_name, "dtype": k_type,"legend": "TBD (Quoted Text)"}

    #### 4. Add var variables to yaml. Type will be determined from the value of the columns. ####
    #### Categorical variables have a special case. An order will be assigned, a color will be assigned ###
    #### For Categorical variables, only those that have categories less than num_categories will be included. ####
    for k,v in var.iteritems():
        var_variables = var_variables + [k]
        k_type = "numerical"
        variable_name = string_cleanup(k)

        if str(v.dtype) == "category":
            list_cats = var[k].cat.categories.tolist()
            if len(list_cats) <= num_categories:
                var[k] = v.astype('str')
                k_type = 'categorical'
                cat_values = {}
                #alpha_numeric sort:
                #use this list in yaml.
                list_cats = sorted(list_cats)
                #colors
                colormap = plt.cm.get_cmap('viridis')
                rgba_colors = [colormap(i) for i in np.linspace(0, 1, len(list_cats))]
                hex_colors = [colors.rgb2hex(c) for c in rgba_colors]
                for ind,cat in enumerate(list_cats):
                    color = hex_colors[ind]
                    cat_name = string_cleanup(cat)
                    cat_values[cat] = {"name":cat_name, "color":color,"order":ind}
                
                var_dict[k] = {"name": variable_name, "dtype": k_type, "values": cat_values, "legend": "TBD (Quoted Text)"}
            else:
                print(f'Skipping categorical variable {k} since its categories ({len(list_cats)}) exceed {num_categories}.')

            var[k] = v.astype('str')
        else:
            var_dict[k] = {"name": variable_name,"dtype":k_type, "legend": "TBD (Quoted Text)"}

    #### 5. Add the variables in obsms table ####
    obsms = []
    for k in adata.obsm_keys():
        if k.startswith('_'):
            continue
        yml['dimred'].append(k)

        oo = pd.DataFrame(adata.obsm[k], index=obs.index)
        oo.columns = '_' + k + '_' + oo.columns.astype(str)
        obsms.append(oo)

    #### 6. prepare final obs table, and final var table. ####
    obs = pd.concat([obs] + obsms, axis=1)
    obs.index.name = '_cell'
    obs = obs.reset_index()

    #transponse.
    var = var.T
    #rename index to 'field'
    var.index.name = 'field'
    #reset indices 0,1..
    var = var.reset_index()

    #### 7. prepare final yaml file ####
    yml["obs_meta"] = obs_dict
    yml["var_meta"] = var_dict


    #### 8. Dump into .prq files and make experiment_new.yaml file ####
    lg.info("Writing output files to:")
    lg.info(" - " + str(outbase.with_suffix('.obs.prq')))
    pl.DataFrame(obs).write_parquet(outbase.with_suffix('.obs.prq'))
    lg.info(" - " + str(outbase.with_suffix('.var.prq')))
    pl.DataFrame(var).write_parquet(outbase.with_suffix('.var.prq'))
    lg.info(" - " + str(outbase.with_suffix('.X.prq')))
    pl.DataFrame(dfx).write_parquet(outbase.with_suffix('.X.prq'))

    new_yaml_file = os.path.join(dir_path, "experiment_new.yaml")
    lg.info(" - " + "Writing new yaml file: " + new_yaml_file)
    with open(new_yaml_file, 'w') as F:
        yaml.dump(yml, F, Dumper=yaml.SafeDumper)
    
    #### 9. Make markdown file. ####
    content = f"""Title: {yml["name"]}'
Date: {yml["year"]}'
Category: Papers
Slug: {yml["author"]}{yml["year"]}'
Tags: TBD 
        
*{yml["full_author_list"]}*

Abstract: TBD

## Gene Expression Views:

* [gene_expression1](./gene_expression?dataset_id={yml["group_id"]}.{yml["version"]}')
* [volcano_plot1](./volcano_plot?dataset_id={yml["group_id"]}.{yml["version"]}')
* [scatter_plot1](./scatter_expression?dataset_id={yml["group_id"]}.{yml["version"]}')
    """
    new_yaml_file = os.path.join(dir_path, f'{yml["author"]}{yml["year"]}.md')
    with open(new_yaml_file, 'w') as f:
        f.write(content)
        
    return
