#!/mnt/home/plachetzki/omt1027/.conda/envs/sc_env/bin/python

import anndata as ad
import argparse
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import os
import pandas as pd
import scanpy as sc
from shiny.express import input, render, ui

def get_cell_annotations(adata):
    df = pd.read_csv("planula_cell_types.csv", index_col=0)

    # convert df index to string to match the leiden_new column of adata.obs, then merge
    df.index = df.index.astype(str)
    adata.obs = adata.obs.merge(df, left_on = "leiden_new", right_index = True, how = "left")

    return adata

def get_cell_resolution(cell_type, dataset):

    if dataset == "planula":
        if cell_type == "broad cell type":
            cell_res = "broad_cell_type"
        else:
            cell_res = "cell_type"
    else:
        if cell_type == "broad cell type":
            cell_res = "annos"
        else:
            cell_res = "annosSub"

    return cell_res

def read_adata(dataset):
    # read in adata and initialize cell type column names and figure output directory for medusa
    if dataset == "medusa":
        adata = sc.read_h5ad("medusa.h5ad")

    # read in adata, initialize cell type column names and figure output directory, and read in cell type annotations table for planula
    if dataset == "planula":
        if not os.path.isfile("planula_processed.h5ad"):
            adata = sc.read_h5ad("planula.h5ad")

            # merge the cell type annotations table into the adata
            adata = get_cell_annotations(adata)
            adata.write_h5ad("planula_processed.h5ad")
        else:
            adata = sc.read_h5ad("planula_processed.h5ad")

    return adata


ui.page_opts(title="Visualize gene expression in Clytia hemisphaerica scRNA-seq data")

# set up sidebar with options for cell type resolution, dataset, and gene of interest
with ui.sidebar():
    ui.input_select("dataset", "Select dataset", choices=["medusa", "planula"])
    ui.input_select("cell_type", "Select cell type resolution", choices=["broad cell type", "cell subtype"])
    ui.input_select("plot_type", "Select plot type", choices=["UMAP", "violin"])
    ui.input_text("gene", "Gene of interest:", "XLOC_045850")  

@render.plot
def umap():
    adata = read_adata(input.dataset())
    cell_res = get_cell_resolution(input.cell_type(), input.dataset()) # returns string for the cell type column name in the adata that matches the specified resolution

    if input.plot_type() == "UMAP":
        fig, axes = plt.subplots(1, 2, figsize=(14, 5))

        sc.pl.umap(adata, color=input.gene(), ax=axes[0], show=False)
        sc.pl.umap(adata, color=cell_res, ax=axes[1], show=False, legend_loc="right margin")

        plt.tight_layout()
        return fig
    
    if input.plot_type() == "violin":
    
        fig, axes = plt.subplots(figsize=(14, 5))
        sc.pl.violin(adata, keys=input.gene(), groupby=cell_res, rotation=30, ax=axes, show=False)

        for label in axes.get_xticklabels():
            label.set_ha("right")

        plt.tight_layout()
        return fig

