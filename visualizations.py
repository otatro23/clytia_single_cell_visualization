#!/mnt/home/plachetzki/omt1027/.conda/envs/sc_env/bin/python

import anndata as ad
import argparse
import matplotlib.pyplot as plt
import os
import pandas as pd
import scanpy as sc

# --dataset is a required argument (despite the flag) that can be either 'medusa' or 'planula' depending on the dataset being used. 
    # 'medusa.h5ad' or 'planula.h5ad' must be present in current working directory.
# genes is a path to a text file that contains the genes of interest, one gene per line
def get_args():
    parser = argparse.ArgumentParser(description = "Visualize the expression of a given gene in an h5ad object")
    parser.add_argument("--dataset", choices = ["planula", "medusa"], required = True, help = "specify using the 'planula' or 'medusa' single cell data. Files must be in current working directory names 'planula.h5ad' and 'medusa.h5ad'.")
    parser.add_argument("genes", help = "path to text file of genes of interest (one gene per line)")
    #parser.add_argument("fasta", help = "transcripts fasta. Header format: '>TCONS_00000000 XLOC_000022")
    args = parser.parse_args()
    return args

# not in use currently, but reads the transcripts fasta to map XLOC gene symbols to TCONS
def get_symbol_mapping(args):
    map_path = "symbol_map.tsv"
    rows = []

    if not os.path.isfile(map_path): 
        with open(args.fasta, "r") as fasta:
            for line in fasta:
                if not line.startswith(">"): # we are only interested in headers
                    continue
                fields = line.rstrip().split()
                tcons_id = fields[0][1:]
                xloc_id = fields[1]
                rows.append({"xloc" : xloc_id})

# reads genes in gene text file to a list
def get_genes(args):
    with open(args.genes, "r") as gene_file:
        genes = []
        for line in gene_file:
            genes.append(line.rstrip())

    return genes

# reads in cell annotation table for planula and uses the 'leiden_new' column to merge annotations onto the cell metadata (adata.obs)
def get_cell_annotations(adata):
    df = pd.read_csv("planula_cell_types.csv", index_col=0)

    # convert df index to string to match the leiden_new column of adata.obs, then merge
    df.index = df.index.astype(str)
    adata.obs = adata.obs.merge(df, left_on = "leiden_new", right_index = True, how = "left")

    return adata

# plots cell type umaps for each cell type column (broad and subtype) and each gene of interest
def plot_umaps(adata, genes, cell_types, output_dir):

    os.makedirs(output_dir + "/umaps/", exist_ok=True)

    for gene in genes:

        if gene not in adata.var_names:
            print(f"Gene {gene} was not found in adata.var_names")
            continue
        
        ax = sc.pl.umap(adata, color=gene, show=False)
        ax.figure.savefig(f"{output_dir}/umaps/{gene}_umap.png", dpi=500, bbox_inches="tight", facecolor="white")
        plt.close()

    for color in cell_types: # make UMAPS colored by cell types
        ax = sc.pl.umap(adata, color=color, show=False)
        ax.figure.savefig(f"{output_dir}/umaps/{color}_umap.png", dpi=500, bbox_inches="tight", facecolor="white")
        plt.close()

# makes violin plots to show the exression of the genes of interest in each cell type in a more concrete way than the Umaps
def plot_violins(adata, genes, cell_types, output_dir):
    os.makedirs(output_dir + "/violins/", exist_ok=True)

    for gene in genes:
    
        if gene not in adata.var_names:
            print(f"Gene {gene} was not found in adata.var_names")
            continue

        for cell_type in cell_types:
            fig, ax = plt.subplots(figsize=(10, 4))
            sc.pl.violin(adata, keys=[gene], groupby=cell_type, rotation=90, use_raw=False, ax=ax, show=False)

            plt.tight_layout()
            plt.savefig(f"{output_dir}/violins/{gene}_{cell_type}_violin.png", dpi=150, bbox_inches='tight')
            plt.show()

def main():
    args = get_args()

    # read in adata and initialize cell type column names and figure output directory for medusa
    if args.dataset == "medusa":
        try: 
            adata = sc.read_h5ad("medusa.h5ad")
        except FileNotFoundError:
            print("Ensure 'medusa.h5ad' is in current working directory.")
            quit()

        cell_types = ["annos", "annosSub"]
        figure_dir = "figures/medusa"

    # read in adata, initialize cell type column names and figure output directory, and read in cell type annotations table for planula
    if args.dataset == "planula":
        if not os.path.isfile("planula_processed.h5ad"):
            try:
                adata = sc.read_h5ad("planula.h5ad")
            except:
                print("Ensure 'planula.h5ad' is in current working directory.")
                quit()
            adata = get_cell_annotations(adata)
            adata.write_h5ad("planula_processed.h5ad")
        else:
            adata = sc.read_h5ad("planula_processed.h5ad")

        cell_types = ["cell_type", "broad_cell_type"]
        figure_dir = "figures/planula"
        
    # get genes of interest
    genes = get_genes(args)

    # plot figures 
    plot_umaps(adata, genes, cell_types, figure_dir)
    plot_violins(adata, genes, cell_types, figure_dir)

if __name__ == "__main__":
    main()