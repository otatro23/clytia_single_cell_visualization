# clytia_single_cell_visualization
Simple script to visualize expression of genes on interest in clytia single cell RNA sequencing data

1. Download medusa data from Chari et al. 2021 (DOI: 10.1126/sciadv.abh1683) and/or planula data from Ramon-mateu et al. 2021 (DOI: 10.1126/sciadv.adv1159). Look for "Annotated Clustered, Kallisto-processed Clytia Starvation Data" for medusa dataset and "scanpy_object_atlas.zip" for planula dataset.
2. Download visualizations.py and planula_cell_types.csv from this github repo.
3. Use visualizations.py to visualize expression of genes of interest on umaps and violin plots. Note that the h5ad files need to be named "medusa.h5ad" and "planula.h5ad".

Example call to visualizations.py:
./visualizations.py medusa example_genes.txt
