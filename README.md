# Single Cell Workshop
### Christopher Large
This repository will walk you through the basic annotation of _C. elegans_ single cell RNA-sequencing data.

To access the data from this repository, you will 'clone' the reposistory on your own computer.

First you need to make sure that you have git installed:

[Git Download](https://git-scm.com/downloads)

Then you will 'clone' this repository in a logical place on your computer with:
```
git clone git@github.com:livinlrg/single_cell_workshop.git
```

Inside of the downloaded repo, you will find:
- single_cell_workshop.R
  - The main set of code to follow
- workshop_subset_matrix.rds
  - The _C. elegans_ single-cell RNA-seq matrix
  - The matrix has been subsetted to just include terminal cell types
- workshop_source_code.R 
  - Some code to run gene set enrichments
- Cel_tissue_markers.xlsx
  - A list of markers for some terminal cell types in the animal
- workshop_gene_data.csv
  - Contains functional annotations of the genes in the dataset
