
'''
Single-cell Workshop
Day two
Chris Large

Using worm data to demonstrate basic data analysis
Pipeline:

1) Load data
2) Normalize the single-cell expression data
3) Annotate large groups of cells using markers
4) Annotate the ciliated neurons
5) Identify markers of the AFD ciliated neuron and the muscle cells
6) Which genes are enriched in those cell types?
'''

# First, let's set the path that you've downloaded the data to
data_path <- "/Path/to/single/cell/workshop/"

library(Seurat)
matrix <- readRDS(paste0(data_path, "workshop_subset_matrix.rds"))

# Create a seurat object from the matrix
# We adjust the minimum number of cells and features to not filter anything out
cel <- CreateSeuratObject(counts = matrix, project = "c.elegans",
                          min.cells = 0, min.features = 0)

# Run the standard Seurat pipeline
# For running the principal component analysis, use 50 dimensions
https://satijalab.org/seurat/articles/essential_commands

# Shows the cells in a reduced dimension plot (UMAP) with the clusters we identified
DimPlot(object = cel, reduction = "umap", label = TRUE)

# Here the clusters are labeled with numbers, but we want to label them with the tissue/cell type

# Our goal is to label the clusters with their tissue or cell identity
# For example, we can do labeling of the intestine cells together
# intestine specific gene markers expression pattern
FeaturePlot(object = cel, features = c("end-1", "elt-2", "elt-7", "pgrn-1"), sort.cell = TRUE)

# We then want to relabel the indets of the cluster to the name of the tissue/cell type
# For example, if cluster 0 was the intestine, we would rename it to "intestine"
cel <- RenameIdents(object = cel,
                    '0' = "intestine")

# Now let's identify the other clusters using the marker sets in the excel table
# By the end, we should have cluster labels for all of the tissue/cell types in the table


# Next let's find the genes that are enriched in the AFD neurons
AFD_markers <- FindMarkers(cel, ident.1 = "AFD", logfc.threshold = 0.25, only.pos = TRUE)

# Now do the same for the muscle cells


# Using the markers, we want to see if there is an enrichment of these genes in a certain cellular processes
# We need to first filter for genes that are enriched in the AFD neuron
AFD_makers_subet <- AFD_markers %>% # Pipes the markers object to the next line
  dplyr::filter(avg_log2FC > 1.5) %>% # Filters for genes with a log2FC > 1.5
  slice_head(n = 200) %>% # Takes the top 200 genes
  rownames() # Returns the gene names

# Now replicate this for the muscle cells


# Next let's use the functional annotations of the genes to see if there
# is an enrichment of certain processes in the AFD neurons and muscle cells
# The functional annotation of the genes is stored in the gene data object
# as the WormCategory columns. There are 3 'tiers' of annotations, going from 
# broad to specific. We will use the basic annotations for this analysis

# First we need to load the gene data object
gene_data <- read.csv(paste0(data_path, "workshop_gene_data.csv"))
rownames(gene_data) <- gene_data$gene_short_name

# Now we need to 'source' the code to run the enrichments
source(paste0(data_path, "workshop_source_code.R"))

# Finally, we can look at the enrichment of certain gene categories in our
# AFD and muscle cell markers
AFD_enrichment <- run_worm_cat(AFD_makers_subet, # Our subsetted markers
             gene_data,                          # The background gene data object
             1)                                  # WormCat tier (1-3)

AFD_enrichment

# Now let's look at the enrichment of the muscle cell markers


# How do the enriched categories differ between the AFD neurons and muscle cells?
