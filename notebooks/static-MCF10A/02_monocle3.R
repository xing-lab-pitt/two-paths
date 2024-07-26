#LOADING AND INSTALLING PACKAGES
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.14")

BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor', 'Matrix.utils'))
install.packages("devtools")
devtools::install_github('cole-trapnell-lab/leidenbase')
devtools::install_github('cole-trapnell-lab/monocle3')

devtools::install_github(repo = "hhoeflin/hdf5r")
devtools::install_github(repo = "mojaveazure/loomR", ref = "develop")

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[2] = "out.txt"
}

#LOAD PACKAGES
library(monocle3)
library(anndata)
library(loomR)
library(readr)
library(Matrix)

#exp_mat <- read_csv("dose_12_exp_mat.csv", col_names=FALSE)
exp_mat <- read_csv(args[1], col_names=FALSE)
exp_mat <- t(exp_mat)
row.names(exp_mat) <- row.names(gene_metadata)

#cell_metadata <- read.table("dose_12_cell_metadata.csv", header=TRUE, row.names=1, sep=",")
#gene_metadata <- read.table("dose_12_gene_metadata.csv", header=TRUE, sep=",", row.names=1)
cell_metadata <- read.table(args[2], header=TRUE, row.names=1, sep=",")
gene_metadata <- read.table(args[3], header=TRUE, sep=",", row.names=1)

exp_dense_matrix <- Matrix(as.matrix(exp_mat), sparse = TRUE)
exp_dense_matrix <- t(exp_dense_matrix)
cds <- new_cell_data_set(exp_dense_matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_metadata)

cds <- preprocess_cds(cds, num_dim = 100)
plot_pc_variance_explained(cds)

cds <- reduce_dimension(cds)
plot_cells(cds)

plot_cells(cds, genes=c('FN1', 'VIM','TOP2A','SNAI1'))

cds = cluster_cells(cds, resolution=1e-5, reduction_method = "UMAP")
cds <- learn_graph(cds)
plot_cells(cds,
           genes=c('FN1'),
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           cell_size=0.9)

plot_cells(cds,
           #color_cells_by = "embryo.time.bin",
           label_cell_groups=FALSE,
           label_leaves=TRUE,
           label_branch_points=TRUE,
           graph_label_size=1.5,
           cell_size=0.75)

cds <- order_cells(cds)

plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=TRUE,
           label_branch_points=TRUE,
           graph_label_size=1.5,
           trajectory_graph_segment_size = 2,
           cell_size=1.5)
