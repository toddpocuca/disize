library(utils)

base_path <- "benchmarks/rnaseq-data/data"
dir.create(base_path, recursive = TRUE, showWarnings = FALSE)

# Process 293T Dataset ---
# Download count matrix file
url <- "https://cf.10xgenomics.com/samples/cell-exp/1.1.0/293t/293t_filtered_gene_bc_matrices.tar.gz"
dest <- file.path(base_path, "293t_data.tar.gz")
download.file(url, destfile = dest, mode = "wb")

# Extract to temporary folder
temp_dir <- file.path(base_path, "temp_293t")
untar(dest, exdir = temp_dir)

# Move folder
from_path <- file.path(temp_dir, "filtered_matrices_mex", "hg19")
to_path <- file.path(base_path, "293t")
file.rename(from_path, to_path)

# Cleanup
unlink(temp_dir, recursive = TRUE)
unlink(dest)


# Process PBMC Dataset ---
# Download count matrix file
url <- "https://cf.10xgenomics.com/samples/cell-exp/4.0.0/SC3_v3_NextGem_SI_PBMC_10K/SC3_v3_NextGem_SI_PBMC_10K_filtered_feature_bc_matrix.tar.gz"
dest <- file.path(base_path, "pbmc_data.tar.gz")
download.file(url, destfile = dest, mode = "wb")

# Extract to temporary folder
temp_dir <- file.path(base_path, "temp_pbmc")
untar(dest, exdir = temp_dir)

# Path in tar is: filtered_feature_bc_matrix
from_path <- file.path(temp_dir, "filtered_feature_bc_matrix")
to_path <- file.path(base_path, "pbmc")
file.rename(from_path, to_path)

# Cleanup
unlink(temp_dir, recursive = TRUE)
unlink(dest)

# Download clustering analysis
url <- "https://cf.10xgenomics.com/samples/cell-exp/4.0.0/SC3_v3_NextGem_SI_PBMC_10K/SC3_v3_NextGem_SI_PBMC_10K_analysis.tar.gz"
dest <- file.path(base_path, "pbmc/analysis.tar.gz")
download.file(url, destfile = dest, mode = "wb")

# Extract to temporary folder
temp_dir <- file.path(base_path, "pbmc/temp_analysis/graph_clust")
untar(dest, exdir = temp_dir)

# Move file
from_path <- file.path(base_path, "pbmc/temp_analysis/analysis/clustering/graphclust/clusters.csv")
to_path <- file.path(base_path, "pbmc/clusters.csv")
file.rename(from_path, to_path)

# --- Cleanup ---
unlink(temp_dir, recursive = TRUE) # Remove temp folders
unlink(dest)
