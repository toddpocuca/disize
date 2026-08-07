library(utils)

process_dataset <- function(dataset) {
    gene_names <- dataset$V1
    counts <- as.matrix(dataset[, -1])
    rownames(counts) <- gene_names
}

base_path <- "benchmarks/scrnaseq-example/data"
dir.create(base_path, recursive = TRUE, showWarnings = FALSE)

# Process Stimulated PBMCs Dataset ---
# Download count matrix file
url <- "https://www.dropbox.com/s/79q6dttg8yl20zg/immune_alignment_expression_matrices.zip?dl=1"

dest <- file.path(base_path, "stim_pbmc_data.zip")
download.file(url, destfile = dest, mode = "wb")

# Extract to temporary folder
temp_dir <- file.path(base_path, "temp_stim_pbmc")
unzip(dest, exdir = temp_dir)

# Concatenate to single dataset
ctrl_data <- data.table::fread(
    file = paste0(base_path, "/temp_stim_pbmc/immune_control_expression_matrix.txt.gz"), sep = "\t"
)
stim_data <- data.table::fread(
    file = paste0(base_path, "/temp_stim_pbmc/immune_stimulated_expression_matrix.txt.gz"), sep = "\t"
)

ctrl_data <- as.matrix(ctrl_data[, -1])

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
