source("benchmarks/utils.R")
source("benchmarks/rnaseq-data/utils.R")

library(Seurat)
library(ggplot2)
library(disize)
library(dplyr)

# Settings ----
set.seed(67)

# Number of (logical) cores available
n_threads <- parallel::detectCores() / 2L

# Number of threads used by disize
disize_threads <- 2L

# Create cluster
cl <- parallel::makeCluster(n_threads %/% disize_threads)

# Configure library paths and export needed functions to cluster
main_lib_paths <- .libPaths()

parallel::clusterExport(
    cl = cl,
    varlist = c(
        "simulate_dataset",
        "get_disize",
        "get_mor",
        "get_tmm",
        "main_lib_paths"
    )
)

# Load needed libraries
parallel::clusterEvalQ(cl, {
    .libPaths(main_lib_paths)

    library(DESeq2)
    library(disize)
    library(edgeR)
    library(dplyr)
})

# Set future plan
future::plan(future::cluster, workers = cl)

# Number of simulations
n_sims <- 8L


# Load dataset
data_pbmc <- CreateSeuratObject(counts = Read10X(data.dir = "benchmarks/rnaseq-data/data/pbmc"), project = "PBMC_10K")

# Insert clustering metadata
clust_pbmc <- read.csv("benchmarks/rnaseq-data/data/pbmc/clusters.csv")
data_pbmc$cluster_id <- as.factor(clust_pbmc$Cluster)

# Partition cells into simulated "batches"/"samples"
data_pbmc@meta.data <- data_pbmc@meta.data |>
    dplyr::as_tibble(rownames = "barcode") |>
    group_by(cluster_id) |>
    mutate(batch_id = as.factor(sample(rep_len(1:5, n())))) |>
    ungroup()

# Construct table of cell numbers per-cluster and -batch
table_pbmc <- data_pbmc@meta.data |>
    count(cluster_id, batch_id)

# Compute weights for binomial thinning
table_pbmc$weights <- min(table_pbmc$n) / table_pbmc$n
data_pbmc@meta.data <- data_pbmc@meta.data |>
    dplyr::left_join(table_pbmc, by = c("cluster_id", "batch_id"))

# Pseudo-bulk count matrix
modmat <- model.matrix(~ 0 + factor(cluster_id):factor(batch_id), data_pbmc@meta.data)
raw_bulk_counts <- t(data_pbmc[["RNA"]]$counts %*% modmat)

# Construct metadata for bulk dataset
bulk_metadata <- data.frame(sample_id = rownames(raw_bulk_counts)) |>
    mutate(
        cluster_id = as.factor(gsub(".*cluster_id\\)([^:]+):.*", "\\1", sample_id)),
        batch_id   = as.factor(gsub(".*batch_id\\)(.+)", "\\1", sample_id))
    ) |>
    left_join(table_pbmc, by = c("cluster_id", "batch_id"))
bulk_weights <- bulk_metadata$weights

# Apply binomial thinning to the pseudo-bulked matrix
thinned_bulk_counts <- matrix(
    rbinom(n = length(raw_bulk_counts), size = as.numeric(raw_bulk_counts), prob = bulk_weights[row(raw_bulk_counts)]),
    nrow = nrow(raw_bulk_counts),
    ncol = ncol(raw_bulk_counts)
)
rownames(thinned_bulk_counts) <- rownames(raw_bulk_counts)
colnames(thinned_bulk_counts) <- colnames(raw_bulk_counts)

# Run benchmark
benchmark <- run_benchmark(n_sims, thinned_bulk_counts, bulk_metadata, ~cluster_id, disize_threads)
write.table(
    x = benchmark,
    file = "benchmarks/rnaseq-data/data/scenario-1.tsv",
    row.names = FALSE,
    sep = "\t"
)
