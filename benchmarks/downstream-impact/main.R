library(dplyr)
library(parallel)
library(future)
library(future.apply)
source("benchmarks/downstream-impact/utils.R")

# Create cluster
cl <- parallel::makeCluster(parallel::detectCores() / 2L - 1L)

# Configure library paths and export needed functions to cluster
main_lib_paths <- .libPaths()

parallel::clusterExport(
    cl = cl,
    varlist = c(
        "split_formula",
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
n_sims <- 30L

# Configure simulation settings
sim_pars <- expand.grid(
    "n_genes" = c(10000L),
    "sparsity" = c(0.15, 0.20, 0.25),
    "mgt" = c(2),
    "avg" = c(10.0, 25, 50, 75, 100.0)
) |>
    tibble::rowid_to_column("setting_id")

# Comparing Two Conditions ----
# Construct experimental design
design_formula <- ~cond_id

n_donors <- 10L
metadata <- data.frame(
    donor_id = factor(1:n_donors),
    cond_id = cut(1:n_donors, 2L),
    batch_id = factor(1:n_donors)
)
levels(metadata$cond_id) <- letters[1:2]

# Run benchmark
benchmark <- run_benchmark(
    n_sims,
    sim_pars,
    design_formula,
    metadata
)
write.table(
    x = benchmark,
    file = "benchmarks/downstream-impact/data/scenario-1.tsv",
    row.names = FALSE,
    sep = "\t"
)
