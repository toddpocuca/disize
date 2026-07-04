library(parallel)
library(future)
library(future.apply)
source("benchmarks/size-factor-accuracy/utils.R")

# Settings ----
set.seed(67)

# Number of (logical) cores available
n_cores <- parallel::detectCores() / 2L

# Number of threads used by disize
disize_threads <- 2L

# Create cluster
cl <- parallel::makeCluster(n_cores %/% disize_threads)

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

    library(disize)
    library(dplyr)
    library(purrr)
    library(tidyr)
})

# Set future plan
future::plan(future::cluster, workers = cl)

# Number of simulations per setting
n_sims <- 100L

# Simulations settings
sim_pars <- expand.grid(
    "n_genes" = c(10000L),
    "sparsity" = c(0.2, 0.8),
    "mgt" = c(0.5, 2.0),
    "avg" = c(10.0, 25, 50, 75, 100, 250, 500, 1000)
)
sim_pars$setting_id <- seq_len(nrow(sim_pars))

# A Trivial Case ----
# Define data generating process
design_formula <- ~ (1 | donor_id)

n_donors <- 8L
metadata <- data.frame(
    donor_id = factor(1:n_donors),
    batch_id = factor(1:n_donors)
)

# Run benchmark
benchmark <- run_benchmark(n_sims, sim_pars, design_formula, metadata, disize_threads)
write.table(
    x = benchmark,
    file = "benchmarks/size-factor-accuracy/data/scenario-1.tsv",
    row.names = FALSE,
    sep = "\t"
)


# Comparing Two Conditions ----
# Define data generating process
design_formula <- ~cond_id

n_donors <- 10L
metadata <- data.frame(
    donor_id = factor(1:n_donors),
    cond_id = cut(1:n_donors, 2L),
    batch_id = factor(1:n_donors)
)
levels(metadata$cond_id) <- letters[1:2]

# Run benchmark
benchmark <- run_benchmark(n_sims, sim_pars, design_formula, metadata, disize_threads)
write.table(
    x = benchmark,
    file = "benchmarks/size-factor-accuracy/data/scenario-2.tsv",
    row.names = FALSE,
    sep = "\t"
)

# Multiple Factors ----
# Define data generating process
design_formula <- ~ cond_id:sex_id

n_donors <- 12L
metadata <- data.frame(
    donor_id = factor(1:n_donors),
    cond_id = cut(1:n_donors, 2L),
    sex_id = rep(cut(1:(n_donors / 2L), 2L), 2L),
    batch_id = factor(1:n_donors)
)
levels(metadata$cond_id) <- letters[1L:2L]
levels(metadata$sex_id) <- c("female", "male")


# Run benchmark
benchmark <- run_benchmark(n_sims, sim_pars, design_formula, metadata, disize_threads)
write.table(
    x = benchmark,
    file = "benchmarks/size-factor-accuracy/data/scenario-3.tsv",
    row.names = FALSE,
    sep = "\t"
)
