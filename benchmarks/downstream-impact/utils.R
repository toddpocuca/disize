source("benchmarks/utils.R")

# Simulate a single dataset given its settings
simulate_dataset <- function(
    design_formula,
    metadata,
    n_genes,
    sparsity,
    mgt,
    avg) {
    # Simulate inverse overdispersion factors
    iodisps <- rlnorm(n_genes, log(10), 0.5)

    # Simulate size factors
    size_factors <- runif(length(levels(metadata$batch_id)), 0.1, 1.0) |>
        setNames(levels(metadata$batch_id))
    batch_effect <- size_factors[metadata$batch_id]

    # Construct sparsity mask
    mask <- rep(1.0, n_genes)
    mask[1:floor(sparsity * n_genes)] <- 0.0

    # Construct design matrix
    design <- model.matrix(
        split_formula(design_formula)$fixed, # Ignore any random-effects
        metadata
    )

    # Simulate counts
    counts <- sapply(1:n_genes, function(g_i) {
        # Simulate effects
        effects <- rnorm(ncol(design), sd = mgt) * mask[g_i]

        # Compute realized magnitude
        log_mu <- as.vector(
            rnorm(1, log(avg)) +
                design %*% effects +
                log(batch_effect)
        )

        # Draw counts from negative binomial
        rnbinom(nrow(design), mu = exp(log_mu), size = iodisps[g_i])
    })

    # Cast to integers
    mode(counts) <- "integer"

    list(
        counts = counts,
        metadata = metadata,
        size_factors = log(
            size_factors / sum(size_factors) * length(size_factors)
        ),
        null = mask == 0.0
    )
}

# Run a benchmark with specified simulation settings
run_benchmark <- function(
    n_sims,
    sim_pars,
    design_formula,
    metadata) {
    benchmark <- sim_pars |>
        dplyr::mutate(
            d = purrr::pmap(
                list(n_genes, sparsity, mgt, avg),
                function(n_genes, sparsity, mgt, avg) {
                    message(
                        "n_genes: ", n_genes,
                        ", sparsity: ", sparsity,
                        ", mgt: ", mgt,
                        ", avg: ", avg
                    )
                    # Simulate
                    res_data <- future.apply::future_replicate(
                        n = n_sims,
                        expr = {
                            # Simulate dataset
                            dataset <- simulate_dataset(
                                design_formula,
                                metadata,
                                n_genes,
                                sparsity,
                                mgt,
                                avg
                            )

                            # Estimate size factors
                            size_factors <- list(
                                ground_truth = dataset$size_factors,
                                disize = get_disize(
                                    dataset,
                                    design_formula,
                                    2L
                                ),
                                mor = get_mor(dataset),
                                tmm = get_tmm(dataset)
                            )

                            # Run DESeq2 with various
                            res_data <- lapply(size_factors, function(sf) {
                                dds <- DESeq2::DESeqDataSetFromMatrix(
                                    countData = Matrix::t(dataset$counts),
                                    colData = dataset$metadata,
                                    design = design_formula
                                )
                                DESeq2::sizeFactors(dds) <- exp(sf)

                                # Run DESeq
                                dds <- DESeq2::DESeq(
                                    dds,
                                    quiet = TRUE,
                                    fitType = "mean"
                                )
                                results <- DESeq2::results(dds)

                                # Format
                                data.frame(
                                    p_values = results$pvalue,
                                    null = dataset$null
                                )
                            }) |>
                                dplyr::bind_rows(.id = "method") |>
                                dplyr::mutate(method = factor(method))

                            res_data
                        },
                        simplify = FALSE
                    ) |>
                        dplyr::bind_rows(.id = "sim_id")

                    # Iterate over methods
                    errors <- lapply(levels(res_data$method), function(method) {
                        cur_data <- res_data[res_data$method == method, ]

                        # Compute conditional type 1 error
                        type_1 <- mean(
                            cur_data$p_values[cur_data$null] < 0.05,
                            na.rm = TRUE
                        )

                        # Compute conditional type 2 error
                        type_2 <- mean(
                            cur_data$p_values[!cur_data$null] > 0.05,
                            na.rm = TRUE
                        )

                        data.frame(
                            type_1 = type_1,
                            type_2 = type_2
                        )
                    })
                    names(errors) <- levels(res_data$method)

                    errors <- errors |>
                        dplyr::bind_rows(.id = "method")

                    errors
                }
            )
        ) |>
        tidyr::unnest(d) |>
        dplyr::group_by(setting_id) |>
        dplyr::mutate(
            type_1_relative = type_1 / type_1[method == "ground_truth"],
            type_2_relative = type_2 / type_2[method == "ground_truth"]
        ) |>
        dplyr::ungroup()

    benchmark
}
