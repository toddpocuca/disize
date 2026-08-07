source("benchmarks/utils.R")

# Simulate a single dataset
simulate_dataset <- function(
    counts,
    metadata,
    size_factors) {
    # Convert relative size factors to binomial thinning weights
    size_factors <- size_factors / max(size_factors)
    weights <- size_factors[metadata$batch_id]

    # Apply binomial thinning to count matrix
    thinned_counts <- base::matrix(
        rbinom(n = length(counts), size = as.numeric(counts), prob = weights[row(counts)]),
        nrow = nrow(counts),
        ncol = ncol(counts)
    )

    # Structure dataset
    list(
        counts = thinned_counts,
        metadata = metadata,
        size_factors = log(
            size_factors / sum(size_factors) * length(size_factors)
        )
    )
}

# Run a benchmark with specified simulation settings
run_benchmark <- function(n_sims, counts, metadata, design_formula, disize_threads = 1L) {
    # Perform differential expression analysis on original counts
    dds <- DESeq2::DESeqDataSetFromMatrix(
        countData = base::t(counts),
        colData = metadata,
        design = design_formula
    )
    DESeq2::sizeFactors(dds) <- rep(1, nrow(counts))

    # Compare clusters 1 & 2
    dds <- DESeq2::DESeq(
        dds,
        quiet = TRUE,
        fitType = "mean"
    )
    true_results <- DESeq2::results(dds, contrast = c("cluster_id", "6", "7"))

    # Define which genes are DE based on significance + log-FC
    true_results$is_de <- true_results$padj < 0.05 & abs(true_results$log2FoldChange) > 1.0
    true_results$is_de[is.na(true_results$is_de)] <- FALSE

    # Run benchmark
    benchmark <- future.apply::future_replicate(
        n = n_sims,
        expr = {
            n_batches <- nlevels(metadata$batch_id)

            # Generate size factors
            size_factors <- runif(n_batches, 0.1, 1.0)
            size_factors <- size_factors / sum(size_factors) * n_batches
            names(size_factors) <- levels(metadata$batch_id)

            # Simulate dataset
            dataset <- simulate_dataset(
                counts,
                metadata,
                size_factors
            )

            # Compute size factors
            sfs <- list(
                disize = get_disize(
                    dataset,
                    design_formula,
                    disize_threads
                )[dataset$metadata$batch_id],
                mor = get_mor(dataset),
                tmm = get_tmm(dataset)
            )

            # Compute new DE results
            de_results <- lapply(sfs, function(sf) {
                dds <- DESeq2::DESeqDataSetFromMatrix(
                    countData = base::t(counts),
                    colData = metadata,
                    design = design_formula
                )
                DESeq2::sizeFactors(dds) <- exp(sf)

                # Run DESeq and compare clusters 1 & 2
                dds <- DESeq2::DESeq(
                    dds,
                    quiet = TRUE,
                    fitType = "mean"
                )
                new_results <- DESeq2::results(dds, contrast = c("cluster_id", "6", "7"))

                new_results
            })

            # Compute absolute error on log-scale
            sf_errors <- tibble::tibble(
                disize = (log(size_factors[dataset$metadata$batch_id]) - sfs$disize)^2,
                mor = (log(size_factors[dataset$metadata$batch_id]) - sfs$mor)^2,
                tmm = (log(size_factors[dataset$metadata$batch_id]) - sfs$tmm)^2
            ) |>
                dplyr::summarise(
                    dplyr::across(
                        .cols = tidyr::everything(),
                        .fns = ~ sqrt(sum(.x))
                    )
                )

            # Compute estimation error
            de_est_errors <- lapply(de_results, function(cur_results) {
                # Evaluate difference
                total_error <- (true_results$log2FoldChange - cur_results$log2FoldChange)^2 |>
                    sum(na.rm = TRUE) |>
                    sqrt()
                avg_error <- total_error / sum(!is.na(cur_results$log2FoldChange))

                avg_error
            })

            # Compute Type 2 error
            de_t2_errors <- lapply(de_results, function(cur_results) {
                # Find new DE genes
                cur_results$is_de <- cur_results$padj < 0.05 & abs(cur_results$log2FoldChange) > 1.0
                cur_results$is_de[is.na(cur_results$is_de)] <- FALSE

                # Compute overlap
                t2_error <- 1 - sum(true_results$is_de & cur_results$is_de) / sum(true_results$is_de)

                t2_error
            })

            bench_df <- rbind(sf_errors, de_est_errors, de_t2_errors)
            bench_df$value <- c("size_factor", "expr_est", "type_2")

            bench_df
        },
        simplify = FALSE
    ) |>
        dplyr::bind_rows(.id = "sim_id")

    # Compute absolute error
    abs_benchmark <- tidyr::pivot_longer(benchmark,
        cols = c(disize, mor, tmm),
        names_to = "method",
        values_to = "error"
    ) |>
        dplyr::mutate(
            type = factor("absolute", c("absolute", "relative"))
        )

    # Compute relative error
    rel_benchmark <- dplyr::mutate(benchmark,
        mor = mor / disize,
        tmm = tmm / disize,
        disize = disize / disize
    ) |>
        tidyr::pivot_longer(
            cols = c(disize, mor, tmm),
            names_to = "method",
            values_to = "error"
        ) |>
        dplyr::mutate(
            type = factor("relative", c("absolute", "relative"))
        )

    benchmark <- rbind(abs_benchmark, rel_benchmark)

    benchmark
}
