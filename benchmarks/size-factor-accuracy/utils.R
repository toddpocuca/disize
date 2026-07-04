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

    # Split design formula
    design <- split_formula(design_formula)

    # Construct fixed-effects model matrix
    if (!is.null(design$fixed)) {
        fe_design <- model.matrix(design$fixed, metadata)
    } else {
        fe_design <- matrix(nrow = nrow(metadata), ncol = 0L)
    }

    # Construct random-effects model matrix
    if (!is.null(design$random)) {
        remm <- reformulas::mkReTrms(
            bars = reformulas::findbars(design$random),
            fr = metadata,
            calc.lambdat = FALSE,
            sparse = TRUE
        )
        re_design <- Matrix::t(remm$Zt)
    } else {
        remm <- list(
            Ztlist = list(matrix(nrow = 0L, ncol = nrow(metadata)))
        )
        re_design <- matrix(nrow = nrow(metadata), ncol = 0L)
    }

    # Simulate size factors
    size_factors <- runif(length(levels(metadata$batch_id)), 0.1, 1.0) |>
        setNames(levels(metadata$batch_id))

    # Simulate counts
    counts <- sapply(1:n_genes, function(g_i) {
        # Simulate fixed-effects
        sparse <- runif(ncol(fe_design)) < sparsity
        fe_coefs <- ifelse(sparse, 0.0, rnorm(ncol(fe_design), 0, mgt))

        # Simulate random-effects
        re_coefs <- lapply(lapply(remm$Ztlist, nrow), function(n_coefs) {
            sparse <- rep(runif(1) < sparsity, n_coefs)

            ifelse(sparse, 0.0, rnorm(n_coefs, sd = mgt))
        }) |>
            unlist()

        # Compute realized magnitude
        log_mu <- as.vector(
            rnorm(1, log(avg)) +
                fe_design %*% fe_coefs +
                re_design %*% re_coefs +
                log(size_factors[metadata$batch_id])
        )

        # Draw counts from negative binomial
        rnbinom(nrow(metadata), mu = exp(log_mu), size = iodisps[g_i])
    })

    # Cast to integers
    mode(counts) <- "integer"

    list(
        counts = counts,
        metadata = metadata,
        size_factors = log(
            size_factors / sum(size_factors) * length(size_factors)
        )
    )
}

# Run a benchmark with specified simulation settings
run_benchmark <- function(n_sims, sim_pars, design_formula, metadata, disize_threads = 1L) {
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
                    errors <- future.apply::future_replicate(
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
                            disize_sf <- get_disize(dataset, design_formula, disize_threads)
                            mor_sf <- get_mor(dataset)
                            tmm_sf <- get_tmm(dataset)

                            # Compute absolute error on log-scale
                            tibble::tibble(
                                disize = (dataset$size_factors - disize_sf)^2,
                                mor = (dataset$size_factors - mor_sf)^2,
                                tmm = (dataset$size_factors - tmm_sf)^2
                            ) |>
                                dplyr::summarise(
                                    dplyr::across(
                                        .cols = tidyr::everything(),
                                        .fns = ~ sqrt(sum(.x))
                                    )
                                )
                        },
                        simplify = FALSE
                    )

                    # Denote simulation ID
                    dplyr::bind_rows(errors) |>
                        tibble::rowid_to_column("sim_id")
                }
            )
        ) |>
        tidyr::unnest(d)

    # Compute absolute error
    abs_benchmark <- tidyr::pivot_longer(benchmark,
        cols = c(disize, mor, tmm),
        names_to = "method",
        values_to = "error"
    ) |>
        dplyr::group_by(n_genes, sparsity, mgt, avg, setting_id, method) |>
        dplyr::summarise(
            q95 = quantile(error, 0.95),
            q75 = quantile(error, 0.75),
            q60 = quantile(error, 0.60),
            q50 = quantile(error, 0.50),
            q40 = quantile(error, 0.40),
            q25 = quantile(error, 0.25),
            q5 = quantile(error, 0.05),
            type = "absolute"
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
        dplyr::group_by(n_genes, sparsity, mgt, avg, setting_id, method) |>
        dplyr::summarise(
            q95 = quantile(error, 0.95),
            q75 = quantile(error, 0.75),
            q60 = quantile(error, 0.60),
            q50 = quantile(error, 0.50),
            q40 = quantile(error, 0.40),
            q25 = quantile(error, 0.25),
            q5 = quantile(error, 0.05),
            type = "relative"
        )


    # Combine dataframes
    benchmark <- dplyr::bind_rows(abs_benchmark, rel_benchmark)
}
