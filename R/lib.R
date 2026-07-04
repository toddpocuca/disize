#' Split design formula into fixed and random effects.
#'
#' @param design_formula The design formula
#'
#' @noRd
split_formula <- function(design_formula) {
    # Extract terms
    terms <- attr(terms(design_formula), "term.labels")

    # Identify random effects
    re <- grepl("\\| ", terms)

    # Separate fixed- and random-effects terms
    fixed <- NULL
    if (any(!re)) {
        fixed <- stats::formula(
            paste0(" ~ 0 + ", paste(terms[!re], collapse = " + "))
        )
    }

    random <- NULL
    if (any(re)) {
        random <- stats::formula(paste0(
            " ~ 0 + ", paste(terms[re], collapse = " + ")
        ))
    }

    list(
        formula = design_formula,
        fixed = fixed,
        random = random
    )
}

#' @title Design-informed size factor estimation
#'
#' @param design_formula The formula describing the experimental design.
#' @param counts A (observation x feature) count matrix.
#' @param metadata A dataframe containing observation-level metadata.
#' @param batch_name The identifier for the batch column in 'metadata'.
#' @param offset An optional offset used to adjust expression estimates for each observation.
#' @param obs_name The identifier for the observation column in 'metadata' (useful if count matrix and metadata do not share ordering).
#' @param n_feats The number of features used during estimation.
#' @param n_iters The number of iterations used for estimation.
#' @param rel_tol The relative tolerance used for convergence.
#' @param init_alpha The initial step-size for the optimizer, lower values
#'  can sometimes make it easier to estimate size factors for more complex
#'  designs.
#' @param n_threads The number of threads to use for parallel processing.
#' @param n_tries The maximum number of times to try fitting.
#' @param verbose The verbosity level (`1`: only errors, `2`: also allows
#'  warnings,`3`: also allows messages, `4`: also prints additional
#' output useful for debugging).
#'
#' @returns A named numeric vector containing the size factor estimates.
#'
#' @export
disize <- function(
    design_formula,
    counts,
    metadata,
    batch_name,
    offset = NULL,
    obs_name = "obs_id",
    n_feats = min(10000L, ncol(counts)),
    n_iters = 10000L,
    rel_tol = 1000,
    init_alpha = 1e-8,
    n_threads = 1L,
    n_tries = 5L,
    verbose = 3L) {
    # Argument Checks ----
    # Check design formula is correct
    if (!methods::is(design_formula, "formula")) {
        stop("'design_formula' should be an R formula")
    } else if (2 < length(design_formula)) {
        stop("'design_formula' should be of the form '~ x + ...'")
    }

    # Check for the same number of samples
    if (nrow(metadata) != nrow(counts)) {
        stop(
            "'counts' and 'metadata' should have the same # of ",
            "observations (rows)."
        )
    }

    # Check offset has the same number of elements
    if (!is.null(offset) && length(offset) != nrow(counts)) {
        stop(
            "'offset' and 'counts' must have the same # of observations."
        )
    }

    # Formatting Data ----
    if (3L <= verbose) {
        message("Formatting data...")
    }

    # Include explicit observation names if not present
    if (is.null(rownames(counts)) && is.null(metadata[[obs_name]])) {
        rownames(counts) <- seq_len(nrow(counts))
        metadata[[obs_name]] <- seq_len(nrow(counts))
    } else if (!is.null(rownames(counts)) && is.null(metadata[[obs_name]])) {
        metadata[[obs_name]] <- rownames(counts)
    } else if (is.null(rownames(counts)) && !is.null(metadata[[obs_name]])) {
        rownames(counts) <- metadata[[obs_name]]
    }

    # Include explicit feature names if not present
    if (is.null(colnames(counts))) {
        colnames(counts) <- seq_len(ncol(counts))
    }

    # Re-order counts to match metadata (if needed)
    counts <- counts[metadata[[obs_name]], ]

    # Filter out features with no counts
    subset <- Matrix::colSums(counts) != 0
    counts <- counts[, subset]

    # Ensure valid number of features selected
    if (3L <= verbose && ncol(counts) < n_feats) {
        message(
            "Insufficient number of features (",
            ncol(counts),
            ") after ",
            "subsetting observations to satisfy n_feats = ",
            n_feats,
            ". "
        )
    }
    n_feats <- min(n_feats, ncol(counts))

    # Preferentially subset features with sufficient counts
    ordering <- order(Matrix::colMeans(counts), decreasing = TRUE)
    counts <- counts[, ordering[1:n_feats]]

    # Cast to dense matrix if needed
    counts <- base::as.matrix(counts)

    # Ensure batch identifier is a factor variable
    metadata[[batch_name]] <- as.factor(metadata[[batch_name]])

    # Allocate data for Stan ----
    stan_data <- list(n_obs = nrow(counts), n_feats = ncol(counts))
    stan_data[["n_batches"]] <- nlevels(metadata[[batch_name]])
    stan_data[["batch_id"]] <- as.integer(metadata[[batch_name]])

    # Include offsets if necessary
    if (!is.null(offset)) {
        stan_data[["n_of"]] <- stan_data$n_obs
        stan_data[["offsets"]] <- offset
    } else {
        stan_data[["n_of"]] <- 0L
        stan_data[["offsets"]] <- array(0, dim = c(0))
    }

    # Split the design formula into fixed- and random-effects
    design <- split_formula(design_formula)

    # Construct fixed-effects model matrix
    if (!base::is.null(design$fixed)) {
        fe_design <- model.matrix(design$fixed, metadata)

        stan_data[["n_fe"]] <- ncol(fe_design)
        stan_data[["fe_design"]] <- fe_design
    } else {
        stan_data[["n_fe"]] <- 0L
        stan_data[["fe_design"]] <- array(0, dim = c(stan_data$n_obs, 0))
    }

    # Construct random-effects model matrix
    if (!is.null(design$random)) {
        remm <- reformulas::mkReTrms(
            bars = reformulas::findbars(design$random),
            fr = metadata,
            calc.lambdat = FALSE,
            sparse = TRUE
        )
        re_design <- Matrix::t(remm$Zt) |> methods::as("RsparseMatrix")

        # Check if all random-effects terms are scalar
        all_scalar <- lapply(remm$cnms, function(b) {
            length(b) == 1
        }) |>
            unlist() |>
            all()
        if (!all_scalar) {
            stop(
                "Only include one predictor on the LHS of a random-effects bar."
            )
        }

        stan_data[["n_re"]] <- ncol(re_design)
        stan_data[["n_nz_re"]] <- length(re_design@x)
        stan_data[["re_design_x"]] <- re_design@x
        stan_data[["re_design_j"]] <- re_design@j + 1L
        stan_data[["re_design_p"]] <- re_design@p + 1L
        stan_data[["n_re_terms"]] <- length(remm$cnms)
        stan_data[["re_id"]] <- rep(seq_along(remm$cnms), times = diff(remm$Gp))
    } else {
        stan_data[["n_re"]] <- 0L
        stan_data[["n_nz_re"]] <- 0L
        stan_data[["re_design_x"]] <- numeric(0)
        stan_data[["re_design_j"]] <- integer(0)
        stan_data[["re_design_p"]] <- integer(stan_data[["n_obs"]] + 1)
        stan_data[["n_re_terms"]] <- 0L
        stan_data[["re_id"]] <- integer(0)
    }

    # Include counts and grainsize for Stan
    stan_data[["counts"]] <- base::t(counts)
    stan_data[["grainsize"]] <- ceiling(
        stan_data[["n_feats"]] / n_threads
    )

    # Construct Stan model
    model <- instantiate::stan_package_model(
        name = "disize",
        package = "disize"
    )

    # TEMPORARY UNTIL CMDSTANR CPP_OPTIONS REFACTOR ----
    cpp_options <- model$.__enclos_env__$private$cpp_options_
    cpp_options$stan_threads <- TRUE
    model$.__enclos_env__$private$cpp_options_ <- cpp_options
    # TEMPORARY UNTIL CMDSTANR CPP_OPTIONS REFACTOR ----

    if (3L <= verbose) {
        message(
            "Optimizing over initialization..."
        )
    }
    # Optimize over initialization
    inits <- list()
    grads <- numeric(3)
    for (i in 1:n_tries) {
        fit <- model$optimize(
            data = stan_data,
            iter = 200L,
            threads = n_threads,
            algorithm = "lbfgs",
            init_alpha = init_alpha,
            history_size = 10L,
            tol_rel_obj = rel_tol,
            tol_rel_grad = rel_tol,
            sig_figs = 16L,
            show_messages = FALSE
        )

        # Grab output
        output <- capture.output(fit$output())
        output <- output[length(output) - 2]

        # Extract final iteration's gradient
        final_gradient <- strsplit(
            x = trimws(gsub("\\s+", " ", output)),
            split = " "
        )[[1]][4] |>
            as.numeric()

        # Store fit and gradient
        inits[[i]] <- fit
        grads[[i]] <- final_gradient
    }

    # Order fits with respect to final gradient
    fit_order <- order(grads)

    # Compute estimated maximum ETA
    max_eta <- mean(sapply(inits, function(x) {
        x$time()$total
    })) / 200 * n_iters

    # Estimate model parameters ----
    options(cmdstanr_warn_inits = FALSE)
    if (3L <= verbose) {
        message(
            "Estimating size factors... (Max ETA: ~",
            round(max_eta, 1),
            "s)"
        )
    }

    # Estimate fit
    for (i in 1:n_tries) {
        # Try fitting model
        fit <- tryCatch(
            expr = {
                model$optimize(
                    data = stan_data,
                    init = inits[[fit_order[i]]],
                    iter = n_iters,
                    threads = n_threads,
                    algorithm = "lbfgs",
                    init_alpha = init_alpha,
                    history_size = 10L,
                    tol_rel_obj = rel_tol,
                    tol_rel_grad = rel_tol,
                    sig_figs = 16L,
                    show_messages = (4L <= verbose),
                    refresh = ceiling(n_iters / 10)
                )
            },
            error = function(err) {
                NULL
            }
        )

        # Handle error case
        if (!is.null(fit)) {
            # Check for convergence
            output <- utils::capture.output(fit$output())
            if (any(grepl("Convergence detected", output))) {
                # Extract size factors
                sf <- fit$mle("sf")
                names(sf) <- levels(metadata[[batch_name]])

                break
            } else if (2L <= verbose && i < n_tries) {
                warning(
                    "Model did not converge, retrying fit with different ",
                    "initialization... (if you see this multiple times try ",
                    "increasing 'n_iters')"
                )
            }
        } else if (3L <= verbose && i < n_tries) {
            message(
                "Error during estimation, retrying fit with different ",
                "initialization..."
            )
        } else {
            stop("Model did not fit without error after ", n_tries, " tries.")
        }
    }

    # Check for convergence
    output <- utils::capture.output(fit$output())
    if (3L <= verbose && any(grepl("Convergence detected", output))) {
        message("Finised in ", round(fit$time()$total, 1), "s!")
    } else if (2L <= verbose) {
        warning("Model did not converge, try increasing 'n_iters'.")
    }

    # Extract size factors
    sf <- fit$mle("sf")
    names(sf) <- levels(metadata[[batch_name]])

    if (4L <= verbose) {
        attr(sf, "fit") <- fit
    }

    sf
}
