# Split formula into fixed- and random-effects
split_formula <- function(design_formula) {
    # Extract terms
    terms <- attr(terms(design_formula), "term.labels")

    # Identify random effects
    re <- grepl("\\| ", terms)

    # Separate fixed- and random-effects terms
    fixed <- NULL
    if (!all(re)) {
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

# Get disize's size factor estimate
get_disize <- function(dataset, design_formula, n_threads = 1L, n_feats = 10000L, verbose = 1L) {
    disize_sf <- disize::disize(
        design_formula,
        dataset$counts,
        dataset$metadata,
        "batch_id",
        n_feats = n_feats,
        n_threads = n_threads,
        rel_tol = 1000,
        verbose = verbose
    )

    disize_sf
}

# Get MoR's size factor estimate
get_mor <- function(dataset) {
    dds <- DESeq2::DESeqDataSetFromMatrix(
        countData = t(dataset$counts),
        colData = dataset$metadata,
        design = ~1
    )
    dds <- DESeq2::estimateSizeFactors(dds)

    # Extract size factors
    deseq2_sf <- DESeq2::sizeFactors(dds)

    # Scale for comparisons
    deseq2_sf <- log(deseq2_sf / sum(deseq2_sf) * length(deseq2_sf))

    deseq2_sf
}

# Get TMM's size factor estimate
get_tmm <- function(dataset) {
    dds <- edgeR::DGEList(counts = t(dataset$counts))
    dds <- edgeR::calcNormFactors(dds, method = "TMM")

    # Extract size factors
    edger_sf <- dds$samples$norm.factors * Matrix::rowSums(dataset$counts)

    # Scale for comparisons
    edger_sf <- log(edger_sf / sum(edger_sf) * length(edger_sf))

    edger_sf
}
