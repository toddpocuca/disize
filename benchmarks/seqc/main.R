source("benchmarks/utils.R")

library(seqc)
library(DESeq2)
library(edgeR)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)


# Construct ground-truth DEG table ----
## Load in taqman data
data("taqman", package = "seqc")

## Construct dataframe with ground-truth differential expression
diffexp_data <- taqman |>
  tidyr::pivot_longer(
    -c(EntrezID, Symbol),
    names_to = c("sample", ".value"),
    names_pattern = "([A-D][0-9])_(value|detection)"
  ) |>
  dplyr::filter(!is.na(EntrezID)) |>
  dplyr::filter(detection != "A") |>
  dplyr::mutate(
    sample_type = substr(sample, 1, 1),
    EntrezID = as.character(EntrezID),
    Symbol = as.character(Symbol)
  ) |>
  dplyr::rename(entrez_id = EntrezID, gene_symbol = Symbol) |>
  dplyr::group_by(entrez_id, gene_symbol, sample_type) |>
  dplyr::summarise(log_value = mean(log2(value)), .groups = "drop") |>
  tidyr::pivot_wider(names_from = sample_type, values_from = log_value)

## Evaluate contrast between A and B
contr_data <- diffexp_data |>
  dplyr::mutate(diff = A - B) |>
  dplyr::filter(!is.na(diff)) |>
  dplyr::select(c(entrez_id, gene_symbol, diff))


# Construct ERCC and endogenous count matrices ----
## Load in count data
data("ILM_refseq_gene_BGI", package = "seqc")

## Grab full count (obs x feat) matrix
counts <- t(as.matrix(ILM_refseq_gene_BGI[, -(1:4)]))
colnames(counts) <- ILM_refseq_gene_BGI$EntrezID

## Grab gene-level metadata
feat_md <- ILM_refseq_gene_BGI[, 1:4]

## Construct sample-level metadata
obs_md <- data.frame(obs_id = colnames(ILM_refseq_gene_BGI)[-(1:4)]) |>
  tidyr::extract(
    obs_id,
    into = c("sample_type", "replicate_id", "lane_id", "flowcell_id"),
    regex = "^([A-F])_([0-9]+)_L([0-9]+)_FlowCell([A-Z])$",
    remove = FALSE
  ) |>
  dplyr::mutate(
    batch_id = interaction(sample_type, replicate_id, flowcell_id)
  )

## Partition counts and metadata into ERCC-derived and endogenous
ercc_counts <- counts[, feat_md$IsERCC]
endo_counts <- counts[, !feat_md$IsERCC]

## Subset observations from samples A and B
subset_obs <- obs_md$sample_type %in% c("A", "B")

obs_md <- obs_md[subset_obs, ]
ercc_counts <- ercc_counts[subset_obs, ]
endo_counts <- endo_counts[subset_obs, ]

## Re-label factors in obs_md
obs_md <- obs_md |>
  dplyr::mutate(
    sample_type = factor(sample_type),
    batch_id = factor(batch_id)
  )


# Construct datasets and estimate size factors ----
## Construct dataset
endo_dataset <- list(
  counts = endo_counts,
  metadata = obs_md
)

## Estimate size factors from endogenous counts
disize_sf <- get_disize(endo_dataset, design_formula = ~sample_type, n_threads = 8L)
mor_sf <- get_mor(endo_dataset)
tmm_sf <- get_tmm(endo_dataset)

sf_list <- list(disize = disize_sf[obs_md$batch_id], mor = mor_sf, tmm = tmm_sf)


# Compute differential expression with size factors ----
run_de <- function(sf) {
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = t(endo_counts),
    colData = obs_md,
    design = ~ 0 + sample_type
  )
  DESeq2::sizeFactors(dds) <- exp(sf)
  dds <- DESeq2::DESeq(dds, quiet = TRUE, fitType = "mean")
  res <- DESeq2::results(dds, contrast = c("sample_type", "A", "B"))

  tibble(
    entrez_id = rownames(res),
    log2fc = res$log2FoldChange,
    padj = res$padj
  )
}

## Merge and filter results
de_results <- lapply(sf_list, run_de) |>
  bind_rows(.id = "method") |>
  dplyr::filter(entrez_id %in% contr_data$entrez_id)


# Benchmark against taqman data ----
## Establish significance thresholds
lfc_threshold <- 1.0
padj_threshold <- 0.05

## Extract the ground-truth DEG pool vector directly from contr_data
taqman_deg_pool <- contr_data$entrez_id[abs(contr_data$diff) >= lfc_threshold]

## Evaluate method performance using vector matching
benchmark_summary <- de_results |>
  dplyr::mutate(
    is_method_deg = !is.na(padj) & padj < padj_threshold & abs(log2fc) > lfc_threshold,
    is_taqman_deg = entrez_id %in% taqman_deg_pool
  ) |>
  dplyr::left_join(dplyr::select(contr_data, entrez_id, gene_symbol), by = "entrez_id") |>
  dplyr::select(method, entrez_id, gene_symbol, log2fc, padj, is_method_deg, is_taqman_deg)

## Export DEG list
write.table(
  x = benchmark_summary,
  file = "benchmarks/seqc/data/deg_list.tsv",
  row.names = FALSE,
  sep = "\t"
)
