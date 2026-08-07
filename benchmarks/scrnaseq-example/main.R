source("benchmarks/utils.R")
source("benchmarks/scrnaseq-example/utils.R")
library(Seurat)
library(SeuratData)
library(ggplot2)

# Settings ----
set.seed(67)

# Load IFNB dataset
InstallData("ifnb")
ifnb <- LoadData("ifnb")

# Download additional metadata
ctrl <- read.table(url("https://raw.githubusercontent.com/yelabucsf/demuxlet_paper_code/master/fig3/ye1.ctrl.8.10.sm.best"), head = T, stringsAsFactors = F)
stim <- read.table(url("https://raw.githubusercontent.com/yelabucsf/demuxlet_paper_code/master/fig3/ye2.stim.8.10.sm.best"), head = T, stringsAsFactors = F)
info <- rbind(ctrl, stim)

# Merge data
info <- rbind(ctrl, stim)
info$BARCODE <- gsub(pattern = "\\-", replacement = "\\.", info$BARCODE)

# Only keep cells with high-confidence donor assignments (Singlets)
info <- info[grep(pattern = "SNG", x = info$BEST), ]

# Subset cells in original object
info <- info[info$BARCODE %in% colnames(ifnb), ]

# Remove cells with multiple barcodes
info <- info[!duplicated(info$BARCODE), ]
ifnb <- subset(ifnb, cells = info$BARCODE)

# Construct metadata for Seurat object
rownames(info) <- info$BARCODE
info <- info[, c("BEST"), drop = F]
names(info) <- c("donor_id")
info$donor_id <- as.factor(info$donor_id)

ifnb <- Seurat::AddMetaData(ifnb, metadata = info)

# Pseudo-bulk
pseudo_ifnb <- Seurat::AggregateExpression(ifnb, assays = "RNA", return.seurat = T, group.by = c("stim", "donor_id", "seurat_annotations"))

# Record number of cells contributing to each pseudo-bulk profile
ifnb$pb_key <- paste(ifnb$stim, ifnb$donor_id, ifnb$seurat_annotations, sep = "_")
pseudo_ifnb$pb_key <- paste(pseudo_ifnb$stim, pseudo_ifnb$donor_id, pseudo_ifnb$seurat_annotations, sep = "_")
pseudo_ifnb@meta.data$seurat_annotations <- as.factor(pseudo_ifnb@meta.data$seurat_annotations)
pseudo_ifnb@meta.data$stim <- as.factor(pseudo_ifnb@meta.data$stim)

# Calculate cell counts per group from the original single-cell data
cell_counts_map <- table(ifnb$pb_key)

# Map the counts directly back into the pseudo-bulk metadata
pseudo_ifnb$ncells <- as.numeric(cell_counts_map[pseudo_ifnb$pb_key])

# Perform normalization
disize_sf <- disize::disize(~ stim * seurat_annotations, t(pseudo_ifnb[["RNA"]]$counts), pseudo_ifnb@meta.data, "donor_id", offset = log(pseudo_ifnb$ncells), n_feats = 5000L, n_threads = 8L)
mor_sf <- get_mor(list(counts = t(pseudo_ifnb[["RNA"]]$counts), metadata = pseudo_ifnb@meta.data))
tmm_sf <- get_tmm(list(counts = t(pseudo_ifnb[["RNA"]]$counts), metadata = pseudo_ifnb@meta.data))

# Perform differential expression analysis with each size factor
sf_list <- list(
  disize = disize_sf[pseudo_ifnb$donor_id] + log(pseudo_ifnb$ncells),
  mor    = mor_sf,
  tmm    = tmm_sf
)

de_results <- lapply(sf_list, function(sf) {
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = pseudo_ifnb[["RNA"]]$counts,
    colData   = pseudo_ifnb@meta.data,
    design    = ~ 0 + stim:seurat_annotations
  )
  DESeq2::sizeFactors(dds) <- exp(sf)

  # Run DESeq
  dds <- DESeq2::DESeq(
    dds,
    quiet   = TRUE,
    fitType = "mean"
  )

  dds
})

cell_types <- gsub(" ", ".", levels(pseudo_ifnb$seurat_annotations))

# Record DEGs across selected comparisons
deg_counts_list <- lapply(de_results, function(cur_dds) {
  sapply(cell_types, function(ct) {
    # Extract results as a standard data frame
    res <- DESeq2::results(cur_dds, contrast = list(c(paste0("stimCTRL.seurat_annotations", ct)), paste0("stimSTIM.seurat_annotations", ct)))

    # Filter for significant DEGs (handling NA values safely)
    deg_subset <- res[!is.na(res$padj) & res$padj < 0.05 & abs(res$log2FoldChange) > 1, ]

    # Return total count of DEGs for this contrast
    nrow(deg_subset)
  })
})

# Construct dataframe with counts of DEGs
deg_counts <- data.frame(
  disize = deg_counts_list$disize,
  mor = deg_counts_list$mor,
  tmm = deg_counts_list$tmm,
  cell_type = cell_types
)

# Construct DFs for relative and absolute measurements
rel_counts <- dplyr::mutate(deg_counts,
  mor = disize - mor,
  tmm = disize - tmm,
  disize = 0,
  type = "relative"
) |>
  tidyr::pivot_longer(
    cols = c(disize, mor, tmm),
    names_to = "method",
    values_to = "num_of_degs"
  )

abs_counts <- dplyr::mutate(deg_counts,
  type = "absolute"
) |>
  tidyr::pivot_longer(
    cols = c(disize, mor, tmm),
    names_to = "method",
    values_to = "num_of_degs"
  )

benchmark <- dplyr::bind_rows(rel_counts, abs_counts)
write.table(
  x = benchmark,
  file = "benchmarks/scrnaseq-example/data/deg_count.tsv",
  row.names = FALSE,
  sep = "\t"
)

# Identify DEGs that disize finds but the other methods dont for CD14 monocytes
res_list <- lapply(de_results, function(cur_dds) {
  # Construct the correct contrast matching your design formula
  res <- DESeq2::results(cur_dds, contrast = list("stimSTIM.seurat_annotationsCD14.Mono", "stimCTRL.seurat_annotationsCD14.Mono"))

  # Convert to a data frame and keep gene names
  df <- as.data.frame(res)
  df$gene <- rownames(df)
  df
})

deg_list <- dplyr::bind_rows(res_list, .id = "method")
write.table(
  x = deg_list,
  file = "benchmarks/scrnaseq-example/data/deg_list.tsv",
  row.names = FALSE,
  sep = "\t"
)
