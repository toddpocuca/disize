library(ggplot2)
library(dplyr)
library(patchwork)
library(ggthemes)
library(ggrepel)

my_palette <- c(
    "#EE7733", # Vibrant Orange
    "#0077BB", # Blue
    "#EE3377", # Magenta
    "#33BBEE", # Cyan
    "#CC3311", # Red
    "#009988", # Teal
    "#BBBBBB" # Grey
)

# Figure 1: Size Factor Accuracy ----
## Read in benchmarking data from the two-condition scenario
benchmark <- read.table(
    "benchmarks/size-factor-accuracy/data/scenario-2.tsv",
    sep = "\t",
    header = TRUE
) |>
    dplyr::mutate(
        method = dplyr::case_match(
            method,
            "mor" ~ "MoR (DESeq/DESeq2)",
            "tmm" ~ "TMM (edgeR)",
            .default = method
        )
    ) |>
    dplyr::filter(sparsity == 0.2 & mgt == 2)

## Construct plot displaying relative error
my_plot <- ggplot(
    data = benchmark |> dplyr::filter(method != "disize", type == "relative"),
    mapping = aes(
        x = avg,
        fill = method
    )
) +
    geom_ribbon(aes(ymin = q5, ymax = q95), alpha = 0.1) +
    geom_ribbon(aes(ymin = q40, ymax = q60), alpha = 0.2) +
    geom_line(aes(x = avg, y = q50, color = method)) +
    scale_y_log10(breaks = c(1, 5, 10, 25)) +
    scale_x_log10() +
    labs(
        x = "Average Baseline Expression Level",
        y = "Error Relative to disize",
        color = "Method",
        fill = "Method",
        linetype = "Method"
    ) +
    geom_hline(yintercept = 1.0, linetype = "dashed") +
    scale_colour_manual(values = my_palette[-1]) +
    scale_fill_manual(values = my_palette[-1]) +
    theme_classic() +
    theme(
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12)
    )
ggplot2::ggsave(
    filename = "paper/figures/fig-1-sfa.png",
    plot = my_plot,
    width = 18,
    height = 15,
    units = "cm"
)

# Figure 2: Downstream Impact ----
## Read in benchmarking data
benchmark <- read.table(
    here::here("benchmarks/downstream-impact/data/scenario-1.tsv"),
    sep = "\t",
    header = TRUE
) |>
    dplyr::mutate(
        method = dplyr::case_match(
            method,
            "mor" ~ "MoR (DESeq/DESeq2)",
            "tmm" ~ "TMM (edgeR)",
            "ground_truth" ~ "Ground Truth",
            .default = method
        )
    )

## Filter data for plot
plot_data <- benchmark |> dplyr::filter(sparsity == 0.25)

## Panel A: Type 1 Error
plot_A <- ggplot(data = plot_data, mapping = aes(x = avg, y = type_1, color = method)) +
    geom_line() +
    labs(
        x = "Average Baseline Expression Level",
        y = "Type I Error",
        color = "Method"
    ) +
    scale_colour_manual(values = my_palette) +
    scale_fill_manual(values = my_palette) +
    theme_classic() +
    theme(
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12)
    )

## Panel B: Type 2 Error
plot_B <- ggplot(data = plot_data, mapping = aes(x = avg, y = type_2, color = method)) +
    geom_line() +
    labs(
        x = "Average Baseline Expression Level",
        y = "Type II Error",
        color = "Method"
    ) +
    scale_colour_manual(values = my_palette) +
    scale_fill_manual(values = my_palette) +
    theme_classic() +
    theme(
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12)
    )

## Combine panels
my_plot <- (plot_A / plot_B) +
    patchwork::plot_layout(guides = "collect") +
    patchwork::plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(size = 16, face = "bold"))
ggplot2::ggsave(
    filename = "paper/figures/fig-2-di.png",
    plot = my_plot,
    width = 18,
    height = 15,
    units = "cm"
)

# Figure 3: SFA Validation With Pseudo-bulk ----
benchmark <- read.table(
    "benchmarks/rnaseq-data/data/scenario-1.tsv",
    sep = "\t",
    header = TRUE
) |>
    dplyr::mutate(
        method = dplyr::case_match(
            method,
            "mor" ~ "MoR (DESeq/DESeq2)",
            "tmm" ~ "TMM (edgeR)",
            .default = method
        )
    ) |>
    dplyr::filter(
        (value == "size_factor" & type == "relative")
    )

base_theme <- theme_classic() +
    theme(
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12)
    )

## Construct panel for size factors
plot_data <- benchmark |> dplyr::filter(method != "disize")
my_plot <- ggplot(plot_data, aes(x = method, y = error, color = method, fill = method)) +
    geom_violin(alpha = 0.5) +
    geom_boxplot(width = 0.1, color = "black", alpha = 0.7, outlier.shape = 16, outlier.size = 1.5) +
    labs(y = "Error Relative to disize") +
    scale_y_log10() +
    scale_colour_manual(values = my_palette[-1]) +
    scale_fill_manual(values = my_palette[-1]) +
    base_theme
ggplot2::ggsave(
    filename = "paper/figures/fig-3-rd_sfa.png",
    plot = my_plot,
    width = 20,
    height = 15,
    units = "cm"
)


# Figure 4: DEA Validation With Pseudo-bulk ----
benchmark <- read.table(
    "benchmarks/rnaseq-data/data/scenario-1.tsv",
    sep = "\t",
    header = TRUE
) |>
    dplyr::mutate(
        method = dplyr::case_match(
            method,
            "mor" ~ "MoR (DESeq/DESeq2)",
            "tmm" ~ "TMM (edgeR)",
            .default = method
        )
    ) |>
    dplyr::filter(
        (value %in% c("expr_est", "type_2") & type == "absolute")
    )

base_theme <- theme_classic() +
    theme(
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12)
    )

## Construct panel for expression estimates
plot_data <- benchmark |> dplyr::filter(value == "expr_est")
plot_A <- ggplot(plot_data, aes(x = method, y = error, color = method, fill = method)) +
    geom_violin(alpha = 0.5) +
    geom_boxplot(width = 0.1, color = "black", alpha = 0.7, outlier.shape = 16, outlier.size = 1.5) +
    labs(y = "Log Fold-Change Error") +
    scale_colour_manual(values = my_palette) +
    scale_fill_manual(values = my_palette) +
    base_theme

## Construct panel for type 2 error
plot_data <- benchmark |> dplyr::filter(value == "type_2")
plot_B <- ggplot(plot_data, aes(x = method, y = error, color = method, fill = method)) +
    geom_violin(alpha = 0.5) +
    geom_boxplot(width = 0.1, color = "black", alpha = 0.7, outlier.shape = 16, outlier.size = 1.5) +
    labs(y = "Type II Error") +
    scale_colour_manual(values = my_palette) +
    scale_fill_manual(values = my_palette) +
    base_theme

## Combine panels
combined_plot <- (plot_A / plot_B) +
    patchwork::plot_layout(axes = "collect") +
    patchwork::plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(size = 16, face = "bold"))
ggplot2::ggsave(
    filename = "paper/figures/fig-4-rd_dea.png",
    plot = combined_plot,
    width = 19,
    height = 15,
    units = "cm"
)


# Figure 5: Counting DEGs For The scRNA-seq Case Study ----
## Read in DEG count data
deg_counts <- read.table(
    "benchmarks/scrnaseq-example/data/deg_count.tsv",
    sep = "\t",
    header = TRUE
) |>
    dplyr::mutate(
        method = dplyr::case_match(
            method,
            "mor" ~ "disize vs. MoR (DESeq/DESeq2)",
            "tmm" ~ "disize vs. TMM (edgeR)",
            .default = method
        )
    ) |>
    dplyr::filter(type == "relative", method != "disize")

## Construct plot for displaying DEGs found
my_plot <- ggplot(deg_counts, aes(x = method, y = num_of_degs, colour = method, fill = method)) +
    geom_violin(alpha = 0.35) +
    ggbeeswarm::geom_beeswarm(size = 2, cex = 1.5) +
    geom_hline(yintercept = 0.0, color = "black", linetype = "dashed") +
    labs(y = "Relative Number of DEGs Found By disize") +
    scale_colour_manual(values = my_palette[-1]) +
    scale_fill_manual(values = my_palette[-1]) +
    base_theme
ggplot2::ggsave(
    filename = "paper/figures/fig-5-sc.png",
    plot = my_plot,
    width = 20,
    height = 15,
    units = "cm"
)


# Figure 6: Volcano Plot For The scRNA-seq Case Study ----
## Read in DEG list for monocytes
deg_list <- read.table(
    "benchmarks/scrnaseq-example/data/deg_list.tsv",
    sep = "\t",
    header = TRUE
) |>
    dplyr::mutate(
        method = dplyr::case_match(
            method,
            "mor" ~ "MoR (DESeq/DESeq2)",
            "tmm" ~ "TMM (edgeR)",
            .default = method
        )
    )

# Categorize genes
gene_classification <- deg_list |>
    # Determine if a row qualifies as a discovered DEG
    mutate(is_deg = !is.na(padj) & padj < 0.05 & abs(log2FoldChange) > 1) |>
    group_by(gene) |>
    summarize(
        discovered_by_disize = any(method == "disize" & is_deg),
        discovered_by_others = any(method != "disize" & is_deg),
        .groups = "drop"
    ) |>
    mutate(
        discovery_status = case_when(
            discovered_by_disize & !discovered_by_others ~ "Only discovered by disize",
            discovered_by_disize & discovered_by_others ~ "Discovered by disize & others",
            TRUE ~ "Not significant"
        )
    )

# Filter data down to just the disize subset for plotting, and merge the classification categories back in
plot_data <- deg_list %>%
    filter(method == "disize") %>%
    left_join(gene_classification, by = "gene")

# Identify the top 10 unique disize discoveries to label on the plot
top_unique_labels <- plot_data %>%
    filter(discovery_status == "Only discovered by disize") %>%
    slice_min(order_by = pvalue, n = 10)

# Construct volcano plo
my_plot <- ggplot(plot_data, aes(x = log2FoldChange, y = -log10(pvalue), color = discovery_status)) +
    geom_point(alpha = 0.6, size = 1.8) +
    scale_color_manual(
        values = c(
            "Only discovered by disize" = "#e91216",
            "Discovered by disize & others" = "#377EB8",
            "Not significant" = "#CCCCCC"
        )
    ) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50", linewidth = 0.5) +
    geom_vline(xintercept = c(-1.0, 1.0), linetype = "dashed", color = "grey50", linewidth = 0.5) +
    geom_text_repel(
        data = top_unique_labels,
        aes(label = gene),
        color = "black",
        size = 3.8,
        max.overlaps = 15,
        box.padding = 0.5,
        point.padding = 0.3,
        force = 2
    ) +
    labs(
        x = expression(log[2] ~ Fold ~ Change),
        y = expression(-log[10] ~ (p[adj])),
        color = "Status"
    ) +
    theme_classic() +
    theme(
        legend.position = "right",
        plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(color = "grey30", size = 11),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12)
    )
ggplot2::ggsave(
    filename = "paper/figures/fig-6-sc-vp.png",
    plot = my_plot,
    width = 19,
    height = 16,
    units = "cm"
)


# Supplement ----
## Functions ----
plot_benchmark <- function(benchmark, title, design_formula) {
    my_plot <- ggplot(
        data = benchmark |> dplyr::filter(method != "disize"),
        mapping = aes(
            x = avg,
            fill = method
        )
    ) +
        geom_ribbon(aes(ymin = q5, ymax = q95), alpha = 0.1) +
        geom_ribbon(aes(ymin = q40, ymax = q60), alpha = 0.2) +
        geom_line(aes(x = avg, y = q50, color = method)) +
        scale_y_log10(breaks = c(1, 5, 10, 25)) +
        scale_x_log10() +
        facet_wrap(
            ncol = 2,
            facets = c("sparsity", "mgt"),
            labeller = labeller(
                sparsity = function(x) paste0("Sparsity: ", x),
                mgt = function(x) paste0("Effect Magnitude: ", x)
            )
        ) +
        labs(
            x = "Average Baseline Expression Level",
            y = "Error Relative To disize",
            color = "Method",
            fill = "Method"
        ) +
        scale_colour_manual(values = my_palette[-1]) +
        scale_fill_manual(values = my_palette[-1]) +
        geom_hline(yintercept = 1.0, linetype = "dashed") +
        theme_classic() +
        theme(
            axis.text = element_text(size = 14),
            axis.title = element_text(size = 16),
            strip.text = element_text(size = 14)
        )

    my_plot
}


## Size Factor Accuracy ----
# Supplement ----
# Scenario 1
benchmark <- read.table(
    "benchmarks/size-factor-accuracy/data/scenario-1.tsv",
    sep = "\t",
    header = TRUE
) |>
    dplyr::mutate(
        method = dplyr::case_match(
            method,
            "mor" ~ "MoR (DESeq/DESeq2)",
            "tmm" ~ "TMM (edgeR)",
            .default = method
        )
    ) |>
    dplyr::filter(type == "relative")
my_plot <- plot_benchmark(benchmark, "A Trivial Setting", ~ (1 | donor_id))
ggplot2::ggsave(
    filename = "paper/figures/suppfig-sfa-1.png",
    plot = my_plot,
    width = 19,
    height = 24,
    units = "cm"
)

# Scenario 1
benchmark <- read.table(
    "benchmarks/size-factor-accuracy/data/scenario-2.tsv",
    sep = "\t",
    header = TRUE
) |>
    dplyr::mutate(
        method = dplyr::case_match(
            method,
            "mor" ~ "MoR (DESeq/DESeq2)",
            "tmm" ~ "TMM (edgeR)",
            .default = method
        )
    ) |>
    dplyr::filter(type == "relative")
my_plot <- plot_benchmark(benchmark, "Comparing Two Conditions", ~cond_id)
ggplot2::ggsave(
    filename = "paper/figures/suppfig-sfa-2.png",
    plot = my_plot,
    width = 19,
    height = 24,
    units = "cm"
)

# Scenario 3
benchmark <- read.table(
    "benchmarks/size-factor-accuracy/data/scenario-3.tsv",
    sep = "\t",
    header = TRUE
) |>
    dplyr::mutate(
        method = dplyr::case_match(
            method,
            "mor" ~ "MoR (DESeq/DESeq2)",
            "tmm" ~ "TMM (edgeR)",
            .default = method
        )
    ) |>
    dplyr::filter(type == "relative")
my_plot <- plot_benchmark(benchmark, "Multifactorial Design", ~ cond_id:sex_id)
ggplot2::ggsave(
    filename = "paper/figures/suppfig-sfa-3.png",
    plot = my_plot,
    width = 19,
    height = 24,
    units = "cm"
)

## Supplement ----
## Read in benchmarking data
benchmark <- read.table(
    here::here("benchmarks/downstream-impact/data/scenario-1.tsv"),
    sep = "\t",
    header = TRUE
) |>
    dplyr::mutate(
        method = dplyr::case_match(
            method,
            "mor" ~ "MoR (DESeq/DESeq2)",
            "tmm" ~ "TMM (edgeR)",
            "ground_truth" ~ "Ground Truth",
            .default = method
        )
    )


# Plot type 1 error
my_plot <- ggplot(
    data = benchmark,
    mapping = aes(x = avg, y = type_1, color = method)
) +
    geom_line() +
    facet_wrap(
        ncol = 1,
        facets = c("sparsity"),
        labeller = labeller(
            sparsity = function(x) paste0("Sparsity: ", x),
            mgt = function(x) paste0("Effect Magnitude: ", x)
        )
    ) +
    labs(
        x = "Average Baseline Expression Level",
        y = "Type I Error",
        color = "Method"
    ) +
    scale_colour_manual(values = my_palette) +
    scale_fill_manual(values = my_palette) +
    theme_classic() +
    theme(
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12)
    )
ggplot2::ggsave(
    filename = "paper/figures/suppfig-di-1.png",
    plot = my_plot,
    width = 18,
    height = 15,
    units = "cm"
)

# Plot type 2 error
my_plot <- ggplot(
    data = benchmark,
    mapping = aes(x = avg, y = type_2, color = method)
) +
    geom_line() +
    facet_wrap(
        ncol = 1,
        facets = c("sparsity"),
        labeller = labeller(
            sparsity = function(x) paste0("Sparsity: ", x),
            mgt = function(x) paste0("Effect Magnitude: ", x)
        )
    ) +
    labs(
        x = "Average Baseline Expression Level",
        y = "Type II Error",
        color = "Method"
    ) +
    scale_colour_manual(values = my_palette) +
    scale_fill_manual(values = my_palette) +
    theme_classic() +
    theme(
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12)
    )
ggplot2::ggsave(
    filename = "paper/figures/suppfig-di-2.png",
    plot = my_plot,
    width = 18,
    height = 15,
    units = "cm"
)
