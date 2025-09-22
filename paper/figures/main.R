library(ggplot2)
library(dplyr)
library(patchwork)

# Functions ----
plot_benchmark <- function(benchmark, title, design_formula) {
    my_plot <- ggplot(
        data = benchmark |> dplyr::filter(comparison != "disize"),
        mapping = aes(
            x = avg,
            fill = comparison
        )
    ) +
        geom_ribbon(aes(ymin = q5, ymax = q95), alpha = 0.1) +
        geom_ribbon(aes(ymin = q40, ymax = q60), alpha = 0.2) +
        geom_line(aes(x = avg, y = q50, color = comparison)) +
        scale_y_log10() +
        facet_wrap(
            ncol = 2,
            facets = c("sparsity", "mgt"),
            labeller = labeller(
                sparsity = function(x) paste0("Sparsity: ", x),
                mgt = function(x) paste0("Effect Magnitude: ", x)
            )
        ) +
        labs(
            title = title,
            subtitle = paste(
                "Design Formula:",
                paste(as.character(design_formula), collapse = " ")
            ),
            x = "Average Baseline Expression Level",
            y = "Relative Error",
            color = "Method",
            fill = "Method"
        ) +
        geom_hline(yintercept = 1.0, linetype = "dashed") +
        theme_classic()

    my_plot
}

# Size Factor Accuracy ----
# Read in benchmarking data
benchmark <- read.table(
    "benchmarks/size-factor-accuracy/data/scenario-2.tsv",
    sep = "\t",
    header = TRUE
) |>
    dplyr::mutate(
        comparison = dplyr::case_match(
            comparison,
            "mor" ~ "MoR (DESeq/DESeq2)",
            "tmm" ~ "TMM (edgeR)",
            .default = comparison
        )
    ) |>
    dplyr::filter(sparsity == 0.2 & mgt == 2)

# Figure 1
my_plot <- ggplot(
    data = benchmark |> dplyr::filter(comparison != "disize"),
    mapping = aes(
        x = avg,
        fill = comparison
    )
) +
    geom_ribbon(aes(ymin = q5, ymax = q95), alpha = 0.1) +
    geom_ribbon(aes(ymin = q40, ymax = q60), alpha = 0.2) +
    geom_line(aes(x = avg, y = q50, color = comparison)) +
    scale_y_log10() +
    labs(
        title = "Comparing Two Conditions",
        subtitle = "Design Formula: ~ cond_id",
        x = "Average Baseline Expression Level",
        y = "Relative Error",
        color = "Method",
        fill = "Method",
        linetype = "Method"
    ) +
    geom_hline(yintercept = 1.0, linetype = "dashed") +
    theme_classic()
ggplot2::ggsave(
    filename = "paper/figures/fig-sfa-1.png",
    plot = my_plot,
    width = 18,
    height = 15,
    units = "cm"
)

# Supplement ----
# Scenario 1
benchmark <- read.table(
    "benchmarks/size-factor-accuracy/data/scenario-1.tsv",
    sep = "\t",
    header = TRUE
) |>
    dplyr::mutate(
        comparison = dplyr::case_match(
            comparison,
            "mor" ~ "MoR (DESeq/DESeq2)",
            "tmm" ~ "TMM (edgeR)",
            .default = comparison
        )
    )
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
        comparison = dplyr::case_match(
            comparison,
            "mor" ~ "MoR (DESeq/DESeq2)",
            "tmm" ~ "TMM (edgeR)",
            .default = comparison
        )
    )
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
        comparison = dplyr::case_match(
            comparison,
            "mor" ~ "MoR (DESeq/DESeq2)",
            "tmm" ~ "TMM (edgeR)",
            .default = comparison
        )
    )
my_plot <- plot_benchmark(benchmark, "Multifactorial Design", ~ cond_id:sex_id)
ggplot2::ggsave(
    filename = "paper/figures/suppfig-sfa-3.png",
    plot = my_plot,
    width = 19,
    height = 24,
    units = "cm"
)


# Downstream Impact ----
# Read in benchmarking data
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
            .default = method
        )
    ) |>
    dplyr::filter(
        method != "ground_truth"
    )

# Figure 2
my_plot <- ((
    ggplot(
        data = benchmark |> dplyr::filter(sparsity == 0.25),
        mapping = aes(x = avg, y = type_1_relative, color = method)
    ) +
        geom_line() +
        labs(
            x = "Average Baseline Expression Level",
            y = "Relative Type 1 Error",
            color = "Method"
        ) +
        geom_hline(yintercept = 1.0, linetype = "dashed") +
        theme_classic()
) / (
    ggplot(
        data = benchmark |> dplyr::filter(sparsity == 0.25),
        mapping = aes(x = avg, y = type_2_relative, color = method)
    ) +
        geom_line() +
        labs(
            x = "Average Baseline Expression Level",
            y = "Relative Type 2 Error",
            color = "Method"
        ) +
        geom_hline(yintercept = 1.0, linetype = "dashed") +
        theme_classic()
)) +
    plot_annotation(
        title = "Comparing Two Conditions",
        subtitle = paste(
            "Design Formula:",
            paste(as.character(~cond_id), collapse = " ")
        )
    )
ggplot2::ggsave(
    filename = "paper/figures/fig-di-1.png",
    plot = my_plot,
    width = 18,
    height = 15,
    units = "cm"
)


## Supplement ----
# Plot type 1 error
my_plot <- ggplot(
    data = benchmark,
    mapping = aes(x = avg, y = type_1_relative, color = method)
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
        title = "Comparing Two Conditions",
        subtitle = paste(
            "Design Formula:",
            paste(as.character(~cond_id), collapse = " ")
        ),
        x = "Average Baseline Expression Level",
        y = "Relative Type 1 Error",
        color = "Method"
    ) +
    geom_hline(yintercept = 1.0, linetype = "dashed") +
    theme_classic()
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
    mapping = aes(x = avg, y = type_2_relative, color = method)
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
        title = "Comparing Two Conditions",
        subtitle = paste(
            "Design Formula:",
            paste(as.character(~cond_id), collapse = " ")
        ),
        x = "Average Baseline Expression Level",
        y = "Relative Type 2 Error",
        color = "Method"
    ) +
    geom_hline(yintercept = 1.0, linetype = "dashed") +
    theme_classic()
ggplot2::ggsave(
    filename = "paper/figures/suppfig-di-2.png",
    plot = my_plot,
    width = 18,
    height = 15,
    units = "cm"
)
