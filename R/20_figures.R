# =============================================================================
# 20_figures.R: Figures 1-4 and the miRNA sensitivity figure (Supplementary 7)
# Figure numbers follow the manuscript.
# =============================================================================

# Locate the project root even if the working directory was not set:
# when this file is run with source("<path>/R/<script>.R"), move to <path>.
if (!file.exists(file.path("R", "config.R"))) {
  this_file <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NULL)
  if (!is.null(this_file)) setwd(dirname(dirname(this_file)))
}
source(file.path("R", "config.R"))
suppressPackageStartupMessages({
  library(ggplot2)
  library(scales)
  library(patchwork)
})
log_step("Figures")
theme_set(theme_minimal(base_size = 12) +
            theme(panel.grid.minor = element_blank(), plot.title.position = "plot"))

save_fig <- function(p, name, width, height) {
  for (ext in c("png", "tiff", "pdf")) {
    args <- list(filename = file.path(PATHS$figures, paste0(name, ".", ext)), plot = p,
                 width = width, height = height, units = "in", bg = "white")
    if (ext != "pdf") args$dpi <- 600
    if (ext == "tiff") args$compression <- "lzw"
    do.call(ggsave, args)
  }
}

# Each figure is drawn only if its input tables exist, so this script can be
# run after any subset of the analysis steps.
have <- function(...) {
  files <- c(...)
  ok <- all(file.exists(res_path(files)))
  if (!ok) message("  skipped (inputs not yet available): ", paste(files, collapse = ", "))
  ok
}

# ---- Figure 1: prioritized genes ---------------------------------------------
if (have("10_longevitymap/gene_evidence_genetic.csv")) {
  ev <- read_csv(res_path("10_longevitymap", "gene_evidence_genetic.csv"), show_col_types = FALSE) |>
    filter(primary_selected) |>
    arrange(replication_score, desc(gene)) |>          # ties listed alphabetically from the top
    mutate(gene = factor(gene, levels = gene))
  f1 <- ggplot(ev, aes(replication_score, gene)) +
    geom_segment(aes(x = 0, xend = replication_score, yend = gene), colour = "grey85") +
    geom_point(aes(size = n_studies, colour = proportion_sig)) +
    scale_colour_viridis_c(labels = percent, limits = c(0.6, 1), name = "Proportion significant") +
    scale_size_area(max_size = 9, breaks = c(2, 5, 10, 20), name = "Study units") +
    labs(x = "Descriptive replication score", y = NULL,
         caption = "Selection: >= 2 significant study units, >= 2 population labels, proportion significant >= 0.60")
  save_fig(f1, "Figure_1", 8, 7)
}

# ---- Figure 2: threshold sensitivity heatmap ---------------------------------
if (have("11_sensitivity/threshold_scenarios_genetic.csv")) {
  sc <- read_csv(res_path("11_sensitivity", "threshold_scenarios_genetic.csv"), show_col_types = FALSE) |>
    mutate(min_pop = factor(paste("Minimum population labels =", min_pop)),
           min_sig = factor(min_sig, levels = rev(sort(unique(min_sig)))),
           min_prop = percent(min_prop, accuracy = 1),
           primary = min_sig == "2" & min_pop == "Minimum population labels = 2" & min_prop == "60%")
  f2 <- ggplot(sc, aes(min_prop, min_sig, fill = n_genes)) +
    geom_tile(colour = "white", linewidth = 0.6) +
    geom_tile(data = filter(sc, primary), colour = "black", fill = NA, linewidth = 0.9) +
    geom_text(aes(label = n_genes), size = 3.2) +
    facet_wrap(~ min_pop, nrow = 1) +
    scale_fill_gradient(low = "#F2F7F8", high = "#275D71", trans = "sqrt", name = "Selected genes") +
    labs(x = "Minimum proportion of significant study units",
         y = "Minimum number of significant study units") +
    theme(panel.grid = element_blank(), strip.text = element_text(face = "bold"))
  save_fig(f2, "Figure_2", 14, 5)
}

# ---- Figure 3: Open Genes evidence and over-representation --------------------
if (have("13_open_genes/primary_genes_open_genes_records.csv", "13_open_genes/open_genes_overrepresentation.csv")) {
  short <- c("Human longevity association" = "Human longevity",
             "Model-organism lifespan experiment" = "Model-organism lifespan",
             "Age-related molecular change" = "Age-related molecular",
             "Age-related process intervention" = "Age-related process",
             "Gene-regulatory evidence" = "Gene regulation",
             "Open Genes database annotation" = "Database annotation")
  prof <- read_csv(res_path("13_open_genes", "primary_genes_open_genes_records.csv"), show_col_types = FALSE) |>
    mutate(category = factor(short[evidence_type], levels = short))
  ord <- prof |> group_by(gene) |> summarise(t = sum(n)) |> arrange(t) |> pull(gene)
  f3a <- ggplot(mutate(prof, gene = factor(gene, ord)), aes(category, gene, fill = n)) +
    geom_tile(colour = "white") +
    geom_text(aes(label = if_else(n > 0, as.character(n), "")), size = 3) +
    scale_fill_gradient(low = "#F4F6F8", high = "#2F5A6F", name = "Records") +
    labs(x = NULL, y = NULL, title = "A") +
    theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 30, hjust = 1))
  enr <- read_csv(res_path("13_open_genes", "open_genes_overrepresentation.csv"), show_col_types = FALSE) |>
    mutate(category = factor(short[evidence_type], levels = rev(short)),
           label = sprintf("%.2f (FDR %.2f)", adjusted_odds_ratio, adjusted_fdr))
  f3b <- ggplot(enr, aes(adjusted_odds_ratio, category)) +
    geom_vline(xintercept = 1, linetype = 2, colour = "grey50") +
    geom_errorbarh(aes(xmin = adjusted_ci_low, xmax = adjusted_ci_high), height = 0.15, colour = "#2F5A6F") +
    geom_point(size = 3, colour = "#C8664A") +
    geom_text(aes(label = label), vjust = -1, size = 3) +
    scale_x_log10() +
    labs(x = "Adjusted odds ratio (log scale)", y = NULL, title = "B")
  save_fig(f3a + f3b + plot_layout(widths = c(1.2, 1)), "Figure_3", 15, 9)
}

# ---- Figure 4: integrated miRNA target coverage -------------------------------
if (have("15_mirna/ranking_integrated.csv", "15_mirna/interactions_integrated.csv")) {
  primary <- read_lines(res_path("10_longevitymap", "primary_genes.txt"))
  top <- read_csv(res_path("15_mirna", "ranking_integrated.csv"), show_col_types = FALSE) |> slice_head(n = 20)
  inter <- read_csv(res_path("15_mirna", "interactions_integrated.csv"), show_col_types = FALSE)
  bars <- inter |>
    filter(mirna %in% top$mirna) |>
    distinct(mirna, gene) |>
    mutate(group = if_else(gene %in% primary, "Primary LongevityMap candidates", "Additional Open Genes candidates")) |>
    count(mirna, group) |>
    mutate(mirna = factor(mirna, levels = rev(top$mirna)))
  f4 <- ggplot(bars, aes(n, mirna, fill = group)) +
    geom_col(width = 0.7) +
    geom_text(data = top |> mutate(mirna = factor(mirna, levels = rev(top$mirna))),
              aes(x = n_targets, y = mirna, label = n_targets), inherit.aes = FALSE, hjust = -0.2, size = 3.2) +
    scale_fill_manual(values = c("Primary LongevityMap candidates" = "#2C7FB8",
                                 "Additional Open Genes candidates" = "#E66101"),
                      breaks = c("Primary LongevityMap candidates", "Additional Open Genes candidates"),
                      name = NULL) +
    labs(x = "Number of distinct validated target genes", y = NULL) +
    theme(legend.position = "bottom")
  save_fig(f4, "Figure_4", 9, 8)
}

# ---- Supplementary 7: assay-label sensitivity and matched null ----------------
if (have("15_mirna/ranking_integrated_assay.csv", "15_mirna/ranking_primary_assay.csv", "15_mirna/null_comparison_primary.csv", "15_mirna/null_comparison_integrated.csv")) {
  ra <- read_csv(res_path("15_mirna", "ranking_integrated_assay.csv"), show_col_types = FALSE) |> slice_head(n = 10)
  rp <- read_csv(res_path("15_mirna", "ranking_primary_assay.csv"), show_col_types = FALSE)
  s6a_data <- bind_rows(
    transmute(ra, mirna, set = "Integrated candidate set", n_targets),
    tibble(mirna = ra$mirna, set = "Primary genes",
           n_targets = coalesce(rp$n_targets[match(ra$mirna, rp$mirna)], 0L))) |>
    mutate(mirna = factor(mirna, levels = rev(ra$mirna)))
  s6a <- ggplot(s6a_data, aes(n_targets, mirna, fill = set)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.75) +
    scale_fill_manual(values = c("Primary genes" = "#2C7FB8", "Integrated candidate set" = "#E66101"),
                      breaks = c("Primary genes", "Integrated candidate set"), name = NULL) +
    labs(title = "A", x = "Distinct target genes (assay-labelled records)", y = NULL)
  focus <- c("hsa-miR-34a-5p", "hsa-miR-17-5p", "hsa-miR-20a-5p", "hsa-let-7a-5p")
  nulls <- bind_rows(
    read_csv(res_path("15_mirna", "null_comparison_primary.csv"), show_col_types = FALSE),
    read_csv(res_path("15_mirna", "null_comparison_integrated.csv"), show_col_types = FALSE)) |>
    filter(mirna %in% focus) |>
    mutate(mirna = factor(mirna, levels = focus),
           gene_set = factor(gene_set, levels = c("primary", "integrated"),
                             labels = c("Primary LongevityMap genes", "Integrated candidate set")),
           label_y = pmax(observed_targets, null_mean + null_sd))
  s6b <- ggplot(nulls, aes(mirna, observed_targets)) +
    geom_col(fill = "grey70", width = 0.6) +
    geom_pointrange(aes(y = null_mean, ymin = null_mean - null_sd, ymax = null_mean + null_sd), colour = "black") +
    geom_text(aes(y = label_y, label = sprintf("FDR %.2f", empirical_fdr)), vjust = -0.6, size = 3) +
    facet_wrap(~ gene_set, scales = "free_y") +
    labs(title = "B", x = NULL, y = "Distinct target genes",
         caption = "Bars: observed; points: null mean +/- SD (coverage-matched random gene sets)") +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
  save_fig(s6a + s6b + plot_layout(widths = c(1, 1.2)), "Supplementary_7", 13, 7)
}

message("Figures written to ", PATHS$figures)
