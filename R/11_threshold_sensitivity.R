# =============================================================================
# 11_threshold_sensitivity.R: 100 threshold combinations and gene stability
# (Methods 2.4; Results 3.3; Supplementary 3; Figure 2)
# =============================================================================

# Locate the project root even if the working directory was not set:
# when this file is run with source("<path>/R/<script>.R"), move to <path>.
if (!file.exists(file.path("R", "config.R"))) {
  this_file <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NULL)
  if (!is.null(this_file)) setwd(dirname(dirname(this_file)))
}
source(file.path("R", "config.R"))
in_dir <- res_path("10_longevitymap")
out_dir <- res_path("11_sensitivity")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
log_step("Threshold sensitivity")

grid <- tidyr::crossing(min_sig = PARAMS$grid_min_sig,
                        min_pop = PARAMS$grid_min_pop,
                        min_prop = PARAMS$grid_min_prop)
stopifnot(nrow(grid) == 100L)

run_grid <- function(evidence, label) {
  long <- tidyr::crossing(grid, gene = evidence$gene) |>
    left_join(select(evidence, gene, n_sig_studies, n_populations, proportion_sig), by = "gene") |>
    # a small tolerance avoids floating-point misses at e.g. 0.6 vs 3/5
    mutate(selected = n_sig_studies >= min_sig & n_populations >= min_pop &
             proportion_sig >= min_prop - 1e-9)
  scenarios <- long |>
    group_by(min_sig, min_pop, min_prop) |>
    summarise(n_genes = sum(selected), selected_genes = paste(sort(gene[selected]), collapse = "; "),
              .groups = "drop") |>
    mutate(evidence_table = label)
  stability <- long |>
    group_by(gene) |>
    summarise(selected_scenarios = sum(selected), total_scenarios = n(), .groups = "drop") |>
    mutate(stability = selected_scenarios / total_scenarios, evidence_table = label) |>
    arrange(desc(stability), gene)
  list(scenarios = scenarios, stability = stability)
}

genetic <- run_grid(read_csv(file.path(in_dir, "gene_evidence_genetic.csv"), show_col_types = FALSE),
                    "source-reviewed genetic")
allmod  <- run_grid(read_csv(file.path(in_dir, "gene_evidence_all_modalities.csv"), show_col_types = FALSE),
                    "all modalities")

write_csv(genetic$scenarios, file.path(out_dir, "threshold_scenarios_genetic.csv"))
write_csv(genetic$stability, file.path(out_dir, "gene_stability_genetic.csv"))
write_csv(allmod$scenarios, file.path(out_dir, "threshold_scenarios_all_modalities.csv"))
write_csv(allmod$stability, file.path(out_dir, "gene_stability_all_modalities.csv"))

g <- genetic$scenarios
pick <- function(s, p, q) g$n_genes[g$min_sig == s & g$min_pop == p & abs(g$min_prop - q) < 1e-9]
check_value("Most permissive (1, 1, 0.50)", pick(1, 1, 0.5), 130L)
check_value("Primary (2, 2, 0.60)", pick(2, 2, 0.6), 21L)
check_value("(2, 2, 0.70)", pick(2, 2, 0.7), 12L)
check_value("(2, 2, 0.80)", pick(2, 2, 0.8), 10L)
a <- allmod$scenarios
check_value("All modalities, primary thresholds",
            a$n_genes[a$min_sig == 2 & a$min_pop == 2 & abs(a$min_prop - 0.6) < 1e-9], 22L)
print(head(genetic$stability, 12))
