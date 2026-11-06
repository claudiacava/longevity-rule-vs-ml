# =============================================================================
# run_all.R: runs the complete pipeline in order.
# Open R (or RStudio) in this folder and run:  source("run_all.R")
# Each step can also be run on its own, in the same order, after the previous
# steps have written their results.
# =============================================================================

steps <- c(
  "R/00_setup.R",                 # packages + sessionInfo
  "R/01_download_open_genes.R",   # uses data/raw/open_genes if present
  "R/02_map_doi_to_pmid.R",       # Europe PMC; cached in data/derived
  "R/10_longevitymap.R",          # Supplementary 1-2, primary genes
  "R/11_threshold_sensitivity.R", # Supplementary 3
  "R/12_go_enrichment.R",         # Supplementary 4
  "R/13_open_genes.R",            # Supplementary 5, additional candidates
  "R/14_multimir_query.R",        # multiMiR; cached in data/raw/multimir
  "R/15_mirna_analysis.R",        # Supplementary 6
  "R/16_downstream_sensitivity.R",# GO + miRNA across alternative candidate sets (Supplementary 3)
  "R/20_figures.R",               # Figures 1-4, Supplementary 7
  "R/21_supplementary_tables.R"   # Excel workbooks
)

dir.create("results", showWarnings = FALSE)
log_file <- file.path("results", paste0("run_log_", format(Sys.time(), "%Y%m%d_%H%M"), ".txt"))

# Messages are shown on screen and also copied to the log file, so the
# observed-vs-manuscript checks can be reviewed after the run.
for (s in steps) {
  t0 <- Sys.time()
  cat("\n########## ", s, " ##########\n", file = log_file, append = TRUE)
  message("\n########## ", s, " ##########")
  withCallingHandlers(
    source(s, local = new.env(), echo = FALSE),
    message = function(m) cat(conditionMessage(m), file = log_file, append = TRUE),
    warning = function(w) cat("WARNING: ", conditionMessage(w), "\n", file = log_file, append = TRUE)
  )
  message("finished in ", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " min")
}
message("\nPipeline completed. Checks and warnings saved in ", log_file)
