# =============================================================================
# run_temporal.R: prospective temporal validation and rule-versus-ML comparison.
# Run after run_all.R (scripts 10 and 12 produce the candidate genes and the
# gene coordinates used here). Open R in this folder and run:
#   source("run_temporal.R")
# Scripts 33-35 download large third-party files and train the models
# (about 1 hour in total).
# =============================================================================

steps <- c(
  "R/30_download_temporal.R",   # GWAS Catalog, NCBI gene2pubmed, Open Genes (current), Europe PMC dates
  "R/30b_fix_doi_map.R",        # corrects three mistyped Open Genes DOIs (Europe PMC)
  "R/31_temporal_validation.R", # outcomes O1 (Open Genes) and O2 (GWAS Catalog) after the cut-off
  "R/32_ml_feasibility.R",      # training labels for the back-tested M1 models
  "R/33_download_ml_features.R",# GO 2017, MSigDB 6.0, BioGRID 3.4.149, GTEx v8, GenAge
  "R/34_build_ml_features.R",   # feature matrices
  "R/35_ml_models.R",           # M1 and M2 models; scores of the 213 candidate genes
  "R/36_rule_vs_ml.R"           # rule vs ML on the post-cut-off outcomes
)

dir.create("results", showWarnings = FALSE)
log_file <- file.path("results", paste0("run_log_temporal_", format(Sys.time(), "%Y%m%d_%H%M"), ".txt"))
run_step <- function(s) {
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
if (exists("ML_PERMUTE")) rm(ML_PERMUTE)
for (s in steps) run_step(s)

# Negative control: models retrained with permuted labels (results/35_ml_models_permuted)
ML_PERMUTE <- TRUE
run_step("R/35_ml_models.R")
rm(ML_PERMUTE)
message("\nDone. Messages saved in ", log_file)
