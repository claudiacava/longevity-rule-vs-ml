# =============================================================================
# config.R: paths, parameters and shared helpers
# Sourced at the top of every analysis script. All paths are relative to the
# project root, so the folder can be moved to any computer without edits.
# =============================================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(tibble)
})

# ---- Project root ------------------------------------------------------------
# The root is the folder that contains run_all.R. When a script is sourced from
# run_all.R the working directory is already the root; when a single script is
# run on its own from the R/ folder, we step up one level.
find_project_root <- function() {
  wd <- normalizePath(getwd(), winslash = "/")
  if (file.exists(file.path(wd, "run_all.R"))) return(wd)
  if (file.exists(file.path(dirname(wd), "run_all.R"))) return(dirname(wd))
  stop("Run the scripts from the project root (the folder containing run_all.R).")
}
ROOT <- find_project_root()

PATHS <- list(
  raw_lm        = file.path(ROOT, "data", "raw", "longevity.csv"),
  raw_og_dir    = file.path(ROOT, "data", "raw", "open_genes"),
  raw_mm_dir    = file.path(ROOT, "data", "raw", "multimir"),
  derived       = file.path(ROOT, "data", "derived"),
  results       = file.path(ROOT, "results"),
  figures       = file.path(ROOT, "figures"),
  supplementary = file.path(ROOT, "supplementary")
)
invisible(lapply(PATHS[c("raw_og_dir", "raw_mm_dir", "derived", "results",
                         "figures", "supplementary")],
                 dir.create, recursive = TRUE, showWarnings = FALSE))

res_path <- function(...) file.path(PATHS$results, ...)

# ---- Analysis parameters (pre-specified) -------------------------------------
PARAMS <- list(
  # Primary LongevityMap selection
  min_sig_studies = 2,
  min_populations = 2,
  min_prop_sig    = 0.60,

  # Threshold-sensitivity grid (5 x 4 x 5 = 100 scenarios)
  grid_min_sig  = 1:5,
  grid_min_pop  = 1:4,
  grid_min_prop = seq(0.50, 0.90, by = 0.10),

  # Gene Ontology (primary analysis: clusterProfiler::enrichGO, propagated
  # annotations). Universe = source-reviewed LongevityMap genes.
  go_min_size   = 10,
  go_max_size   = 100,   # ~half of the 213-gene universe; avoids near-universal terms
  go_fdr        = 0.05,
  go_simplify_cutoff = 0.7,

  # Leave-region-out
  region_gap_bp = 1e6,
  predefined_groups = list(
    APOE_TOMM40_APOC1 = c("APOE", "TOMM40", "APOC1"),
    HSPA_MHC          = c("HSPA1A", "HSPA1B", "HSPA1L", "HLA-DRB1")
  ),

  # Open Genes additional candidates
  og_min_prop_sig = 0.60,

  # miRNA analysis
  assay_regex  = "reporter|luciferase|western[ -]?blot|q(rt)?[- ]?pcr|elisa",
  pmid_truncation_nchar = 240,
  null_iterations = 1000,
  null_strata     = 10,
  null_universe   = "combined",  # "combined" = LongevityMap + Open Genes genes
  multimir_chunk  = 50,

  seed = 20261001
)

# LongevityMap source-review decisions (documented in Supplementary 1)
SOURCE_REVIEW <- list(
  # Passtoors et al. (PMID 23061800, record G461): pathway-level mTOR result.
  passtoors_genes = c("FOXO1", "RPTOR", "EIF4EBP1", "LAMTOR2", "AKT1S1", "PRR5L", "RHOA"),
  # Sebastiani et al. (PMID 22279548, record G469): 281-SNP predictive model;
  # only rs2075650 reached genome-wide significance. rs2075650 lies in an intron
  # of TOMM40 and is assigned to TOMM40 only (as in LongevityMap record 286).
  sebastiani_genes = c("TOMM40"),
  # APOE alleles mis-mapped to genes whose aliases are E2/E3/E4
  apoe_alias_pmids = c("24534555", "24746514"),
  apoe_alias_genes = c("DBT", "DLD", "UBE4B"),
  # Conneely (PMID 22340368) and Sebastiani (PMID 22279548) share NECS samples
  lmna_overlap_pmids = c("22340368", "22279548"),
  lmna_overlap_variant = "rs915179"
)

OG_FILES <- tribble(
  ~file,                                     ~evidence_type,                         ~endpoint,
  "human_longevity_associations.csv",        "Human longevity association",          "research/associations-with-lifespan",
  "model_organism_lifespan_experiments.csv", "Model-organism lifespan experiment",   "research/lifespan-change",
  "age_related_molecular_changes.csv",       "Age-related molecular change",         "research/age-related-changes",
  "age_related_process_interventions.csv",   "Age-related process intervention",     "research/gene-activity-change-impact",
  "gene_regulation.csv",                     "Gene-regulatory evidence",             "research/gene-regulation",
  "open_genes_gene_annotations.csv",         "Open Genes database annotation",       "gene/search"
)
OG_COMPLEMENTARY <- c("Model-organism lifespan experiment", "Age-related molecular change",
                      "Age-related process intervention", "Gene-regulatory evidence")

# ---- Shared helpers ----------------------------------------------------------
norm_text <- function(x) {
  y <- as.character(x)
  y[is.na(y)] <- ""
  str_squish(y)
}

norm_gene <- function(x) {
  y <- str_to_upper(str_remove_all(norm_text(x), "\\s+"))
  y[y %in% c("", "NA", "N/A", "NULL")] <- NA_character_
  y
}

split_genes <- function(x) {
  str_split(str_replace_all(norm_text(x), "[;|]", ","), "\\s*,\\s*|\\s+/\\s+")
}

extract_pmids <- function(x) str_extract_all(norm_text(x), "\\b[0-9]{6,9}\\b")

read_chr_csv <- function(path) {
  read_csv(path, col_types = cols(.default = col_character()),
           na = c("", "NA"), show_col_types = FALSE, progress = FALSE)
}

# Prints observed vs expected values. Expected values are those reported in the
# manuscript; a mismatch is a message, not an error, because some steps are
# deliberately corrected in this pipeline.
check_value <- function(label, observed, expected = NULL) {
  if (is.null(expected)) {
    message(sprintf("  %-55s %s", label, format(observed)))
  } else {
    flag <- if (isTRUE(all.equal(observed, expected))) "OK" else "DIFFERS"
    message(sprintf("  %-55s %s (manuscript: %s) [%s]", label, format(observed),
                    format(expected), flag))
  }
}

log_step <- function(...) message("\n== ", ..., " ==")

# Bioconductor packages (AnnotationDbi, clusterProfiler, multiMiR) mask several
# dplyr verbs. Defining them here (global environment) keeps dplyr's versions
# in every script, whatever was attached before.
select    <- dplyr::select
filter    <- dplyr::filter
rename    <- dplyr::rename
mutate    <- dplyr::mutate
arrange   <- dplyr::arrange
slice     <- dplyr::slice
summarise <- dplyr::summarise
first     <- dplyr::first
