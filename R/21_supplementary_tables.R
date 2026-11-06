# =============================================================================
# 21_supplementary_tables.R: export Supplementary 1-6 as Excel workbooks
# Each workbook starts with a ReadMe sheet; sheet names match the manuscript.
# =============================================================================

# Locate the project root even if the working directory was not set:
# when this file is run with source("<path>/R/<script>.R"), move to <path>.
if (!file.exists(file.path("R", "config.R"))) {
  this_file <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NULL)
  if (!is.null(this_file)) setwd(dirname(dirname(this_file)))
}
source(file.path("R", "config.R"))
suppressPackageStartupMessages(library(openxlsx))
log_step("Supplementary tables")

write_book <- function(file, readme, sheets) {
  # Skip a workbook whose input tables have not been produced yet, so this
  # script can be run after any subset of the analysis steps.
  if (any(vapply(sheets, is.null, logical(1)))) {
    message("  skipped (inputs not yet available): ", file)
    return(invisible(NULL))
  }
  wb <- createWorkbook()
  addWorksheet(wb, "ReadMe")
  writeData(wb, "ReadMe", tibble(Item = names(readme), Description = unname(unlist(readme))))
  setColWidths(wb, "ReadMe", cols = 1:2, widths = c(28, 120))
  for (nm in names(sheets)) {
    addWorksheet(wb, nm)
    writeData(wb, nm, sheets[[nm]], withFilter = TRUE)
    freezePane(wb, nm, firstRow = TRUE)
  }
  saveWorkbook(wb, file.path(PATHS$supplementary, file), overwrite = TRUE)
  message("  ", file)
}
rd <- function(...) {
  f <- res_path(...)
  if (file.exists(f)) read_csv(f, show_col_types = FALSE) else NULL
}

write_book("Supplementary_1_LongevityMap_source_review.xlsx",
  list(Content = "Record-level source-review audit of LongevityMap Build 3 (Methods 2.1; Results 3.1).",
       Record_audit = "One row per analysis-ready record (547). gene_s = original LongevityMap assignment; gene_s_corrected = after source review.",
       Source_decisions = "Decisions for records that were reconstructed, reassigned or excluded. Gene counts refer to genes, record counts to records.",
       Manual_source_evidence = "Gene-level evidence derived from the original publications of G461 and G469.",
       LMNA_cohort_overlap = "Records flagged for overlapping NECS samples (Conneely et al.; Sebastiani et al.).",
       Summary = "Counts reported in the manuscript.",
       Script = "R/10_longevitymap.R"),
  list(Record_audit = rd("10_longevitymap", "record_level_audit.csv"),
       Source_decisions = rd("10_longevitymap", "source_review_decisions.csv"),
       Manual_source_evidence = rd("10_longevitymap", "manual_source_evidence.csv"),
       LMNA_cohort_overlap = rd("10_longevitymap", "lmna_cohort_overlap.csv"),
       Summary = rd("10_longevitymap", "summary.csv")))

write_book("Supplementary_2_gene_level_evidence.xlsx",
  list(Content = "Gene-level evidence for the 213-gene source-reviewed universe (Methods 2.2; Results 3.2).",
       Counting = "Study units = one gene per publication; mixed units counted as non-significant. proportion_sig = n_sig_studies / n_studies.",
       Replication_score = "log2(1 + n_studies) x (2 x proportion_sig - 1); descriptive ranking only.",
       Script = "R/10_longevitymap.R"),
  list(Gene_evidence_genetic = rd("10_longevitymap", "gene_evidence_genetic.csv"),
       Gene_evidence_all_modalities = rd("10_longevitymap", "gene_evidence_all_modalities.csv")))

write_book("Supplementary_3_threshold_sensitivity.xlsx",
  list(Content = "Threshold sensitivity over 100 combinations and gene stability (Methods 2.4; Results 3.3).",
       Primary_thresholds = "min significant study units = 2; min population labels = 2; min proportion = 0.60.",
       Stability = "Fraction of the 100 scenarios in which the gene is selected; descriptive only.",
       Criteria_sensitivity = "Selection under alternative definitions: populations counted only among significant study units; population labels collapsed by removing parenthetical qualifiers; both; multi-gene composite records split and included.",
       Downstream_sensitivity = "GO enrichment and miRNA ranking/coverage-matched null repeated for every distinct alternative candidate set (threshold grid, alternative criteria, stable core); one row per distinct gene set. Downstream_GO_significant lists every GO term with FDR <= 0.05 in any of these sets.",
       Script = "R/10_longevitymap.R; R/11_threshold_sensitivity.R; R/16_downstream_sensitivity.R"),
  list(Criteria_sensitivity = rd("10_longevitymap", "criteria_sensitivity_genes.csv"),
       Criteria_summary = rd("10_longevitymap", "criteria_sensitivity_summary.csv"),
       Scenarios_genetic = rd("11_sensitivity", "threshold_scenarios_genetic.csv"),
       Stability_genetic = rd("11_sensitivity", "gene_stability_genetic.csv"),
       Scenarios_all_modalities = rd("11_sensitivity", "threshold_scenarios_all_modalities.csv"),
       Stability_all_modalities = rd("11_sensitivity", "gene_stability_all_modalities.csv"),
       Downstream_sensitivity = rd("16_downstream_sensitivity", "downstream_sensitivity_summary.csv"),
       Downstream_GO_significant = rd("16_downstream_sensitivity", "downstream_GO_significant_terms.csv")))

write_book("Supplementary_5_Open_Genes.xlsx",
  list(Content = "Open Genes integration: overlap with LongevityMap, over-representation and additional candidates (Methods 2.5-2.8; Results 3.5-3.7).",
       Overlap = "Human-longevity records matching LongevityMap by gene and PMID (PMID field or DOI mapped to PMID via Europe PMC) were removed.",
       Direction = "Direction inferred by keyword rules from effect/phenotype/comment text; negative keywords take precedence.",
       Models = "One row per gene of the 213-gene universe; logistic model: candidate ~ evidence + log(1 + LongevityMap studies).",
       Script = "R/02_map_doi_to_pmid.R; R/13_open_genes.R"),
  list(Additional_candidates = rd("13_open_genes", "additional_candidates.csv"),
       Overrepresentation = rd("13_open_genes", "open_genes_overrepresentation.csv"),
       Primary_gene_records = rd("13_open_genes", "primary_genes_open_genes_records.csv"),
       Human_overlap = rd("13_open_genes", "human_records_overlap_with_longevitymap.csv")))

write_book("Supplementary_6_miRNA.xlsx",
  list(Content = "Validated miRNA-target interactions and miRNA rankings (Methods 2.9; Results 3.8).",
       Source = "multiMiR validated tables (miRTarBase, TarBase, miRecords). Versions in data/raw/multimir/.",
       support_type = "TarBase positive/negative = direction of regulation reported by TarBase (negative = repression), not a failed validation.",
       Score = "n_targets x log2(1 + n_databases) x log2(2 + n_publications); descriptive only.",
       Assay_label = paste("Records whose experiment/support fields match:", PARAMS$assay_regex),
       Null = paste0(PARAMS$null_iterations, " coverage-matched random gene sets drawn from the queried LongevityMap + Open Genes universe."),
       Script = "R/14_multimir_query.R; R/15_mirna_analysis.R"),
  list(Interactions_primary = rd("15_mirna", "interactions_primary.csv"),
       Ranking_primary = rd("15_mirna", "ranking_primary.csv"),
       Interactions_integrated = rd("15_mirna", "interactions_integrated.csv"),
       Ranking_integrated = rd("15_mirna", "ranking_integrated.csv"),
       Ranking_primary_assay = rd("15_mirna", "ranking_primary_assay.csv"),
       Ranking_integrated_assay = rd("15_mirna", "ranking_integrated_assay.csv"),
       Null_primary = rd("15_mirna", "null_comparison_primary.csv"),
       Null_integrated = rd("15_mirna", "null_comparison_integrated.csv")))

write_book("Supplementary_4_GO_enrichment.xlsx",
  list(Content = "GO Biological Process over-representation (Methods 2.3; Results 3.4).",
       Primary = "clusterProfiler::enrichGO, propagated annotations, universe = 213 source-reviewed genes, term size 10-100, BH.",
       Sensitivity = "Direct annotations only, and leave-region-out analyses.",
       Script = "R/12_go_enrichment.R"),
  list(GO_propagated = rd("12_go", "GO_BP_propagated_primary.csv"),
       GO_leave_region_out = rd("12_go", "GO_BP_leave_region_out_summary.csv"),
       GO_direct = rd("12_go", "GO_BP_direct_primary.csv"),
       GO_direct_leave_region_out = rd("12_go", "GO_BP_direct_leave_region_out.csv")))
