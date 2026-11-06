# =============================================================================
# 10_longevitymap.R: quality control, source review and gene-level evidence
#
# Input : data/raw/longevity.csv (LongevityMap Build 3, 550 records)
# Output: results/10_longevitymap/*  (record audit, gene evidence, candidates)
#
# Rules (Methods 2.1-2.2):
#   * records without id, association or gene symbol are excluded
#   * APOE alleles mis-assigned to DBT/DLD/UBE4B are reassigned to APOE; a record
#     is classified as single- or multi-gene AFTER this correction
#   * G461 (Passtoors) -> seven single-gene expression results, used only in the
#     all-modalities sensitivity analysis
#   * G469 (Sebastiani) -> only rs2075650 retained, assigned to TOMM40 (the gene
#     in which it lies); APOE is not credited with this signal
#   * other multi-gene records are kept for audit but not counted
#   * exact duplicates removed on study, gene, population, association, variant
#   * one gene-study unit per gene and publication; mixed units = non-significant
# =============================================================================

# Locate the project root even if the working directory was not set:
# when this file is run with source("<path>/R/<script>.R"), move to <path>.
if (!file.exists(file.path("R", "config.R"))) {
  this_file <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NULL)
  if (!is.null(this_file)) setwd(dirname(dirname(this_file)))
}
source(file.path("R", "config.R"))
out_dir <- res_path("10_longevitymap")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
log_step("LongevityMap source review")

# ---- 1. Read and quality control -------------------------------------------
raw <- read_chr_csv(PATHS$raw_lm)
raw <- raw[, 1:6]
names(raw) <- c("id", "association", "population", "variant_s", "gene_s", "pub_med")
check_value("LongevityMap records downloaded", nrow(raw), 550L)

dat <- raw |>
  mutate(across(everything(), norm_text),
         association = str_to_lower(association),
         pub_med = str_remove(pub_med, "\\.0$")) |>
  filter(id != "", association != "", gene_s != "") |>
  filter(association %in% c("significant", "non-significant"))

excluded_ids <- setdiff(raw$id, dat$id)
check_value("Analysis-ready records", nrow(dat), 547L)
check_value("Excluded record ids", paste(excluded_ids, collapse = ", "))

# ---- 2. Record-level corrections -------------------------------------------
gene_alias_pattern <- paste0("(^|,)\\s*(", paste(SOURCE_REVIEW$apoe_alias_genes, collapse = "|"),
                             ")\\s*(,|$)")
dat <- dat |>
  mutate(
    study_id = if_else(pub_med == "", paste0("record_", id), paste0("PMID_", pub_med)),
    n_genes_original = lengths(split_genes(gene_s)),
    apoe_alias = pub_med %in% SOURCE_REVIEW$apoe_alias_pmids &
      str_detect(str_to_upper(variant_s), "(^|[,; /])E[234]([,; /]|$)") &
      str_detect(str_to_upper(gene_s), gene_alias_pattern),
    gene_s_corrected = if_else(apoe_alias, "APOE", gene_s),
    n_genes_corrected = lengths(split_genes(gene_s_corrected)),
    multi_gene_record = n_genes_corrected > 1,
    record_type = case_when(
      id == "G469" ~ "model-level predictive model",
      id == "G461" ~ "pathway-level mTOR result",
      apoe_alias & n_genes_original > 1 ~ "APOE alleles mis-mapped to several genes",
      multi_gene_record ~ "multi-gene composite",
      TRUE ~ "single-gene record"
    ),
    source_review_action = case_when(
      id == "G469" ~ "replaced by rs2075650, assigned to TOMM40",
      id == "G461" ~ "replaced by seven single-gene expression results (sensitivity analysis only)",
      id == "G412" ~ "composite GSTT1-TP53 genotype analysis; excluded from gene-level counting",
      multi_gene_record ~ "excluded pending source-specific decomposition",
      apoe_alias ~ "reassigned to APOE and counted as single-gene record",
      TRUE ~ "retained as single-gene record"
    ),
    alias_correction = if_else(apoe_alias, "APOE reassigned from DBT/DLD/UBE4B (allele labels E2/E3/E4)", "none")
  )

check_value("Significant records", sum(dat$association == "significant"), 273L)
check_value("Non-significant records", sum(dat$association == "non-significant"), 274L)
check_value("Distinct PMIDs", n_distinct(dat$pub_med[dat$pub_med != ""]), 270L)
check_value("Distinct population labels", n_distinct(dat$population), 81L)
check_value("Records with >1 original gene symbol", sum(dat$n_genes_original > 1), 66L)
check_value("Records reassigned to APOE", sum(dat$apoe_alias), 5L)
check_value("Multi-gene records after alias correction", sum(dat$multi_gene_record), 63L)

all_symbols <- unique(unlist(split_genes(dat$gene_s_corrected)))
all_symbols <- all_symbols[all_symbols != ""]
check_value("Distinct corrected gene symbols (descriptive audit)", length(all_symbols), 880L)

record_audit <- dat |>
  select(id, study_id, pub_med, association, population, variant_s, gene_s,
         gene_s_corrected, n_genes_original, n_genes_corrected, multi_gene_record,
         record_type, source_review_action, alias_correction)
write_csv(record_audit, file.path(out_dir, "record_level_audit.csv"))

# ---- 3. Manual source evidence ----------------------------------------------
manual_evidence <- bind_rows(
  tibble(record_id = "G461", study_id = "PMID_23061800", pub_med = "23061800",
         association = "significant", population = "Dutch",
         variant_s = "gene-expression single-gene result",
         gene = SOURCE_REVIEW$passtoors_genes,
         evidence_modality = "gene-level expression", shared_locus = NA_character_,
         source_note = "Passtoors et al.: seven genes with corrected single-gene expression differences >= 5%"),
  tibble(record_id = "G469", study_id = "PMID_22279548", pub_med = "22279548",
         association = "significant", population = "American (Caucasian)",
         variant_s = "rs2075650", gene = SOURCE_REVIEW$sebastiani_genes,
         evidence_modality = "genetic association", shared_locus = "TOMM40-APOE-APOC1 region",
         source_note = paste("Sebastiani et al.: rs2075650 (TOMM40 intron) was the only genome-wide-significant",
                             "SNP of the 281-SNP model; it duplicates LongevityMap record 286 and is removed",
                             "as an exact duplicate"))
)
write_csv(manual_evidence, file.path(out_dir, "manual_source_evidence.csv"))

# ---- 4. Gene-level evidence --------------------------------------------------
single_gene <- dat |>
  filter(!multi_gene_record, !id %in% c("G461", "G469")) |>
  transmute(record_id = id, study_id, pub_med, association, population, variant_s,
            gene = str_trim(gene_s_corrected), evidence_modality = "genetic association",
            shared_locus = NA_character_)

gene_records_all <- bind_rows(single_gene, select(manual_evidence, -source_note)) |>
  distinct(study_id, gene, population, association, variant_s, .keep_all = TRUE) |>
  mutate(
    overlap_group = if_else(
      pub_med %in% SOURCE_REVIEW$lmna_overlap_pmids &
        str_detect(variant_s, paste0("(^|,)", SOURCE_REVIEW$lmna_overlap_variant, "(,|$)")),
      "NECS_LMNA_rs915179", NA_character_),
    study_unit = coalesce(overlap_group, study_id)
  )

summarise_genes <- function(records) {
  units <- records |>
    group_by(gene, study_unit) |>
    summarise(has_sig = any(association == "significant"),
              has_non = any(association == "non-significant"), .groups = "drop") |>
    mutate(study_class = if_else(has_non, "non-significant", "significant"),
           mixed = has_sig & has_non)
  pops <- records |>
    group_by(gene) |>
    summarise(n_populations = n_distinct(population),
              populations = paste(sort(unique(population)), collapse = "; "), .groups = "drop")
  units |>
    group_by(gene) |>
    summarise(n_studies = n(),
              n_sig_studies = sum(study_class == "significant"),
              n_nonsig_studies = sum(study_class == "non-significant"),
              n_mixed_studies = sum(mixed), .groups = "drop") |>
    left_join(pops, by = "gene") |>
    mutate(proportion_sig = n_sig_studies / n_studies,
           replication_score = log2(1 + n_studies) * (2 * proportion_sig - 1),
           primary_selected = n_sig_studies >= PARAMS$min_sig_studies &
             n_populations >= PARAMS$min_populations &
             proportion_sig >= PARAMS$min_prop_sig) |>
    arrange(desc(primary_selected), desc(replication_score), gene) |>
    relocate(populations, .after = last_col())
}

gene_records_genetic <- filter(gene_records_all, evidence_modality == "genetic association")
gene_evidence <- summarise_genes(gene_records_genetic)
gene_evidence_all_mod <- summarise_genes(gene_records_all)
primary_genes <- gene_evidence$gene[gene_evidence$primary_selected]

check_value("Source-reviewed genetic gene universe", nrow(gene_evidence), 213L)
check_value("All-modalities gene universe", nrow(gene_evidence_all_mod), 219L)
check_value("Primary candidate genes", length(primary_genes), 21L)
check_value("All-modalities candidates at primary thresholds",
            sum(gene_evidence_all_mod$primary_selected), 22L)
apoe <- filter(gene_evidence, gene == "APOE")
check_value("APOE study units / significant / populations",
            paste(apoe$n_studies, apoe$n_sig_studies, apoe$n_populations, sep = "/"))

write_csv(gene_records_all, file.path(out_dir, "gene_level_records_all_modalities.csv"))
write_csv(gene_evidence, file.path(out_dir, "gene_evidence_genetic.csv"))
write_csv(gene_evidence_all_mod, file.path(out_dir, "gene_evidence_all_modalities.csv"))
write_lines(primary_genes, file.path(out_dir, "primary_genes.txt"))

# ---- 4b. Sensitivity of the selection criteria --------------------------------
# (i)   population labels counted only among significant study units
# (ii)  population labels collapsed by removing parenthetical qualifiers
#       (e.g. "Italian (Southern)" -> "Italian")
# (iii) (i) and (ii) combined
# (iv)  multi-gene composite records included, split into their listed genes
country_label <- function(x) str_squish(str_remove_all(x, "\\s*\\([^)]*\\)"))

select_genes <- function(records, pop_fun = identity, sig_pops_only = FALSE) {
  records <- mutate(records, pop = pop_fun(population))
  units <- records |>
    group_by(gene, study_unit) |>
    summarise(sig_unit = !any(association == "non-significant"), .groups = "drop")
  pops <- records |>
    left_join(units, by = c("gene", "study_unit")) |>
    filter(!sig_pops_only | sig_unit) |>
    group_by(gene) |>
    summarise(n_pop_criterion = n_distinct(pop), .groups = "drop")
  units |>
    group_by(gene) |>
    summarise(n_studies = n(), n_sig_studies = sum(sig_unit), .groups = "drop") |>
    left_join(pops, by = "gene") |>
    mutate(n_pop_criterion = coalesce(n_pop_criterion, 0L),
           proportion_sig = n_sig_studies / n_studies,
           selected = n_sig_studies >= PARAMS$min_sig_studies &
             n_pop_criterion >= PARAMS$min_populations &
             proportion_sig >= PARAMS$min_prop_sig - 1e-9)
}

composite_records <- dat |>
  filter(multi_gene_record, !id %in% c("G461", "G469")) |>
  transmute(record_id = id, study_id, pub_med, association, population, variant_s,
            gene = split_genes(gene_s_corrected)) |>
  unnest_longer(gene) |>
  filter(gene != "") |>
  mutate(study_unit = study_id)
records_with_composites <- bind_rows(gene_records_genetic, composite_records) |>
  distinct(study_id, gene, population, association, variant_s, .keep_all = TRUE) |>
  mutate(study_unit = coalesce(overlap_group, study_id))

criteria_runs <- list(
  primary                         = select_genes(gene_records_genetic),
  significant_populations_only    = select_genes(gene_records_genetic, sig_pops_only = TRUE),
  country_level_labels            = select_genes(gene_records_genetic, pop_fun = country_label),
  significant_and_country_level   = select_genes(gene_records_genetic, country_label, TRUE),
  composite_records_included      = select_genes(records_with_composites)
)
criteria_long <- imap_dfr(criteria_runs, ~ mutate(.x, analysis = .y))
criteria_wide <- criteria_long |>
  filter(gene %in% unique(criteria_long$gene[criteria_long$selected])) |>
  select(gene, analysis, selected) |>
  pivot_wider(names_from = analysis, values_from = selected, values_fill = FALSE) |>
  mutate(n_analyses_selected = rowSums(across(-gene))) |>
  arrange(desc(primary), desc(n_analyses_selected), gene)
criteria_summary <- criteria_long |>
  group_by(analysis) |>
  summarise(genes_in_universe = n(), genes_selected = sum(selected),
            selected_genes = paste(sort(gene[selected]), collapse = "; "), .groups = "drop") |>
  mutate(lost_vs_primary = map_chr(selected_genes, ~ paste(setdiff(primary_genes, str_split(.x, "; ")[[1]]), collapse = "; ")),
         gained_vs_primary = map_chr(selected_genes, ~ paste(setdiff(str_split(.x, "; ")[[1]], primary_genes), collapse = "; ")))
write_csv(criteria_long, file.path(out_dir, "criteria_sensitivity_long.csv"))
write_csv(criteria_wide, file.path(out_dir, "criteria_sensitivity_genes.csv"))
write_csv(criteria_summary, file.path(out_dir, "criteria_sensitivity_summary.csv"))
for (i in seq_len(nrow(criteria_summary))) {
  check_value(paste0("Criteria sensitivity: ", criteria_summary$analysis[i]),
              paste0(criteria_summary$genes_selected[i], " genes; lost: ",
                     criteria_summary$lost_vs_primary[i], "; gained: ", criteria_summary$gained_vs_primary[i]))
}
check_value("Genes selected in all five analyses", sum(criteria_wide$n_analyses_selected == 5))

# ---- 5. Audit tables for Supplementary 1 -------------------------------------
composite <- dat |> filter(multi_gene_record, !id %in% c("G461", "G469"))
composite_genes <- unique(unlist(split_genes(composite$gene_s_corrected)))
source_decisions <- tibble(
  decision_unit = c("G461 (Passtoors et al.)", "G469 (Sebastiani et al.)",
                    "APOE allele records (PMIDs 24534555, 24746514)",
                    "Other multi-gene composite records"),
  n_records = c(1L, 1L, sum(dat$apoe_alias), nrow(composite)),
  original_assignment = c(paste(dat$n_genes_original[dat$id == "G461"], "mTOR-pathway genes"),
                          paste(dat$n_genes_original[dat$id == "G469"],
                                "genes from a 281-SNP predictive model"),
                          "DBT, DLD and/or UBE4B (aliases E2, E3, E4)",
                          paste(sum(lengths(split_genes(composite$gene_s_corrected))),
                                "gene assignments;", length(composite_genes), "distinct genes")),
  genes_counted_primary_genetic = c(0L, 1L, 1L, 0L),
  genes_counted_all_modalities = c(7L, 1L, 1L, 0L),
  decision = c(
    "Pathway-level result; seven single-gene expression results kept for the all-modalities sensitivity analysis only",
    "Model membership not treated as gene-level evidence; only rs2075650 retained and assigned to TOMM40 (already present as record 286 from the same publication, so no additional study unit is created)",
    "Reassigned to APOE; counted as single-gene APOE records",
    paste0("Not decomposed; retained for audit only. ",
           length(setdiff(composite_genes, gene_evidence$gene)),
           " of these genes have no single-gene record and are outside the gene universe")
  )
)
write_csv(source_decisions, file.path(out_dir, "source_review_decisions.csv"))

lmna_overlap <- dat |>
  filter(pub_med %in% SOURCE_REVIEW$lmna_overlap_pmids,
         str_detect(variant_s, SOURCE_REVIEW$lmna_overlap_variant) | id == "G469") |>
  transmute(id, pub_med, association, population, gene_s, source_review_action,
            note = paste("Conneely et al. report that Sebastiani et al. included most of the NECS",
                         "cases and controls. Both records were excluded or replaced during source",
                         "review, so no LMNA evidence enters the gene-level counts."))
write_csv(lmna_overlap, file.path(out_dir, "lmna_cohort_overlap.csv"))

summary_tbl <- tibble(
  metric = c("records_downloaded", "records_valid", "records_excluded", "significant", "non_significant",
             "distinct_pmids", "population_labels", "records_multi_gene_original",
             "records_apoe_alias_corrected", "records_multi_gene_after_correction",
             "corrected_gene_symbols", "gene_universe_genetic", "gene_universe_all_modalities",
             "primary_candidates", "all_modalities_candidates"),
  value = c(nrow(raw), nrow(dat), length(excluded_ids), sum(dat$association == "significant"),
            sum(dat$association == "non-significant"), n_distinct(dat$pub_med[dat$pub_med != ""]),
            n_distinct(dat$population), sum(dat$n_genes_original > 1), sum(dat$apoe_alias),
            sum(dat$multi_gene_record), length(all_symbols), nrow(gene_evidence),
            nrow(gene_evidence_all_mod), length(primary_genes), sum(gene_evidence_all_mod$primary_selected))
)
write_csv(summary_tbl, file.path(out_dir, "summary.csv"))
message("Primary genes: ", paste(primary_genes, collapse = ", "))
