# =============================================================================
# 13_open_genes.R: Open Genes integration (Methods 2.5-2.8; Results 3.5-3.7)
#
#  1. Import the six Open Genes datasets as separate evidence categories.
#  2. Remove Open Genes human-longevity records that duplicate LongevityMap
#     records (same gene and same publication). Publications are matched by
#     PMID, using the DOI -> PMID map from 02_map_doi_to_pmid.R because Open
#     Genes reports DOIs for almost all human records.
#  3. Gene-level over-representation of each category among the primary genes
#     versus the rest of the source-reviewed universe: one-sided Fisher test and
#     logistic regression adjusted for log(1 + number of LongevityMap studies).
#     Each model uses one row per gene (213 rows).
#  4. Additional exploratory candidates.
# =============================================================================

# Locate the project root even if the working directory was not set:
# when this file is run with source("<path>/R/<script>.R"), move to <path>.
if (!file.exists(file.path("R", "config.R"))) {
  this_file <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NULL)
  if (!is.null(this_file)) setwd(dirname(dirname(this_file)))
}
source(file.path("R", "config.R"))
suppressPackageStartupMessages(library(broom))

lm_dir <- res_path("10_longevitymap")
out_dir <- res_path("13_open_genes")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
log_step("Open Genes integration")

# ---- 1. Import ---------------------------------------------------------------
first_col <- function(dat, candidates) {
  hit <- intersect(candidates, names(dat))
  if (length(hit)) norm_text(dat[[hit[1]]]) else rep("", nrow(dat))
}
norm_sig <- function(x) {
  y <- str_to_lower(norm_text(x))
  case_when(
    y %in% c("1", "true", "yes", "significant") ~ "significant",
    y %in% c("0", "false", "no", "non-significant", "non significant", "nonsignificant") ~ "non-significant",
    TRUE ~ "unknown"
  )
}
# Direction is inferred from the effect/phenotype/comment text with keyword
# rules (negative keywords take precedence). This is a heuristic and is
# reported as such in Methods; the full classification is exported for audit.
NEG_PATTERN <- "non[- ]?associated|decreas|shorter|deteriorat|negative|reduc|lower|risk"
POS_PATTERN <- "longevity|increase|extend|improv|protect|longer|beneficial|positive"

read_og <- function(file, evidence_type) {
  d <- read_chr_csv(file.path(PATHS$raw_og_dir, file))
  effect <- first_col(d, c("longevityEffect", "result", "interventionResultForLifespan", "changeType",
                           "regulationType", "interventionImproves", "interventionDeteriorates"))
  phenotype <- first_col(d, c("longevityEffect", "vitalProcess", "result", "comment"))
  comment <- first_col(d, c("comment", "description"))
  text <- str_to_lower(paste(effect, phenotype, comment, sep = " | "))
  tibble(
    evidence_type = evidence_type,
    record_id = paste(file, seq_len(nrow(d)), sep = "::"),
    gene = norm_gene(first_col(d, c("geneSymbol", "symbol"))),
    pmid_raw = first_col(d, c("pmid")),
    doi = str_remove(str_to_lower(first_col(d, c("doi"))), "^https?://(dx\\.)?doi\\.org/"),
    variant = first_col(d, c("polymorphismId", "associatedAllele")),
    population = first_col(d, c("ethnicity", "sample", "organismLine")),
    effect_text = effect,
    significance_raw = first_col(d, c("significance", "pValue", "lMeanChangeStatSignificance")),
    significance = norm_sig(first_col(d, c("significance", "pValue", "lMeanChangeStatSignificance"))),
    direction = case_when(str_detect(text, NEG_PATTERN) ~ "negative",
                          str_detect(text, POS_PATTERN) ~ "positive",
                          TRUE ~ "unclear")
  ) |>
    filter(!is.na(gene))
}

og <- map2_dfr(OG_FILES$file, OG_FILES$evidence_type, read_og)
check_value("Open Genes usable records", nrow(og), 13055L)
check_value("Open Genes distinct genes", n_distinct(og$gene), 2404L)
check_value("Human-longevity records", sum(og$evidence_type == "Human longevity association"), 1460L)

# ---- 2. Publication overlap with LongevityMap --------------------------------
doi_map_file <- file.path(PATHS$derived, "doi_to_pmid_map.csv")
if (file.exists(doi_map_file)) {
  doi_map <- read_chr_csv(doi_map_file) |> filter(!is.na(pmid)) |> distinct(doi, pmid_from_doi = pmid)
} else {
  warning("DOI -> PMID map not found: run R/02_map_doi_to_pmid.R. Overlap will use PMIDs only.")
  doi_map <- tibble(doi = character(), pmid_from_doi = character())
}

human <- og |>
  filter(evidence_type == "Human longevity association") |>
  left_join(doi_map, by = "doi") |>
  mutate(pmid_list = map2(extract_pmids(pmid_raw), pmid_from_doi,
                          ~ unique(c(.x, if (!is.na(.y)) .y))))

lm_pairs <- read_csv(file.path(lm_dir, "record_level_audit.csv"),
                     col_types = cols(.default = col_character())) |>
  transmute(lm_id = id, pmid = str_extract(pub_med, "[0-9]{6,9}"), gene = split_genes(gene_s_corrected)) |>
  unnest_longer(gene) |>
  mutate(gene = norm_gene(gene)) |>
  filter(!is.na(gene), !is.na(pmid)) |>
  distinct()

overlap <- human |>
  select(record_id, gene, pmid_raw, doi, pmid_from_doi, pmid_list) |>
  unnest_longer(pmid_list, values_to = "pmid", keep_empty = TRUE) |>
  left_join(lm_pairs, by = c("gene", "pmid"), relationship = "many-to-many") |>
  group_by(record_id, gene, pmid_raw, doi, pmid_from_doi) |>
  summarise(overlap_any = any(!is.na(lm_id)),
            overlap_via = case_when(any(!is.na(lm_id) & pmid %in% extract_pmids(first(pmid_raw))[[1]]) ~ "PMID field",
                                    any(!is.na(lm_id)) ~ "DOI mapped to PMID",
                                    TRUE ~ "none"),
            matched_lm_ids = paste(sort(unique(na.omit(lm_id))), collapse = "; "),
            .groups = "drop")
write_csv(overlap, file.path(out_dir, "human_records_overlap_with_longevitymap.csv"))

dup_ids <- overlap$record_id[overlap$overlap_any]
check_value("Human records with a PMID", sum(human$pmid_raw != ""), 48L)
check_value("Human records with a PMID after DOI mapping", sum(lengths(human$pmid_list) > 0))
check_value("Overlapping human records (PMID field only)", sum(overlap$overlap_via == "PMID field"), 5L)
check_value("Overlapping human records (PMID + DOI mapping)", length(dup_ids))

og <- og |>
  mutate(retained = !(evidence_type == "Human longevity association" & record_id %in% dup_ids))
write_csv(og, file.path(out_dir, "open_genes_records_reconciled.csv"))

# ---- 3. Over-representation among the primary genes --------------------------
lm_evidence <- read_csv(file.path(lm_dir, "gene_evidence_genetic.csv"), show_col_types = FALSE) |>
  mutate(gene = norm_gene(gene))
primary <- lm_evidence$gene[lm_evidence$primary_selected]
check_value("Universe genes present in Open Genes", sum(lm_evidence$gene %in% og$gene), 134L)
check_value("Combined LongevityMap + Open Genes universe",
            length(union(lm_evidence$gene, og$gene)), 2483L)

present <- og |> filter(retained) |> distinct(gene, evidence_type)
gene_cat <- tidyr::crossing(select(lm_evidence, gene, n_studies), evidence_type = OG_FILES$evidence_type) |>
  left_join(mutate(present, evidence_present = TRUE), by = c("gene", "evidence_type")) |>
  mutate(evidence_present = coalesce(evidence_present, FALSE),
         primary_candidate = gene %in% primary)
write_csv(gene_cat, file.path(out_dir, "gene_by_category_presence.csv"))

enrich <- map_dfr(OG_FILES$evidence_type, function(cat) {
  d <- filter(gene_cat, evidence_type == cat)          # one row per gene
  tab <- table(factor(d$primary_candidate, c(TRUE, FALSE)), factor(d$evidence_present, c(TRUE, FALSE)))
  ft <- fisher.test(tab, alternative = "greater")
  fit <- glm(primary_candidate ~ evidence_present + log1p(n_studies), data = d, family = binomial())
  est <- tidy(fit, exponentiate = TRUE, conf.int = TRUE) |> filter(term == "evidence_presentTRUE")
  tibble(evidence_type = cat,
         candidates_with_evidence = tab["TRUE", "TRUE"], candidates_total = sum(tab["TRUE", ]),
         others_with_evidence = tab["FALSE", "TRUE"], others_total = sum(tab["FALSE", ]),
         fisher_odds_ratio = unname(ft$estimate), fisher_p = ft$p.value,
         adjusted_odds_ratio = est$estimate, adjusted_ci_low = est$conf.low,
         adjusted_ci_high = est$conf.high, adjusted_p = est$p.value)
}) |>
  mutate(fisher_fdr = p.adjust(fisher_p, "BH"), adjusted_fdr = p.adjust(adjusted_p, "BH"))
write_csv(enrich, file.path(out_dir, "open_genes_overrepresentation.csv"))
print(mutate(enrich, across(where(is.numeric), ~ signif(.x, 3))), width = Inf)

# ---- 4. Additional exploratory candidates ------------------------------------
human_summary <- og |>
  filter(retained, evidence_type == "Human longevity association") |>
  group_by(gene) |>
  summarise(og_human_records = n(),
            og_human_significant = sum(significance == "significant"),
            og_human_nonsignificant = sum(significance == "non-significant"),
            og_human_unknown_significance = sum(significance == "unknown"),
            og_human_positive = sum(direction == "positive"),
            og_human_negative = sum(direction == "negative"),
            og_human_unclear_direction = sum(direction == "unclear"),
            og_human_positive_significant = sum(direction == "positive" & significance == "significant"),
            og_human_negative_significant = sum(direction == "negative" & significance == "significant"),
            .groups = "drop") |>
  mutate(og_human_classified = og_human_significant + og_human_nonsignificant,
         og_human_significant_proportion = if_else(og_human_classified > 0,
                                                   og_human_significant / og_human_classified, NA_real_))

complementary <- og |>
  filter(retained, evidence_type %in% OG_COMPLEMENTARY) |>
  count(gene, evidence_type) |>
  mutate(evidence_type = paste0("og_", str_replace_all(evidence_type, "[^A-Za-z0-9]+", "_"))) |>
  pivot_wider(names_from = evidence_type, values_from = n, values_fill = 0)
comp_cols <- setdiff(names(complementary), "gene")

candidates <- human_summary |>
  full_join(complementary, by = "gene") |>
  mutate(across(where(is.numeric) & !og_human_significant_proportion, ~ replace_na(.x, 0))) |>
  left_join(select(lm_evidence, gene, n_studies, n_sig_studies, n_nonsig_studies,
                   proportion_sig, n_populations), by = "gene") |>
  mutate(in_lm_gene_universe = !is.na(n_studies),
         in_primary_set = gene %in% primary,
         n_complementary_categories = rowSums(across(all_of(comp_cols), ~ .x > 0)),
         additional_candidate = !in_primary_set &
           og_human_positive_significant >= 1 &
           coalesce(og_human_significant_proportion >= PARAMS$og_min_prop_sig, FALSE) &
           n_complementary_categories >= 1,
         lm_status = case_when(!in_lm_gene_universe ~ "not in LongevityMap gene universe",
                               proportion_sig >= PARAMS$min_prop_sig ~ "LongevityMap proportion >= 0.60",
                               TRUE ~ "LongevityMap proportion < 0.60")) |>
  arrange(desc(additional_candidate), desc(og_human_positive_significant),
          desc(n_complementary_categories), gene)
write_csv(candidates, file.path(out_dir, "candidate_screen_all_genes.csv"))
additional <- filter(candidates, additional_candidate)
write_csv(additional, file.path(out_dir, "additional_candidates.csv"))

check_value("Additional candidates", nrow(additional), 69L)
check_value("  ... in LongevityMap universe", sum(additional$in_lm_gene_universe), 26L)
check_value("  ... with proportion = 1", sum(additional$og_human_significant_proportion == 1), 57L)
check_value("  ... with LongevityMap proportion < 0.60",
            sum(additional$lm_status == "LongevityMap proportion < 0.60"))
check_value("Complementary categories (1/2/3/4)",
            paste(tabulate(additional$n_complementary_categories, 4), collapse = "/"), "18/19/16/16")

primary_profile <- og |>
  filter(retained, gene %in% primary) |>
  count(gene, evidence_type) |>
  complete(gene = primary, evidence_type = OG_FILES$evidence_type, fill = list(n = 0))
write_csv(primary_profile, file.path(out_dir, "primary_genes_open_genes_records.csv"))
check_value("Retained Open Genes records on primary genes", sum(primary_profile$n), 797L)
six <- primary_profile |> group_by(gene) |> summarise(k = sum(n > 0)) |> filter(k == 6)
check_value("Primary genes present in all six categories", paste(sort(six$gene), collapse = ", "))
write_lines(sort(union(primary, additional$gene)), file.path(out_dir, "integrated_genes.txt"))
