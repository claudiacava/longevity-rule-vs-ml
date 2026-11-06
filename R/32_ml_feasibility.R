# =============================================================================
# 32_ml_feasibility.R: can a machine-learning model be trained on the same
# LongevityMap evidence without using the future?
#
# Analysis plan: 3 November 2026.
# The model on the same evidence is trained with a pseudo-cut-off P:
#   features = LongevityMap gene-study units published before P
#   label    = confirmation published in [P, 2017-06-01) (Open Genes, GWAS)
# and then applied to all LongevityMap evidence to predict confirmation after
# 2017-06-01 (the outcomes of script 31).
#
# This script only counts. For each candidate pseudo-cut-off it reports
#   1. genes with LongevityMap evidence before P, and genes the published rule
#      would have selected at P
#   2. confirmed genes in the training window [P, 2017-06-01): Open Genes (O1),
#      GWAS Catalog +/- 500 kb (O2), either
#   3. how many training labels come from publications already in LongevityMap
#   4. coverage of GO Biological Process annotations (current org.Hs.eg.db)
#      and its correlation with research attention
# Decision rule written in the plan: about 15 confirmed genes or more in the
# training window are needed to train the model.
#
# Inputs : results of scripts 10, 12, 30; data/raw/temporal/<snapshot>/
# Output : results/32_ml_feasibility/*.csv
# =============================================================================

if (!file.exists(file.path("R", "config.R"))) {
  this_file <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NULL)
  if (!is.null(this_file)) setwd(dirname(dirname(this_file)))
}
source(file.path("R", "config.R"))

FS <- list(
  pseudo_cutoffs = as.Date(c("2010-01-01", "2011-01-01", "2012-01-01", "2013-01-01")),
  cutoff         = as.Date("2017-06-01"),
  gws            = 5e-8,
  window_bp      = 5e5,
  min_events     = 15
)
out_dir <- res_path("32_ml_feasibility")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
snapshots <- sort(list.dirs(file.path(ROOT, "data", "raw", "temporal"), recursive = FALSE, full.names = TRUE))
snap_dir <- tail(snapshots[grepl("[0-9]{4}-[0-9]{2}-[0-9]{2}$", snapshots)], 1)
stopifnot(length(snap_dir) == 1)
log_step("ML feasibility (snapshot ", basename(snap_dir), ")")

# ---- Identifiers and dates (as in script 31) ---------------------------------------
upd <- read_chr_csv(res_path("12_go", "symbol_updates.csv"))
entrez <- read_chr_csv(res_path("12_go", "symbol_to_entrez.csv"))
official <- function(g) { o <- upd$symbol[match(g, upd$original)]; str_to_upper(coalesce(o, g)) }
coords <- read_csv(res_path("12_go", "gene_coordinates_regions.csv"), show_col_types = FALSE,
                   col_types = cols(SYMBOL = "c", chr = "c", .default = col_guess()))
dates <- read_chr_csv(file.path(PATHS$derived, "pmid_publication_dates.csv"))
pub_date <- function(p) {
  d <- dates[match(p, dates$pmid), ]
  as.Date(coalesce(d$first_publication_date, ifelse(is.na(d$pub_year), NA_character_, paste0(d$pub_year, "-01-01"))), optional = TRUE)
}

# ---- 1. LongevityMap evidence before each pseudo-cut-off ------------------------------
rec <- read_chr_csv(res_path("10_longevitymap", "gene_level_records_all_modalities.csv")) |>
  filter(evidence_modality == "genetic association") |>
  mutate(date = pub_date(pub_med))
check_value("LongevityMap genetic records without a date", sum(is.na(rec$date)), 0L)
universe <- sort(unique(rec$gene))
check_value("Genes in the universe", length(universe), 213L)

# Same definitions as summarise_genes() in script 10
evidence_at <- function(r) {
  units <- r |> group_by(gene, study_unit) |>
    summarise(has_non = any(association == "non-significant"), .groups = "drop")
  pops <- r |> group_by(gene) |> summarise(n_populations = n_distinct(population), .groups = "drop")
  units |> group_by(gene) |>
    summarise(n_studies = n(), n_sig_studies = sum(!has_non), .groups = "drop") |>
    left_join(pops, by = "gene") |>
    mutate(proportion_sig = n_sig_studies / n_studies,
           selected = n_sig_studies >= PARAMS$min_sig_studies & n_populations >= PARAMS$min_populations &
             proportion_sig >= PARAMS$min_prop_sig - 1e-9)
}
full <- evidence_at(rec)
check_value("Rule on all LongevityMap evidence: selected genes", sum(full$selected), 21L)

# ---- 2. Outcomes in a time window ------------------------------------------------------
norm_doi <- function(x) norm_text(x) |> str_to_lower() |> str_remove("^https?://(dx\\.)?doi\\.org/")
doi_map <- read_chr_csv(file.path(PATHS$derived, "doi_to_pmid_map_temporal.csv"))
og <- read_chr_csv(file.path(PATHS$raw_og_dir, "human_longevity_associations.csv")) |>
  mutate(pmid_n = coalesce(str_extract(norm_text(pmid), "[0-9]{6,9}"), doi_map$pmid[match(norm_doi(doi), doi_map$doi)]),
         date = pub_date(pmid_n), gene_key = str_to_upper(norm_text(geneSymbol)),
         sig = significance %in% c("Significant", "1"), nonsig = significance %in% c("Non-significant", "0"))
og_window <- function(from, to) {
  d <- filter(og, !is.na(date), date >= from, date < to, !is.na(pmid_n))
  units <- d |> group_by(gene_key, pmid_n) |> summarise(s = any(sig), n = any(nonsig), .groups = "drop") |>
    filter(s, !n)
  tibble(gene = universe, key = official(universe)) |>
    mutate(og = key %in% units$gene_key,
           og_pmids = vapply(key, function(k) paste(sort(unique(units$pmid_n[units$gene_key == k])), collapse = ";"), "")) |>
    select(-key)
}

traits <- read_chr_csv(file.path(PATHS$derived, "temporal_trait_selection.csv"))
gw <- readr::read_tsv(file.path(snap_dir, "gwas_catalog_longevity_candidate_associations.tsv"),
                      col_types = cols(.default = col_character()), quote = "", progress = FALSE, na = c("", "NA")) |>
  left_join(select(traits, mapped_trait, disease_trait, include),
            by = c("MAPPED_TRAIT" = "mapped_trait", "DISEASE/TRAIT" = "disease_trait")) |>
  mutate(include = as.logical(include) %in% TRUE, date = as.Date(DATE, optional = TRUE),
         p = suppressWarnings(as.numeric(`P-VALUE`)), chr = CHR_ID, pos = suppressWarnings(as.numeric(CHR_POS)))
gene_win <- coords |> transmute(gene = SYMBOL, chr = as.character(chr),
                                lo = start - FS$window_bp, hi = end + FS$window_bp)
gwas_window <- function(from, to) {
  d <- filter(gw, include, !is.na(date), date >= from, date < to, !is.na(p), p <= FS$gws)
  mapped <- unique(str_to_upper(unlist(split_genes(str_replace_all(coalesce(d$MAPPED_GENE, ""), "\\s+-\\s+|\\s+x\\s+", ",")))))
  hits <- filter(d, !is.na(pos))
  win_hit <- vapply(seq_len(nrow(gene_win)), function(i)
    any(hits$chr == gene_win$chr[i] & hits$pos >= gene_win$lo[i] & hits$pos <= gene_win$hi[i]), logical(1))
  tibble(gene = universe) |>
    mutate(gwas = official(gene) %in% mapped | official(gene) %in% str_to_upper(gene_win$gene[win_hit]))  # official symbols (fix 5 Nov 2026)
}

# Check against script 31 (outcomes after the cut-off)
chk <- left_join(og_window(FS$cutoff, as.Date("2100-01-01")), gwas_window(FS$cutoff, as.Date("2100-01-01")), by = "gene")
check_value("Confirmed after 2017-06-01: Open Genes / GWAS", paste(sum(chk$og), sum(chk$gwas), sep = " / "), "27 / 31")

# ---- 3. One row per pseudo-cut-off -----------------------------------------------------
lm_pmids <- unique(rec$pub_med)
gene_rows <- list()
summary_rows <- lapply(seq_along(FS$pseudo_cutoffs), function(i) {
  P <- FS$pseudo_cutoffs[i]
  ev <- evidence_at(filter(rec, date < P))
  lab <- left_join(og_window(P, FS$cutoff), gwas_window(P, FS$cutoff), by = "gene") |>
    mutate(either = og | gwas)
  g <- tibble(gene = universe, pseudo_cutoff = P) |>
    left_join(select(ev, gene, n_studies_before = n_studies, n_sig_before = n_sig_studies, selected_at_P = selected), by = "gene") |>
    mutate(across(c(n_studies_before, n_sig_before), ~ coalesce(.x, 0L)), selected_at_P = coalesce(selected_at_P, FALSE)) |>
    left_join(lab, by = "gene")
  # training labels whose only confirming publications are already in LongevityMap
  # are circular (LongevityMap predicting LongevityMap): the label used for
  # training is "independent" = GWAS, or Open Genes from a publication not in LongevityMap
  og_only_lm <- g$og & vapply(strsplit(g$og_pmids, ";"), function(p) length(p) > 0 && all(p %in% lm_pmids), logical(1))
  g <- mutate(g, independent = gwas | (og & !og_only_lm))
  gene_rows[[i]] <<- g
  with_ev <- g$n_studies_before > 0
  tibble(pseudo_cutoff = P,
         genes_with_evidence_before = sum(with_ev),
         genes_without_evidence_before = sum(!with_ev),
         rule_selected_at_P = sum(g$selected_at_P),
         confirmed_open_genes = sum(g$og),
         confirmed_gwas = sum(g$gwas),
         confirmed_either = sum(g$either),
         confirmed_either_among_genes_with_evidence = sum(g$either & with_ev),
         open_genes_labels_only_from_longevitymap_pmids = sum(og_only_lm),
         independent_among_genes_with_evidence = sum(g$independent & with_ev),
         independent_gwas_among_genes_with_evidence = sum(g$gwas & with_ev),
         rule_at_P_independent_hits = sum(g$independent & g$selected_at_P),
         rule_at_P_precision_either = if (sum(g$selected_at_P) > 0) sum(g$either & g$selected_at_P) / sum(g$selected_at_P) else NA_real_,
         trainable = sum(g$independent & with_ev) >= FS$min_events)
})
fs_sum <- bind_rows(summary_rows)
write_csv(fs_sum, file.path(out_dir, "feasibility_summary.csv"), na = "")
write_csv(bind_rows(gene_rows), file.path(out_dir, "feasibility_gene_level.csv"), na = "")

# ---- 4. GO coverage (current annotations) and research attention ------------------------
go_cov <- tryCatch({
  suppressPackageStartupMessages(library(org.Hs.eg.db))
  ids <- entrez$ENTREZID[match(official(universe), str_to_upper(entrez$SYMBOL))]
  go <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys = unique(na.omit(ids)), keytype = "ENTREZID",
                                               columns = c("GOALL", "ONTOLOGYALL")))
  n_terms <- go |> filter(ONTOLOGYALL == "BP") |> distinct(ENTREZID, GOALL) |> count(ENTREZID, name = "go_bp_terms")
  out <- tibble(gene = universe, entrez = ids) |> left_join(n_terms, by = c("entrez" = "ENTREZID")) |>
    mutate(go_bp_terms = coalesce(go_bp_terms, 0L))
  g31 <- read_csv(res_path("30_temporal", "gene_outcomes.csv"), show_col_types = FALSE)
  out <- left_join(out, select(g31, gene, pubs_pre), by = "gene")
  write_csv(out, file.path(out_dir, "go_coverage.csv"), na = "")
  tibble(genes_without_go_bp = sum(out$go_bp_terms == 0), median_go_bp_terms = median(out$go_bp_terms),
         spearman_go_terms_vs_pubs_pre = suppressWarnings(cor(out$go_bp_terms, out$pubs_pre, method = "spearman", use = "complete.obs")))
}, error = function(e) { message("  GO coverage skipped: ", conditionMessage(e)); NULL })
if (!is.null(go_cov)) write_csv(go_cov, file.path(out_dir, "go_coverage_summary.csv"), na = "")

# ---- Report --------------------------------------------------------------------------
log_step("Training window [pseudo-cut-off, 2017-06-01)")
print(as.data.frame(select(fs_sum, pseudo_cutoff, genes_with_evidence_before, rule_selected_at_P,
                           confirmed_either_among_genes_with_evidence, open_genes_labels_only_from_longevitymap_pmids,
                           independent_among_genes_with_evidence, rule_at_P_independent_hits, trainable)), row.names = FALSE)
message("Decision (plan, 3 Nov 2026): primary pseudo-cut-off 2013-01-01, sensitivity 2011-01-01; ",
        "training label = independent confirmation only")
if (!is.null(go_cov)) { log_step("GO Biological Process (current annotations)"); print(as.data.frame(go_cov)) }
message("\nWritten: ", out_dir)
