# =============================================================================
# 31_temporal_validation.R: were the genes prioritised from LongevityMap
# confirmed by evidence published afterwards?
#
# Pre-specified plan (agreed 2 November 2026)
#   Outcomes (per gene of the LongevityMap universe)
#     O1 Open Genes: >= 1 significant gene-study unit in human longevity
#        studies published on/after the cut-off (unit = gene x publication;
#        a unit with any non-significant record is non-significant)
#     O2 GWAS Catalog: genome-wide significant association (P <= 5e-8) for a
#        selected longevity trait, published on/after the cut-off, within
#        +/- 500 kb of the gene (GRCh38); genes without coordinates use the
#        catalogue's mapped gene. Also analysed at region level (1-Mb regions
#        of the article), so APOE-TOMM40-APOC1 counts once.
#   Cut-off: LongevityMap Build 3 release (1 June 2017); sensitivity: day after
#   the latest LongevityMap publication (6 August 2015).
#   A  Primary: 21 selected vs 192 other genes; Fisher test and logistic
#      regression (Firth, penalised) adjusted for research attention (log PubMed publications per
#      gene before and after the cut-off, NCBI gene2pubmed); Open Genes also
#      conditional on the gene being re-tested.
#   B  Selection rules: 100 thresholds, alternative definitions, LongevityMap
#      without source review, simple baselines; precision and recall, compared
#      with random gene sets of equal size matched on research attention.
#   C  Effect of source review (reviewed vs raw LongevityMap).
#   D  Sensitivity: earlier cut-off, own vs parental longevity, healthspan
#      traits added, P <= 1e-5, variants inside the gene only, sex (Open Genes).
#
# Inputs : results of scripts 10-12, data/raw/temporal/<snapshot>/ (script 30)
# Output : results/30_temporal/*.csv
# =============================================================================

if (!file.exists(file.path("R", "config.R"))) {
  this_file <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NULL)
  if (!is.null(this_file)) setwd(dirname(dirname(this_file)))
}
source(file.path("R", "config.R"))

TV <- list(
  cutoffs   = c(release = as.Date("2017-06-01"), last_lm_publication = as.Date("2015-08-06")),
  gws       = 5e-8,
  suggestive = 1e-5,
  window_bp = 5e5,
  n_null    = 1000,
  n_strata  = 5,
  min_events = 3   # minimum confirmed genes per group to fit the adjusted model
)
if (!requireNamespace("logistf", quietly = TRUE)) install.packages("logistf")
out_dir <- res_path("30_temporal")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
snapshots <- sort(list.dirs(file.path(ROOT, "data", "raw", "temporal"), recursive = FALSE, full.names = TRUE))
snap_dir <- tail(snapshots[grepl("[0-9]{4}-[0-9]{2}-[0-9]{2}$", snapshots)], 1)
stopifnot(length(snap_dir) == 1)
log_step("Temporal validation (snapshot ", basename(snap_dir), ")")
set.seed(PARAMS$seed)

# ---- Gene identifiers ----------------------------------------------------------------
# LongevityMap symbols -> current official symbols (script 12) -> Entrez IDs
upd <- read_chr_csv(res_path("12_go", "symbol_updates.csv"))
entrez <- read_chr_csv(res_path("12_go", "symbol_to_entrez.csv"))
official <- function(g) { o <- upd$symbol[match(g, upd$original)]; str_to_upper(coalesce(o, g)) }
coords <- read_csv(res_path("12_go", "gene_coordinates_regions.csv"), show_col_types = FALSE,
                   col_types = cols(SYMBOL = "c", chr = "c", .default = col_guess()))

# ---- LongevityMap evidence: reviewed, alternative definitions, raw ---------------------
ev <- read_csv(res_path("10_longevitymap", "gene_evidence_genetic.csv"), show_col_types = FALSE)
crit <- read_csv(res_path("10_longevitymap", "criteria_sensitivity_long.csv"), show_col_types = FALSE)
select_at <- function(d, s, p, q, popcol = "n_populations")
  d$gene[d$n_sig_studies >= s & d[[popcol]] >= p & d$proportion_sig >= q - 1e-9]
primary <- select_at(ev, PARAMS$min_sig_studies, PARAMS$min_populations, PARAMS$min_prop_sig)
check_value("Primary genes", length(primary), 21L)

# LongevityMap without source review: every record with one gene symbol, as
# written, no corrections; same units and rule
raw <- read_chr_csv(PATHS$raw_lm)[, 1:6]
names(raw) <- c("id", "association", "population", "variant_s", "gene_s", "pub_med")
raw <- raw |>
  mutate(across(everything(), norm_text), association = str_to_lower(association),
         pub_med = str_remove(pub_med, "\\.0$")) |>
  filter(id != "", gene_s != "", association %in% c("significant", "non-significant")) |>
  mutate(study_id = if_else(pub_med == "", paste0("record_", id), paste0("PMID_", pub_med)),
         n_genes = lengths(split_genes(gene_s))) |>
  filter(n_genes == 1) |>
  mutate(gene = str_trim(gene_s)) |>
  distinct(study_id, gene, population, association, variant_s, .keep_all = TRUE)
raw_ev <- raw |>
  group_by(gene, study_id) |> summarise(sig = !any(association == "non-significant"), .groups = "drop") |>
  group_by(gene) |> summarise(n_studies = n(), n_sig_studies = sum(sig), .groups = "drop") |>
  left_join(raw |> group_by(gene) |> summarise(n_populations = n_distinct(population), .groups = "drop"), by = "gene") |>
  mutate(proportion_sig = n_sig_studies / n_studies)
raw_selected <- select_at(raw_ev, PARAMS$min_sig_studies, PARAMS$min_populations, PARAMS$min_prop_sig)
check_value("LongevityMap without source review: universe / selected",
            paste(nrow(raw_ev), length(raw_selected), sep = " / "), "214 / 22")

all_genes <- unique(c(ev$gene, crit$gene, raw_ev$gene))

# ---- Publication dates ---------------------------------------------------------------
dates <- read_chr_csv(file.path(PATHS$derived, "pmid_publication_dates.csv"))
pub_date <- function(p) {
  d <- dates[match(p, dates$pmid), ]
  as.Date(coalesce(d$first_publication_date, ifelse(is.na(d$pub_year), NA_character_, paste0(d$pub_year, "-01-01"))), optional = TRUE)
}

# ---- Outcome O1: Open Genes ----------------------------------------------------------
norm_doi <- function(x) norm_text(x) |> str_to_lower() |> str_remove("^https?://(dx\\.)?doi\\.org/")
doi_map <- read_chr_csv(file.path(PATHS$derived, "doi_to_pmid_map_temporal.csv"))
og <- read_chr_csv(file.path(PATHS$raw_og_dir, "human_longevity_associations.csv")) |>
  mutate(pmid_n = coalesce(str_extract(norm_text(pmid), "[0-9]{6,9}"), doi_map$pmid[match(norm_doi(doi), doi_map$doi)]),
         date = pub_date(pmid_n), gene_key = str_to_upper(norm_text(geneSymbol)),
         sig = significance %in% c("Significant", "1"), nonsig = significance %in% c("Non-significant", "0"),
         sex = str_to_lower(norm_text(sex)))
og_outcome <- function(cut, sex_filter = NULL) {
  d <- filter(og, !is.na(date), date >= cut, !is.na(pmid_n))
  if (!is.null(sex_filter)) d <- filter(d, sex %in% sex_filter)
  units <- d |> group_by(gene_key, pmid_n) |> summarise(s = any(sig), n = any(nonsig), .groups = "drop")
  tibble(gene = all_genes, key = official(all_genes)) |>
    mutate(retested = key %in% d$gene_key,
           confirmed = key %in% units$gene_key[units$s & !units$n]) |>
    select(-key)
}

# ---- Outcome O2: GWAS Catalog --------------------------------------------------------
traits <- read_chr_csv(file.path(PATHS$derived, "temporal_trait_selection.csv"))
gw <- readr::read_tsv(file.path(snap_dir, "gwas_catalog_longevity_candidate_associations.tsv"),
                      col_types = cols(.default = col_character()), quote = "", progress = FALSE, na = c("", "NA")) |>
  left_join(select(traits, mapped_trait, disease_trait, include),
            by = c("MAPPED_TRAIT" = "mapped_trait", "DISEASE/TRAIT" = "disease_trait")) |>
  mutate(include = as.logical(include) %in% TRUE,
         date = as.Date(DATE, optional = TRUE),
         p = suppressWarnings(as.numeric(`P-VALUE`)),
         chr = CHR_ID, pos = suppressWarnings(as.numeric(CHR_POS)),
         parental = str_detect(MAPPED_TRAIT, "parental") & !str_detect(MAPPED_TRAIT, "life span determination"),
         healthspan = str_detect(MAPPED_TRAIT, "healthspan"),
         intergenic = str_detect(coalesce(MAPPED_GENE, ""), "\\s-\\s"))
gene_win <- coords |> transmute(gene = SYMBOL, key = official(SYMBOL), chr = as.character(chr),
                                lo = start - TV$window_bp, hi = end + TV$window_bp)
gwas_outcome <- function(cut, pmax = TV$gws, traitset = "main", in_gene_only = FALSE) {
  d <- filter(gw, !is.na(date), date >= cut, !is.na(p), p <= pmax)
  d <- switch(traitset,
              main = filter(d, include),
              own = filter(d, include, !parental),
              parental = filter(d, include, parental),
              with_healthspan = filter(d, include | healthspan))
  if (in_gene_only) d <- filter(d, !intergenic)
  mapped <- unique(str_to_upper(unlist(split_genes(str_replace_all(coalesce(d$MAPPED_GENE, ""), "\\s+-\\s+|\\s+x\\s+", ",")))))
  hits <- filter(d, !is.na(pos))
  win_hit <- vapply(seq_len(nrow(gene_win)), function(i)
    any(hits$chr == gene_win$chr[i] & hits$pos >= gene_win$lo[i] & hits$pos <= gene_win$hi[i]), logical(1))
  win_genes <- if (in_gene_only) character() else gene_win$gene[win_hit]
  tibble(gene = all_genes, key = official(all_genes)) |>
    mutate(confirmed_mapped = key %in% mapped,
           confirmed = confirmed_mapped | key %in% str_to_upper(win_genes)) |>  # official symbols (fix 5 Nov 2026)
    select(-key)
}

# ---- Research attention (PubMed publications per gene) ---------------------------------
g2p_file <- file.path(snap_dir, "gene2pubmed_human.tsv")
pmid_cut <- function(cut) {
  # PMIDs are assigned in time order; the threshold that best separates the
  # dated PMIDs of this study around the cut-off is used for gene2pubmed
  cal <- dates |> filter(status == "found") |> mutate(p = as.numeric(pmid), d = pub_date(pmid)) |> filter(!is.na(d))
  cand <- sort(unique(cal$p))
  acc <- vapply(cand, function(t) mean((cal$p >= t) == (cal$d >= cut)), numeric(1))
  cand[which.max(acc)]
}
attention_for <- function(cut) {
  base <- tibble(gene = all_genes, key = official(all_genes),
                 entrez = entrez$ENTREZID[match(official(all_genes), str_to_upper(entrez$SYMBOL))])
  if (!file.exists(g2p_file)) {
    message("  gene2pubmed not available: research attention = LongevityMap study units")
    return(base |> left_join(select(ev, gene, n_studies), by = "gene") |>
             transmute(gene, pubs_pre = coalesce(n_studies, 1L), pubs_post = NA_real_, source = "LongevityMap units"))
  }
  g2p <- readr::read_tsv(g2p_file, col_types = cols(.default = col_character()), progress = FALSE)
  t <- pmid_cut(cut)
  counts <- g2p |> filter(gene_id %in% base$entrez) |>
    group_by(gene_id) |> summarise(pubs_pre = sum(as.numeric(pmid) < t), pubs_post = sum(as.numeric(pmid) >= t), .groups = "drop")
  base |> left_join(counts, by = c("entrez" = "gene_id")) |>
    transmute(gene, pubs_pre = coalesce(pubs_pre, 0), pubs_post = coalesce(pubs_post, 0),
              source = if_else(is.na(entrez), "no Entrez ID", paste0("gene2pubmed, PMID threshold ", t)))
}

# ---- A. Primary comparison ----------------------------------------------------------
compare_sets <- function(outcome, att, universe, selected, label, conditional_on = NULL) {
  d <- tibble(gene = universe) |>
    left_join(outcome, by = "gene") |> left_join(att, by = "gene") |>
    mutate(selected = gene %in% selected)
  if (!is.null(conditional_on)) d <- filter(d, .data[[conditional_on]])
  a <- sum(d$selected & d$confirmed); b <- sum(d$selected & !d$confirmed)
  c <- sum(!d$selected & d$confirmed); dd <- sum(!d$selected & !d$confirmed)
  ft <- fisher.test(matrix(c(a, b, c, dd), 2), alternative = "greater")
  # Adjusted model: Firth's penalised logistic regression (stable with rare
  # outcomes and separation); not fitted with fewer than TV$min_events
  # confirmed genes in either group, where any estimate would be unreliable
  f <- if (all(is.na(d$pubs_post))) confirmed ~ selected + log1p(pubs_pre) else confirmed ~ selected + log1p(pubs_pre) + log1p(pubs_post)
  adj_note <- if (min(a, c) < TV$min_events) paste("not estimated: fewer than", TV$min_events, "confirmed genes in a group") else "Firth logistic regression"
  adj <- if (min(a, c) < TV$min_events) c(or = NA, lo = NA, hi = NA, p = NA) else tryCatch({
    m <- logistf::logistf(f, data = mutate(d, selected = as.numeric(selected), confirmed = as.numeric(confirmed)))
    k <- which(names(coef(m)) == "selected")
    c(or = exp(unname(coef(m)[k])), lo = exp(unname(m$ci.lower[k])), hi = exp(unname(m$ci.upper[k])), p = unname(m$prob[k]))
  }, error = function(e) { adj_note <<- paste("model failed:", conditionMessage(e)); c(or = NA, lo = NA, hi = NA, p = NA) })
  tibble(analysis = label, n_selected = a + b, confirmed_selected = a, n_other = c + dd, confirmed_other = c,
         pct_selected = a / (a + b), pct_other = c / (c + dd),
         fisher_or = unname(ft$estimate), fisher_p = ft$p.value,
         adjusted_or = adj[["or"]], adjusted_ci_low = adj[["lo"]], adjusted_ci_high = adj[["hi"]], adjusted_p = adj[["p"]],
         adjusted_model = adj_note,
         confirmed_genes_selected = paste(sort(d$gene[d$selected & d$confirmed]), collapse = "; "))
}
region_level <- function(outcome, universe, selected, label) {
  reg <- tibble(gene = universe) |>
    mutate(region = coalesce(coords$region[match(official(gene), str_to_upper(coords$SYMBOL))], paste0("unplaced_", gene))) |>
    left_join(outcome, by = "gene") |>
    group_by(region) |> summarise(selected = any(gene %in% selected), confirmed = any(confirmed),
                                  genes = paste(gene, collapse = "; "), .groups = "drop")
  a <- sum(reg$selected & reg$confirmed); b <- sum(reg$selected & !reg$confirmed)
  c <- sum(!reg$selected & reg$confirmed); dd <- sum(!reg$selected & !reg$confirmed)
  ft <- fisher.test(matrix(c(a, b, c, dd), 2), alternative = "greater")
  tibble(analysis = label, n_selected = a + b, confirmed_selected = a, n_other = c + dd, confirmed_other = c,
         pct_selected = a / (a + b), pct_other = c / (c + dd), fisher_or = unname(ft$estimate), fisher_p = ft$p.value,
         confirmed_genes_selected = paste(reg$genes[reg$selected & reg$confirmed], collapse = " | "))
}

cut0 <- TV$cutoffs[["release"]]
att0 <- attention_for(cut0)
o1 <- og_outcome(cut0)
o2 <- gwas_outcome(cut0)
gene_out <- tibble(gene = ev$gene) |>
  mutate(selected = gene %in% primary, stable_core = gene %in% Reduce(intersect, split(crit$gene[crit$selected], crit$analysis[crit$selected]))) |>
  left_join(rename(o1, og_retested = retested, og_confirmed = confirmed), by = "gene") |>
  left_join(rename(o2, gwas_confirmed_mapped_gene = confirmed_mapped, gwas_confirmed = confirmed), by = "gene") |>
  left_join(att0, by = "gene") |>
  mutate(region = coords$region[match(official(gene), str_to_upper(coords$SYMBOL))])
write_csv(gene_out, file.path(out_dir, "gene_outcomes.csv"))

U <- ev$gene
primary_tbl <- bind_rows(
  compare_sets(o1, att0, U, primary, "O1 Open Genes: significant unit"),
  compare_sets(o1, att0, U, primary, "O1 Open Genes: significant unit, among re-tested genes", "retested"),
  compare_sets(o2, att0, U, primary, "O2 GWAS: P <= 5e-8 within 500 kb"),
  region_level(o2, U, primary, "O2 GWAS: P <= 5e-8, region level")
)
write_csv(primary_tbl, file.path(out_dir, "primary_comparison.csv"))
log_step("A. Primary comparison (cut-off ", cut0, ")")
for (i in seq_len(nrow(primary_tbl))) with(primary_tbl[i, ], check_value(analysis, sprintf(
  "%d/%d vs %d/%d; Fisher P = %.2g; adjusted OR = %s (%s-%s), P = %s", confirmed_selected, n_selected, confirmed_other, n_other,
  fisher_p, signif(adjusted_or, 3), signif(adjusted_ci_low, 3), signif(adjusted_ci_high, 3), signif(adjusted_p, 2))))
check_value("Research attention measure", unique(att0$source)[1])

# ---- B. Selection rules -------------------------------------------------------------
log_step("B. Selection rules")
rules <- list()
for (s in PARAMS$grid_min_sig) for (p in PARAMS$grid_min_pop) for (q in PARAMS$grid_min_prop)
  rules[[sprintf("grid: sig>=%d, pop>=%d, prop>=%.2f", s, p, q)]] <- list(universe = U, selected = select_at(ev, s, p, q), family = "threshold grid")
for (a in unique(crit$analysis)) {
  d <- filter(crit, analysis == a)
  rules[[paste("definition:", a)]] <- list(universe = d$gene, selected = d$gene[d$selected], family = "alternative definition")
}
rules[["LongevityMap without source review"]] <- list(universe = raw_ev$gene, selected = raw_selected, family = "no source review")
rules[["baseline: >= 1 significant unit"]] <- list(universe = U, selected = ev$gene[ev$n_sig_studies >= 1], family = "baseline")
rules[["baseline: all genes"]] <- list(universe = U, selected = U, family = "baseline")

strata_of <- function(att, genes) {
  x <- log1p(att$pubs_pre[match(genes, att$gene)] + coalesce(att$pubs_post[match(genes, att$gene)], 0))
  x[is.na(x)] <- 0
  dplyr::ntile(x, TV$n_strata)
}
null_test <- function(universe, selected, outcome, att) {
  conf <- outcome$confirmed[match(universe, outcome$gene)] %in% TRUE
  st <- strata_of(att, universe)
  need <- table(factor(st[universe %in% selected], levels = 1:TV$n_strata))
  obs <- sum(conf[universe %in% selected])
  sims <- replicate(TV$n_null, {
    pick <- unlist(lapply(names(need), function(k) {
      pool <- which(st == as.integer(k)); n <- min(need[[k]], length(pool))
      if (n > 0) pool[sample.int(length(pool), n)] else integer()
    }))
    sum(conf[pick])
  })
  c(null_mean = mean(sims), p_emp = (1 + sum(sims >= obs)) / (TV$n_null + 1))
}
perf <- function(outcome, outcome_label, att) {
  bind_rows(lapply(names(rules), function(nm) {
    r <- rules[[nm]]; sel <- intersect(r$selected, r$universe)
    conf <- outcome$confirmed[match(r$universe, outcome$gene)] %in% TRUE
    n_sel <- length(sel); hit <- sum(conf[r$universe %in% sel])
    nt <- if (n_sel > 0 && n_sel < length(r$universe) && all(r$universe %in% att$gene)) null_test(r$universe, sel, outcome, att) else c(null_mean = NA, p_emp = NA)
    tibble(outcome = outcome_label, rule = nm, family = r$family, universe = length(r$universe), n_selected = n_sel,
           confirmed_selected = hit, precision = if (n_sel) hit / n_sel else NA_real_,
           recall = if (sum(conf)) hit / sum(conf) else NA_real_, prevalence = mean(conf),
           enrichment = if (n_sel) (hit / n_sel) / mean(conf) else NA_real_,
           null_mean_confirmed = nt[["null_mean"]], null_p = nt[["p_emp"]])
  }))
}
rule_tbl <- bind_rows(perf(o1, "O1 Open Genes", att0), perf(o2, "O2 GWAS", att0))
write_csv(rule_tbl, file.path(out_dir, "selection_rules_performance.csv"))
for (o in unique(rule_tbl$outcome)) {
  d <- filter(rule_tbl, outcome == o)
  pr <- filter(d, rule == "grid: sig>=2, pop>=2, prop>=0.60")
  check_value(paste(o, "- primary rule precision / recall / null P"),
              sprintf("%.2f / %.2f / %.3f", pr$precision, pr$recall, pr$null_p))
  check_value(paste(o, "- grid rules with null P < 0.05"),
              sprintf("%d of %d (with >= 3 genes)", sum(d$family == "threshold grid" & d$n_selected >= 3 & d$null_p < 0.05, na.rm = TRUE),
                      sum(d$family == "threshold grid" & d$n_selected >= 3)))
}

# ---- C. Effect of source review -----------------------------------------------------------
log_step("C. Source review")
only_raw <- setdiff(raw_selected, primary); only_rev <- setdiff(primary, raw_selected)
o1all <- og_outcome(cut0); o2all <- gwas_outcome(cut0)
review_tbl <- tibble(gene = c(only_raw, only_rev),
                     selected_in = c(rep("raw LongevityMap only", length(only_raw)), rep("source-reviewed only", length(only_rev)))) |>
  left_join(rename(o1all, og_retested = retested, og_confirmed = confirmed), by = "gene") |>
  left_join(select(rename(o2all, gwas_confirmed = confirmed), gene, gwas_confirmed), by = "gene")
write_csv(review_tbl, file.path(out_dir, "source_review_effect.csv"))
check_value("Selected only without review", paste(only_raw, collapse = ", "))
check_value("Selected only after review", paste(only_rev, collapse = ", "))

# ---- D. Sensitivity analyses --------------------------------------------------------------
log_step("D. Sensitivity analyses")
cut1 <- TV$cutoffs[["last_lm_publication"]]
att1 <- attention_for(cut1)
sens <- bind_rows(
  compare_sets(og_outcome(cut1), att1, U, primary, "O1, cut-off 2015-08-06"),
  compare_sets(og_outcome(cut0, "female"), att0, U, primary, "O1, female-only records"),
  compare_sets(og_outcome(cut0, "male"), att0, U, primary, "O1, male-only records"),
  compare_sets(gwas_outcome(cut1), att1, U, primary, "O2, cut-off 2015-08-06"),
  compare_sets(gwas_outcome(cut0, TV$suggestive), att0, U, primary, "O2, P <= 1e-5"),
  compare_sets(gwas_outcome(cut0, in_gene_only = TRUE), att0, U, primary, "O2, variants inside the gene only"),
  compare_sets(gwas_outcome(cut0, traitset = "own"), att0, U, primary, "O2, own longevity traits only"),
  compare_sets(gwas_outcome(cut0, traitset = "parental"), att0, U, primary, "O2, parental longevity traits only"),
  compare_sets(gwas_outcome(cut0, traitset = "with_healthspan"), att0, U, primary, "O2, healthspan traits added"),
  compare_sets(o2, att0, U, intersect(primary, gene_out$gene[gene_out$stable_core]), "O2, stable core (13 genes)"),
  compare_sets(o1, att0, U, intersect(primary, gene_out$gene[gene_out$stable_core]), "O1, stable core (13 genes)")
)
write_csv(sens, file.path(out_dir, "sensitivity_analyses.csv"))
for (i in seq_len(nrow(sens))) with(sens[i, ], check_value(analysis, sprintf(
  "%d/%d vs %d/%d; Fisher P = %.2g; adjusted OR = %s, P = %s", confirmed_selected, n_selected, confirmed_other, n_other,
  fisher_p, signif(adjusted_or, 3), signif(adjusted_p, 2))))
message("\nResults written to ", out_dir)
