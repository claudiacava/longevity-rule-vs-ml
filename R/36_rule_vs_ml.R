# =============================================================================
# 36_rule_vs_ml.R: the published rule against the machine-learning models,
# judged on the same outcomes published after the cut-off (2017-06-01).
#
# Analysis plan: 3 November 2026.
# This is the first and only script that reads the future outcomes together
# with the model scores (scores fixed by script 35).
#
#   Competitors (each gives a score to the 213 evaluation genes; top 21 kept,
#   as many as the rule selects)
#     R            published rule (score: selected first, then replication score)
#     M1_*         same LongevityMap evidence (script 35); primary = M1_LR_2013
#     M2_*         positive-unlabelled models on biological features
#     M3_*         average rank of R with M1_LR_2013 / with M2_ANN
#     B_pop        PubMed publications before the cut-off (gene2pubmed)
#     B_GO         number of GO Biological Process terms (2017; Fabris et al. 2020)
#     B_nsig       number of significant LongevityMap units
#   Outcomes      O1 Open Genes, O2 GWAS Catalog +/- 500 kb (script 31),
#                 O1 among genes re-tested in Open Genes after the cut-off
#   Metrics       confirmed genes and precision in the top 21 (primary: M1_LR_2013
#                 vs R on O1), precision in the top 10/30/50, average precision,
#                 Firth odds ratio of the top 21 adjusted for log publications,
#                 attention-matched random null (5 strata, 1000 draws),
#                 difference vs R with bootstrap 95% CI (10,000) and exact test
#                 on discordant genes; Holm correction for secondary comparisons;
#                 Spearman correlation of each score with publications
#   Ties at the top-21 boundary are broken by one fixed random order (seed).
#
# Inputs : results/35_ml_models, results/30_temporal, results/34_ml_features
# Output : results/36_rule_vs_ml/*.csv, figures/36_precision_at_k.pdf/.png
# =============================================================================

if (!file.exists(file.path("R", "config.R"))) {
  this_file <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NULL)
  if (!is.null(this_file)) setwd(dirname(dirname(this_file)))
}
source(file.path("R", "config.R"))
if (!requireNamespace("logistf", quietly = TRUE)) install.packages("logistf")

CMP <- list(seed = 20261104, k = 21, k_extra = c(10, 30, 50), n_null = 1000, n_boot = 10000,
            n_strata = 5, min_events = 3, primary = "M1_LR_2013")
out_dir <- res_path("36_rule_vs_ml"); dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
log_step("Rule versus machine learning on outcomes after 2017-06-01")
set.seed(CMP$seed)

# ---- Data ----------------------------------------------------------------------------------
sc  <- read_csv(res_path("35_ml_models", "scores_evaluation_genes.csv"), show_col_types = FALSE, col_types = cols(entrez = "c"))
oc  <- read_csv(res_path("30_temporal", "gene_outcomes.csv"), show_col_types = FALSE)
fc  <- read_csv(res_path("34_ml_features", "evaluation_gene_feature_counts.csv"), show_col_types = FALSE, col_types = cols(entrez = "c"))
d <- sc |>
  left_join(select(oc, gene, og_confirmed, og_retested, gwas_confirmed, pubs_pre), by = "gene") |>
  left_join(select(fc, gene, go_terms), by = "gene")
check_value("Evaluation genes", nrow(d), 213L)
check_value("Rule-selected / O1 confirmed / O2 confirmed", paste(sum(d$rule_selected), sum(d$og_confirmed), sum(d$gwas_confirmed), sep = " / "), "21 / 27 / 31")
rk <- function(x) rank(coalesce(x, -Inf), ties.method = "average")
d <- d |> mutate(R = as.numeric(rule_selected) * 100 + replication_score,
                 B_pop = pubs_pre, B_GO = go_terms, B_nsig = n_sig,
                 M3_R_M1 = (rk(R) + rk(M1_LR_2013)) / 2,
                 M3_R_M2ANN = (rk(R) + rk(M2_ANN)) / 2)
competitors <- tribble(
  ~name,        ~group,                                 ~label,
  "R",          "rule",                                 "Published rule",
  "M1_LR_2013", "same evidence",                        "M1 ridge logistic, 2013 (primary)",
  "M1_RF_2013", "same evidence",                        "M1 balanced random forest, 2013",
  "M1_LR_2011", "same evidence",                        "M1 ridge logistic, 2011",
  "M1_RF_2011", "same evidence",                        "M1 balanced random forest, 2011",
  "M2_ANN",     "biological, annotation-dependent",     "M2-A GO + pathways + PPI",
  "M2_GO",      "biological, annotation-dependent",     "M2 GO 2017",
  "M2_PATH",    "biological, annotation-dependent",     "M2 pathways (MSigDB 6.0)",
  "M2_PPI",     "biological, annotation-dependent",     "M2 BioGRID 3.4.149",
  "M2_EXPR",    "biological, less annotation-dependent","M2-B GTEx v8",
  "M3_R_M1",    "combination",                          "Rule + M1 (average rank)",
  "M3_R_M2ANN", "combination",                          "Rule + M2-A (average rank)",
  "B_pop",      "baseline",                             "Publications before cut-off",
  "B_GO",       "baseline",                             "Number of GO terms",
  "B_nsig",     "baseline",                             "Significant LongevityMap units")
stopifnot(all(competitors$name %in% names(d)))
outcomes <- list(
  O1 = list(y = d$og_confirmed, mask = rep(TRUE, nrow(d)), label = "Open Genes"),
  O2 = list(y = d$gwas_confirmed, mask = rep(TRUE, nrow(d)), label = "GWAS Catalog +/- 500 kb"),
  O1_retested = list(y = d$og_confirmed, mask = d$og_retested, label = "Open Genes, re-tested genes only"))

# ---- Helpers --------------------------------------------------------------------------------
tiebreak <- runif(nrow(d))
ord <- function(s) order(-coalesce(s, -Inf), tiebreak)                  # genes without data last
top <- function(s, k) seq_along(s) %in% ord(s)[seq_len(k)]
avg_precision <- function(s, y, mask) {
  o <- ord(s); o <- o[mask[o]]; yy <- y[o]
  if (sum(yy) == 0) return(NA_real_)
  sum((cumsum(yy) / seq_along(yy))[yy]) / sum(yy)
}
ties_at_boundary <- function(s, k) {
  v <- coalesce(s, -Inf); kth <- v[ord(s)[k]]
  if (sum(v == kth) > 1) sum(v == kth) else 0L
}
strata <- cut(rank(d$pubs_pre, ties.method = "first"), CMP$n_strata, labels = FALSE)
null_p <- function(set, y, mask) {
  obs <- sum(set & mask & y)
  idx <- split(seq_along(set), strata); need <- tabulate(strata[set], CMP$n_strata)
  sims <- replicate(CMP$n_null, {
    s <- unlist(lapply(seq_len(CMP$n_strata), function(h) idx[[h]][sample.int(length(idx[[h]]), need[h])]))
    sum(mask[s] & y[s])
  })
  c(null_mean = mean(sims), null_p = (1 + sum(sims >= obs)) / (1 + CMP$n_null))
}
firth <- function(set, y, mask) {
  dd <- tibble(y = as.numeric(y), inset = as.numeric(set), lp = log1p(d$pubs_pre))[mask, ]
  a <- sum(dd$inset == 1 & dd$y == 1); c <- sum(dd$inset == 0 & dd$y == 1)
  if (min(a, c) < CMP$min_events) return(c(or = NA, lo = NA, hi = NA, p = NA))
  tryCatch({
    m <- logistf::logistf(y ~ inset + lp, data = dd)
    k <- which(names(coef(m)) == "inset")
    c(or = exp(unname(coef(m)[k])), lo = exp(unname(m$ci.lower[k])), hi = exp(unname(m$ci.upper[k])), p = unname(m$prob[k]))
  }, error = function(e) c(or = NA, lo = NA, hi = NA, p = NA))
}

# ---- Main table -----------------------------------------------------------------------------
sets <- lapply(setNames(competitors$name, competitors$name), function(n) top(d[[n]], CMP$k))
rows <- list()
for (o in names(outcomes)) {
  y <- outcomes[[o]]$y; mask <- outcomes[[o]]$mask
  for (n in competitors$name) {
    s <- d[[n]]; set <- sets[[n]]
    nl <- null_p(set, y, mask); fr <- firth(set, y, mask)
    extra <- sapply(CMP$k_extra, function(k) { t <- top(s, k); sum(t & mask & y) / sum(t & mask) })
    rows[[length(rows) + 1]] <- tibble(
      outcome = o, outcome_label = outcomes[[o]]$label, competitor = n,
      genes_scored = sum(!is.na(s)), top21_in_outcome_set = sum(set & mask),
      top21_confirmed = sum(set & mask & y), top21_precision = sum(set & mask & y) / sum(set & mask),
      base_rate = sum(mask & y) / sum(mask),
      precision_top10 = extra[1], precision_top30 = extra[2], precision_top50 = extra[3],
      average_precision = avg_precision(s, y, mask),
      firth_or = fr[["or"]], firth_ci_low = fr[["lo"]], firth_ci_high = fr[["hi"]], firth_p = fr[["p"]],
      null_mean_confirmed = nl[["null_mean"]], null_p = nl[["null_p"]],
      ties_at_top21_boundary = ties_at_boundary(s, CMP$k),
      spearman_with_publications = suppressWarnings(cor(s, d$pubs_pre, method = "spearman", use = "complete.obs")))
  }
}
main <- bind_rows(rows) |> left_join(competitors, by = c("competitor" = "name")) |> relocate(group, label, .after = competitor)
write_csv(main, file.path(out_dir, "comparison_main.csv"), na = "")

# ---- Paired comparison with the rule ----------------------------------------------------------
paired <- list()
for (o in names(outcomes)) {
  y <- outcomes[[o]]$y; mask <- outcomes[[o]]$mask
  for (n in setdiff(competitors$name, "R")) {
    inM <- sets[[n]] & mask; inR <- sets[["R"]] & mask
    dlt <- (as.numeric(inM) - as.numeric(inR)) * y
    boot <- replicate(CMP$n_boot, sum(dlt[sample.int(length(dlt), replace = TRUE)]))
    m_only <- sum(inM & !inR & y); r_only <- sum(inR & !inM & y)
    bt <- if (m_only + r_only > 0) binom.test(m_only, m_only + r_only)$p.value else NA_real_
    paired[[length(paired) + 1]] <- tibble(
      outcome = o, competitor = n, confirmed_competitor = sum(inM & y), confirmed_rule = sum(inR & y),
      difference = sum(dlt), boot_ci_low = unname(quantile(boot, 0.025)), boot_ci_high = unname(quantile(boot, 0.975)),
      genes_shared_with_rule = sum(sets[[n]] & sets[["R"]]),
      confirmed_only_competitor = m_only, confirmed_only_rule = r_only, discordant_exact_p = bt,
      primary = (o == "O1" & n == CMP$primary))
  }
}
paired <- bind_rows(paired) |> group_by(outcome) |>
  mutate(holm_p = { p <- discordant_exact_p; p[primary] <- NA; p.adjust(p, "holm") }) |> ungroup()
write_csv(paired, file.path(out_dir, "paired_vs_rule.csv"), na = "")

# ---- Top-21 genes per competitor --------------------------------------------------------------
tops <- bind_rows(lapply(competitors$name, function(n) {
  i <- ord(d[[n]])[seq_len(CMP$k)]
  tibble(competitor = n, rank = seq_along(i), gene = d$gene[i], score = d[[n]][i], rule_selected = d$rule_selected[i],
         O1_confirmed = d$og_confirmed[i], O2_confirmed = d$gwas_confirmed[i], pubs_pre = d$pubs_pre[i])
}))
write_csv(tops, file.path(out_dir, "top21_genes.csv"), na = "")

# ---- Precision at k (figure) ----------------------------------------------------------------
fig_comp <- c("R", "M1_LR_2013", "M2_ANN", "M2_EXPR", "M3_R_M1", "B_pop")
pk <- bind_rows(lapply(c("O1", "O2"), function(o) bind_rows(lapply(fig_comp, function(n) {
  oo <- ord(d[[n]]); y <- outcomes[[o]]$y[oo]
  tibble(outcome = outcomes[[o]]$label, competitor = competitors$label[competitors$name == n], k = seq_along(y), precision = cumsum(y) / seq_along(y))
}))))
write_csv(pk, file.path(out_dir, "precision_at_k.csv"), na = "")
if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)
  base <- tibble(outcome = c(outcomes$O1$label, outcomes$O2$label), rate = c(mean(d$og_confirmed), mean(d$gwas_confirmed)))
  g <- ggplot(filter(pk, k <= 100), aes(k, precision, colour = competitor)) +
    geom_hline(data = base, aes(yintercept = rate), linetype = "dashed", colour = "grey50") +
    geom_vline(xintercept = CMP$k, colour = "grey80") +
    geom_line(linewidth = 0.7) + facet_wrap(~ outcome) +
    labs(x = "Genes ranked (top k)", y = "Precision (confirmed after 2017-06-01)", colour = NULL) +
    theme_minimal(base_size = 11) + theme(legend.position = "bottom") + guides(colour = guide_legend(ncol = 2))
  dir.create(file.path(ROOT, "figures"), showWarnings = FALSE)
  ggsave(file.path(ROOT, "figures", "36_precision_at_k.pdf"), g, width = 9, height = 5)
  ggsave(file.path(ROOT, "figures", "36_precision_at_k.png"), g, width = 9, height = 5, dpi = 200)
}

# ---- Report ----------------------------------------------------------------------------------
log_step("Top 21 genes: confirmed (precision) and comparison with the rule")
rep <- main |> filter(outcome %in% c("O1", "O2")) |>
  select(outcome, competitor, top21_confirmed, top21_precision, average_precision, null_p, spearman_with_publications) |>
  left_join(select(paired, outcome, competitor, difference, boot_ci_low, boot_ci_high, discordant_exact_p), by = c("outcome", "competitor")) |>
  mutate(across(where(is.numeric), ~ round(.x, 3)))
print(as.data.frame(rep), row.names = FALSE)
pr <- filter(paired, primary)
message(sprintf("\nPrimary comparison (O1, %s vs rule): %d vs %d confirmed, difference %+d (95%% CI %+.0f to %+.0f), discordant exact P = %.3f",
                CMP$primary, pr$confirmed_competitor, pr$confirmed_rule, pr$difference, pr$boot_ci_low, pr$boot_ci_high, pr$discordant_exact_p))
message("Written: ", out_dir)
