# =============================================================================
# 16_downstream_sensitivity.R: do the GO and miRNA conclusions depend on the
# selection thresholds? (Methods 2.4; Results 3.3; Supplementary 3)
#
# The GO enrichment (script 12) and the miRNA ranking + coverage-matched null
# (script 15) are repeated for every alternative candidate set produced by the
# selection-sensitivity analyses:
#   - the distinct gene sets of the 100-scenario threshold grid (script 11);
#   - the four alternative criteria definitions (script 10);
#   - the 13 genes selected under all criteria definitions ("stable core").
# Identical gene sets are analysed once. Sets with fewer than 3 genes are skipped.
# GO uses the same settings as the primary analysis (propagated annotations,
# term size 10-100, BH). The background is the source-reviewed 213-gene universe,
# except for the composite-records analysis, whose own universe is used.
# The miRNA null uses the same interaction table, coverage strata, background
# and number of iterations as script 15.
# =============================================================================

# Locate the project root even if the working directory was not set:
# when this file is run with source("<path>/R/<script>.R"), move to <path>.
if (!file.exists(file.path("R", "config.R"))) {
  this_file <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NULL)
  if (!is.null(this_file)) setwd(dirname(dirname(this_file)))
}
source(file.path("R", "config.R"))
suppressPackageStartupMessages({
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(GO.db)
  library(clusterProfiler)
})
# Bioconductor packages mask dplyr verbs; restore them explicitly.
select <- dplyr::select; filter <- dplyr::filter; rename <- dplyr::rename

out_dir <- res_path("16_downstream_sensitivity")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
log_step("Threshold sensitivity of the GO and miRNA results")
set.seed(PARAMS$seed)

min_set_size <- 3
focus_mirna  <- "hsa-miR-34a-5p"
split_set    <- function(x) { s <- str_split(x, ";\\s*")[[1]]; sort(unique(s[s != ""])) }

# ---- Candidate sets ------------------------------------------------------------
primary  <- read_lines(res_path("10_longevitymap", "primary_genes.txt"))
evidence <- read_csv(res_path("10_longevitymap", "gene_evidence_genetic.csv"), show_col_types = FALSE)
universe_213 <- sort(unique(evidence$gene))

grid <- read_csv(res_path("11_sensitivity", "threshold_scenarios_genetic.csv"), show_col_types = FALSE) |>
  mutate(definition = sprintf("grid: sig>=%d, pop>=%d, prop>=%.2f", min_sig, min_pop, min_prop),
         genes = map(coalesce(selected_genes, ""), split_set), universe = "213")

crit_long <- read_csv(res_path("10_longevitymap", "criteria_sensitivity_long.csv"), show_col_types = FALSE)
universe_composite <- sort(unique(crit_long$gene[crit_long$analysis == "composite_records_included"]))
crit <- crit_long |>
  filter(analysis != "primary") |>
  group_by(analysis) |>
  summarise(genes = list(sort(gene[selected])), .groups = "drop") |>
  mutate(definition = paste0("criteria: ", analysis),
         universe = if_else(analysis == "composite_records_included", "composite", "213"))

stable_core <- read_csv(res_path("10_longevitymap", "criteria_sensitivity_genes.csv"), show_col_types = FALSE) |>
  filter(n_analyses_selected == max(n_analyses_selected)) |> pull(gene) |> sort()

all_defs <- bind_rows(
  tibble(definition = "primary (sig>=2, pop>=2, prop>=0.60)", genes = list(sort(primary)), universe = "213"),
  tibble(definition = "stable core (selected under all criteria definitions)", genes = list(stable_core), universe = "213"),
  select(crit, definition, genes, universe),
  select(grid, definition, genes, universe)
) |>
  mutate(n_genes = lengths(genes), key = paste(universe, map_chr(genes, paste, collapse = ";")))

# One analysis per distinct gene set (and background)
sets <- all_defs |>
  group_by(key) |>
  summarise(definitions = paste(definition, collapse = " | "), genes = genes[1],
            universe = universe[1], n_genes = n_genes[1], .groups = "drop") |>
  mutate(is_primary = map_lgl(genes, ~ setequal(.x, primary))) |>
  arrange(desc(is_primary), n_genes) |>
  mutate(set_id = sprintf("S%02d", row_number()))
write_csv(all_defs |> left_join(select(sets, key, set_id), by = "key") |>
            transmute(definition, set_id, n_genes, genes = map_chr(genes, paste, collapse = "; ")),
          file.path(out_dir, "definition_to_set.csv"))
check_value("Candidate-set definitions / distinct gene sets",
            paste(nrow(all_defs), n_distinct(sets$key), sep = " / "))
check_value("Distinct sets skipped (< 3 genes)", sum(sets$n_genes < min_set_size))
sets <- filter(sets, n_genes >= min_set_size)

# ---- GO (same settings as script 12) -------------------------------------------
official <- keys(org.Hs.eg.db, keytype = "SYMBOL")
update_symbol <- function(sym) {
  if (sym %in% official) return(sym)
  if (paste0("MT-", sym) %in% official) return(paste0("MT-", sym))
  hit <- tryCatch(suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys = sym, keytype = "ALIAS",
                                                         columns = "SYMBOL")), error = function(e) NULL)
  cand <- unique(stats::na.omit(hit$SYMBOL))
  if (length(cand) == 1) cand else sym
}
all_symbols <- sort(unique(c(universe_213, universe_composite, unlist(sets$genes))))
to_current <- setNames(vapply(all_symbols, update_symbol, character(1)), all_symbols)
cur <- function(x) sort(unique(unname(to_current[x])))
backgrounds <- list("213" = cur(universe_213), composite = cur(universe_composite))

go_one <- function(genes, universe_key) {
  bg <- backgrounds[[universe_key]]
  ego <- suppressMessages(enrichGO(gene = cur(genes), universe = bg, OrgDb = org.Hs.eg.db,
                                   keyType = "SYMBOL", ont = "BP", pAdjustMethod = "BH",
                                   minGSSize = PARAMS$go_min_size, maxGSSize = PARAMS$go_max_size,
                                   pvalueCutoff = 1, qvalueCutoff = 1, readable = FALSE))
  if (is.null(ego) || nrow(ego@result) == 0)
    return(tibble(go_terms_tested = 0L, go_terms_fdr05 = 0L, go_terms_nominal = 0L,
                  go_top_term = NA_character_, go_top_p = NA_real_, go_top_fdr = NA_real_,
                  go_cholesterol_efflux_p = NA_real_, go_cholesterol_efflux_fdr = NA_real_))
  r <- arrange(as_tibble(ego@result), pvalue)
  sig <- filter(r, p.adjust <= PARAMS$go_fdr)
  if (nrow(sig)) go_sig_terms[[length(go_sig_terms) + 1]] <<-
    transmute(sig, set_id = current_set, ID, Description, GeneRatio, BgRatio, pvalue, p.adjust, geneID)
  ce <- filter(r, Description == "cholesterol efflux")
  tibble(go_terms_tested = nrow(r),
         go_terms_fdr05 = sum(r$p.adjust <= PARAMS$go_fdr, na.rm = TRUE),
         go_terms_nominal = sum(r$pvalue < 0.05, na.rm = TRUE),
         go_top_term = r$Description[1], go_top_p = r$pvalue[1], go_top_fdr = r$p.adjust[1],
         go_cholesterol_efflux_p = if (nrow(ce)) ce$pvalue[1] else NA_real_,
         go_cholesterol_efflux_fdr = if (nrow(ce)) ce$p.adjust[1] else NA_real_)
}
message("  GO enrichment for ", nrow(sets), " gene sets ...")
go_sig_terms <- list(); current_set <- NA_character_
go_res <- pmap_dfr(list(sets$genes, sets$universe, sets$set_id), function(g, u, id) {
  current_set <<- id
  go_one(g, u)
})
go_sig_tbl <- bind_rows(go_sig_terms)
if (nrow(go_sig_tbl) == 0) go_sig_tbl <- tibble(set_id = character(), Description = character())
go_sig_tbl <- left_join(go_sig_tbl, select(sets, set_id, n_genes, definitions), by = "set_id")
write_csv(go_sig_tbl, file.path(out_dir, "downstream_GO_significant_terms.csv"))
for (id in unique(go_sig_tbl$set_id)) {
  t <- filter(go_sig_tbl, set_id == id)
  message("  ", id, " (", t$n_genes[1], " genes) significant GO terms: ",
          paste(sprintf("%s (FDR %.3g)", t$Description, t$p.adjust), collapse = "; "))
}

# ---- miRNA ranking and coverage-matched null (same settings as script 15) ------
edges <- read_csv(res_path("15_mirna", "interactions_deduplicated_all_queried_genes.csv.gz"),
                  col_types = cols(.default = col_character()))
edges <- mutate(edges, interaction_id = as.integer(interaction_id))
edge_pmids <- edges |>
  transmute(interaction_id, mirna, pmid = extract_pmids(pubmed_id_raw)) |>
  unnest_longer(pmid, keep_empty = TRUE) |>
  filter(!is.na(pmid))
coverage <- read_csv(res_path("15_mirna", "null_background_coverage.csv"), show_col_types = FALSE)

mirna_levels <- sort(unique(edges$mirna))
pairs <- distinct(edges, mirna, gene)
targets_by_gene <- split(match(pairs$mirna, mirna_levels), pairs$gene)
count_targets <- function(genes) {
  idx <- unlist(targets_by_gene[intersect(genes, names(targets_by_gene))], use.names = FALSE)
  tabulate(idx, nbins = length(mirna_levels))
}

rank_mirnas <- function(genes) {
  base <- filter(edges, gene %in% genes)
  if (nrow(base) == 0) return(tibble())
  pubs <- edge_pmids |>
    filter(interaction_id %in% base$interaction_id) |>
    group_by(mirna) |> summarise(n_publications = n_distinct(pmid), .groups = "drop")
  base |>
    group_by(mirna) |>
    summarise(n_targets = n_distinct(gene), n_databases = n_distinct(database), .groups = "drop") |>
    left_join(pubs, by = "mirna") |>
    mutate(n_publications = coalesce(n_publications, 0L),
           score = n_targets * log2(1 + n_databases) * log2(2 + n_publications)) |>
    arrange(desc(score), desc(n_targets), mirna) |>
    mutate(rank = row_number())
}

null_one <- function(genes) {
  target_cov <- filter(coverage, gene %in% genes)
  need <- count(target_cov, stratum, name = "n")
  pools <- lapply(need$stratum, function(s) coverage$gene[coverage$stratum == s & !coverage$gene %in% genes])
  null_mat <- vapply(seq_len(PARAMS$null_iterations), function(i) {
    sampled <- unlist(Map(function(pool, n) sample(pool, n, replace = length(pool) < n), pools, need$n))
    count_targets(sampled)
  }, numeric(length(mirna_levels)))
  # computed outside tibble(): inside it, a new column named "observed" would
  # mask this vector and break the comparison below
  obs_all <- count_targets(genes)
  keep <- obs_all > 0
  obs <- obs_all[keep]
  null_keep <- null_mat[keep, , drop = FALSE]
  tibble(mirna = mirna_levels[keep], observed = obs,
         null_mean = rowMeans(null_keep),
         p = (1 + rowSums(null_keep >= obs)) / (PARAMS$null_iterations + 1)) |>
    mutate(fdr = p.adjust(p, "BH")) |>
    arrange(p, desc(observed))
}

mirna_one <- function(genes) {
  rk <- rank_mirnas(genes)
  if (nrow(rk) == 0) return(tibble(mirna_genes_with_interactions = 0L))
  nl <- null_one(genes)
  f_rank <- filter(rk, mirna == focus_mirna)
  f_null <- filter(nl, mirna == focus_mirna)
  tibble(mirna_genes_with_interactions = n_distinct(edges$gene[edges$gene %in% genes]),
         mirna_top = rk$mirna[1], mirna_top_targets = rk$n_targets[1],
         mir34a_rank = if (nrow(f_rank)) f_rank$rank else NA_integer_,
         mir34a_targets = if (nrow(f_rank)) f_rank$n_targets else 0L,
         mir34a_null_mean = if (nrow(f_null)) f_null$null_mean else NA_real_,
         mir34a_null_fdr = if (nrow(f_null)) f_null$fdr else NA_real_,
         null_min_fdr = min(nl$fdr), null_top_mirna = nl$mirna[1],
         null_mirnas_fdr05 = sum(nl$fdr < 0.05))
}
message("  miRNA ranking and ", PARAMS$null_iterations, "-iteration null for ", nrow(sets), " gene sets ...")
mirna_res <- map_dfr(sets$genes, mirna_one)

# ---- Summary --------------------------------------------------------------------
summary_tbl <- bind_cols(
  transmute(sets, set_id, n_genes, background = if_else(universe == "213", "213 source-reviewed genes",
                                                        "composite-records universe"),
            n_shared_with_primary = map_int(genes, ~ length(intersect(.x, primary))),
            definitions, genes = map_chr(genes, paste, collapse = "; ")),
  go_res, mirna_res)
write_csv(summary_tbl, file.path(out_dir, "downstream_sensitivity_summary.csv"))

check_value("Gene sets analysed", nrow(summary_tbl))
check_value("Sets with >= 1 GO term at FDR <= 0.05",
            paste(sum(summary_tbl$go_terms_fdr05 > 0), "of", nrow(summary_tbl)))
check_value("Sets with cholesterol efflux as the top GO term",
            paste(sum(summary_tbl$go_top_term == "cholesterol efflux", na.rm = TRUE), "of", nrow(summary_tbl)))
check_value(paste("Sets with", focus_mirna, "ranked first"),
            paste(sum(summary_tbl$mirna_top == focus_mirna, na.rm = TRUE), "of", nrow(summary_tbl)))
check_value("Sets with >= 1 miRNA above the coverage-matched null (FDR < 0.05)",
            paste(sum(summary_tbl$null_mirnas_fdr05 > 0, na.rm = TRUE), "of", nrow(summary_tbl)))
check_value(paste(focus_mirna, "null FDR range across sets"),
            paste(signif(range(summary_tbl$mir34a_null_fdr, na.rm = TRUE), 2), collapse = " - "))
for (i in which(!grepl("^grid", summary_tbl$definitions) | grepl("primary", summary_tbl$definitions))) {
  with(summary_tbl[i, ], message(sprintf("  %s (%d genes; %s): GO FDR<=0.05 = %d, top %s (P = %.2g, FDR = %.2g); top miRNA %s; miR-34a rank %s, null FDR %.2f",
                                         set_id, n_genes, substr(definitions, 1, 60), go_terms_fdr05, go_top_term,
                                         go_top_p, go_top_fdr, mirna_top, mir34a_rank, mir34a_null_fdr)))
}
