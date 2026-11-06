# =============================================================================
# 15_mirna_analysis.R: miRNA rankings, assay-label sensitivity, matched null
# (Methods 2.9; Results 3.8; Supplementary 6-7; Figure 4)
#
# Interaction table : all multiMiR "validated" records (miRTarBase, TarBase,
#   miRecords). In TarBase the support_type field "positive/negative" describes
#   the direction of regulation, not whether the experiment confirmed the
#   interaction (e.g. miR-138-5p -> TERT, luciferase, PMID 18201269 is labelled
#   "negative" in TarBase and "Functional MTI" in miRTarBase), so TarBase records
#   are retained. miRTarBase records labelled "Non-Functional MTI" (experiments
#   that did not support a functional interaction) are excluded.
# Deduplication     : miRNA, gene, database, experiment, PubMed field. Records
#   differing only by an empty vs non-empty support_type are collapsed.
# Score             : n_targets * log2(1 + n_databases) * log2(2 + n_publications)
# Null              : for each candidate set, random gene sets of the same size
#   are drawn from the background (all queried genes not in the set), matched on
#   coverage strata (zero records, then deciles of log1p(records)).
# =============================================================================

# Locate the project root even if the working directory was not set:
# when this file is run with source("<path>/R/<script>.R"), move to <path>.
if (!file.exists(file.path("R", "config.R"))) {
  this_file <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NULL)
  if (!is.null(this_file)) setwd(dirname(dirname(this_file)))
}
source(file.path("R", "config.R"))
out_dir <- res_path("15_mirna")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
log_step("miRNA analysis")
set.seed(PARAMS$seed)

primary <- read_lines(res_path("10_longevitymap", "primary_genes.txt"))
additional <- read_csv(res_path("13_open_genes", "additional_candidates.csv"), show_col_types = FALSE)$gene
integrated <- sort(union(primary, additional))
background_genes <- read_lines(file.path(PATHS$raw_mm_dir, "queried_genes.txt"))
check_value("Primary / additional / integrated genes",
            paste(length(primary), length(additional), length(integrated), sep = " / "))

# ---- Interaction table -------------------------------------------------------
raw <- read_chr_csv(file.path(PATHS$raw_mm_dir, "multimir_validated_interactions.csv.gz"))
edges <- raw |>
  transmute(mirna = norm_text(mature_mirna_id), gene = norm_gene(target_symbol),
            database = norm_text(database), experiment = norm_text(experiment),
            support_type = norm_text(support_type), pubmed_id_raw = norm_text(pubmed_id)) |>
  filter(mirna != "", !is.na(gene)) |>
  filter(str_detect(mirna, "^hsa-")) |>      # human miRNAs only
  arrange(desc(support_type != "")) |>
  distinct(mirna, gene, database, experiment, pubmed_id_raw, .keep_all = TRUE)
check_value("Deduplicated rows before excluding non-functional records", nrow(edges))
check_value("miRTarBase 'Non-Functional MTI' rows excluded",
            sum(str_detect(edges$support_type, regex("non-functional", ignore_case = TRUE))))
edges <- edges |>
  filter(!str_detect(support_type, regex("non-functional", ignore_case = TRUE))) |>
  mutate(interaction_id = row_number(),
         truncated_pubmed_cell = nchar(pubmed_id_raw) >= PARAMS$pmid_truncation_nchar,
         assay_label = str_detect(str_to_lower(paste(experiment, support_type)), PARAMS$assay_regex))
write_csv(edges, file.path(out_dir, "interactions_deduplicated_all_queried_genes.csv.gz"))

edge_pmids <- edges |>
  transmute(interaction_id, mirna, pmid = extract_pmids(pubmed_id_raw)) |>
  unnest_longer(pmid, keep_empty = TRUE)

check_value("Retained interaction rows (all queried genes)", nrow(edges))
check_value("Rows labelled TarBase 'negative' (retained)",
            sum(edges$database == "tarbase" & edges$support_type == "negative"))

# ---- Rankings ----------------------------------------------------------------
rank_mirnas <- function(edge_data, genes, label) {
  base <- filter(edge_data, gene %in% genes)
  pubs <- edge_pmids |>
    filter(interaction_id %in% base$interaction_id, !is.na(pmid)) |>
    group_by(mirna) |> summarise(n_publications = n_distinct(pmid), .groups = "drop")
  base |>
    group_by(mirna) |>
    summarise(n_targets = n_distinct(gene), n_databases = n_distinct(database),
              n_interactions = n(), n_truncated_pubmed_cells = sum(truncated_pubmed_cell),
              targets = paste(sort(unique(gene)), collapse = "; "), .groups = "drop") |>
    left_join(pubs, by = "mirna") |>
    mutate(n_publications = coalesce(n_publications, 0L),
           miRNA_score = n_targets * log2(1 + n_databases) * log2(2 + n_publications),
           gene_set = label) |>
    arrange(desc(miRNA_score), desc(n_targets), mirna) |>
    mutate(rank = row_number(), .before = 1)
}

sets <- list(primary = primary, integrated = integrated)
assay_edges <- filter(edges, assay_label)
rankings <- list(
  primary          = rank_mirnas(edges, primary, "primary"),
  integrated       = rank_mirnas(edges, integrated, "integrated"),
  primary_assay    = rank_mirnas(assay_edges, primary, "primary, assay-labelled"),
  integrated_assay = rank_mirnas(assay_edges, integrated, "integrated, assay-labelled")
)
iwalk(rankings, ~ write_csv(.x, file.path(out_dir, paste0("ranking_", .y, ".csv"))))
iwalk(sets, ~ write_csv(filter(edges, gene %in% .x), file.path(out_dir, paste0("interactions_", .y, ".csv"))))

for (nm in names(sets)) {
  e <- filter(edges, gene %in% sets[[nm]])
  check_value(paste0(nm, ": interaction rows / miRNAs / genes with interactions"),
              paste(nrow(e), n_distinct(e$mirna), n_distinct(e$gene), sep = " / "))
  check_value(paste0(nm, ": genes without validated interactions"),
              paste(setdiff(sets[[nm]], e$gene), collapse = ", "))
  check_value(paste0(nm, ": top miRNA"),
              with(rankings[[nm]][1, ], paste0(mirna, " (", n_targets, " targets, score ", round(miRNA_score, 2), ")")))
}
integ_edges <- filter(edges, gene %in% integrated)
check_value("Assay-labelled rows within the integrated set",
            sprintf("%d of %d (%.2f%%)", sum(integ_edges$assay_label), nrow(integ_edges),
                    100 * mean(integ_edges$assay_label)))
check_value("Top assay-restricted integrated miRNA",
            with(rankings$integrated_assay[1, ], paste0(mirna, " (", n_targets, " targets)")))

# ---- Coverage-matched null -----------------------------------------------------
coverage <- tibble(gene = background_genes) |>
  left_join(count(edges, gene, name = "n_records"), by = "gene") |>
  mutate(n_records = coalesce(n_records, 0L))
nonzero <- coverage$n_records > 0
coverage$stratum <- "zero"
coverage$stratum[nonzero] <- paste0("decile_", dplyr::ntile(log1p(coverage$n_records[nonzero]), PARAMS$null_strata))
write_csv(coverage, file.path(out_dir, "null_background_coverage.csv"))

mirna_levels <- sort(unique(edges$mirna))
pairs <- distinct(edges, mirna, gene)
targets_by_gene <- split(match(pairs$mirna, mirna_levels), pairs$gene)

count_targets <- function(genes) {
  idx <- unlist(targets_by_gene[intersect(genes, names(targets_by_gene))], use.names = FALSE)
  tabulate(idx, nbins = length(mirna_levels))
}

null_test <- function(target_genes, label) {
  target_cov <- filter(coverage, gene %in% target_genes)
  need <- count(target_cov, stratum, name = "n")
  pools <- lapply(need$stratum, function(s) coverage$gene[coverage$stratum == s & !coverage$gene %in% target_genes])
  short <- need$n > lengths(pools)
  if (any(short)) warning(label, ": background too small in strata ", paste(need$stratum[short], collapse = ", "),
                          "; sampling with replacement there.")
  null_mat <- vapply(seq_len(PARAMS$null_iterations), function(i) {
    sampled <- unlist(Map(function(pool, n) sample(pool, n, replace = length(pool) < n), pools, need$n))
    count_targets(sampled)
  }, numeric(length(mirna_levels)))
  observed <- count_targets(target_genes)
  keep <- observed > 0
  tibble(mirna = mirna_levels[keep], observed_targets = observed[keep],
         null_mean = rowMeans(null_mat[keep, , drop = FALSE]),
         null_sd = apply(null_mat[keep, , drop = FALSE], 1, sd),
         empirical_p = (1 + rowSums(null_mat[keep, , drop = FALSE] >= observed[keep])) /
           (PARAMS$null_iterations + 1)) |>
    mutate(empirical_fdr = p.adjust(empirical_p, "BH"), gene_set = label,
           n_iterations = PARAMS$null_iterations,
           background_size = length(setdiff(coverage$gene, target_genes))) |>
    arrange(empirical_p, desc(observed_targets))
}

null_results <- list(primary = null_test(primary, "primary"),
                     integrated = null_test(integrated, "integrated"))
iwalk(null_results, ~ write_csv(.x, file.path(out_dir, paste0("null_comparison_", .y, ".csv"))))
for (nm in names(null_results)) {
  r <- filter(null_results[[nm]], mirna == "hsa-miR-34a-5p")
  if (nrow(r)) check_value(paste0(nm, ": miR-34a-5p observed / null mean / FDR"),
                           sprintf("%d / %.1f / %.3g", r$observed_targets, r$null_mean, r$empirical_fdr))
}
