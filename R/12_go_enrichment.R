# =============================================================================
# 12_go_enrichment.R: GO Biological Process over-representation
# (Methods 2.3; Results 3.4)
#
# Primary analysis : clusterProfiler::enrichGO (annotations propagated to parent
#                    terms), universe = source-reviewed LongevityMap genes,
#                    term size 10-100 genes within the universe, BH correction.
# Sensitivity      : (a) direct annotations only;
#                    (b) leave-region-out for every 1-Mb proximity region that
#                    contains a candidate and for two predefined gene groups.
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

in_dir <- res_path("10_longevitymap")
out_dir <- res_path("12_go")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
log_step("GO Biological Process enrichment")

evidence <- read_csv(file.path(in_dir, "gene_evidence_genetic.csv"), show_col_types = FALSE)

# ---- Update outdated gene symbols ---------------------------------------------
# LongevityMap (2017) uses some symbols that have since been renamed (e.g. CARS
# -> CARS1, H2AFX -> H2AX) and mitochondrial genes without the "MT-" prefix.
# A symbol that is not a current official symbol is replaced by (i) "MT-" + symbol
# if that is an official symbol, otherwise (ii) the official symbol of a unique
# alias match in org.Hs.eg.db. Ambiguous or unknown symbols are left unchanged.
official <- keys(org.Hs.eg.db, keytype = "SYMBOL")
update_symbol <- function(sym) {
  if (sym %in% official) return(tibble(original = sym, symbol = sym, rule = "official"))
  if (paste0("MT-", sym) %in% official)
    return(tibble(original = sym, symbol = paste0("MT-", sym), rule = "mitochondrial prefix"))
  hit <- tryCatch(suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys = sym, keytype = "ALIAS",
                                                         columns = "SYMBOL")),
                  error = function(e) NULL)
  cand <- unique(stats::na.omit(hit$SYMBOL))
  if (length(cand) == 1) return(tibble(original = sym, symbol = cand, rule = "unique alias"))
  tibble(original = sym, symbol = sym,
         rule = if (length(cand) > 1) "ambiguous alias (unchanged)" else "not found (unchanged)")
}
symbol_map <- map_dfr(sort(unique(evidence$gene)), update_symbol)
write_csv(symbol_map, file.path(out_dir, "symbol_updates.csv"))
print(filter(symbol_map, rule != "official"))
to_current <- setNames(symbol_map$symbol, symbol_map$original)

universe <- sort(unique(unname(to_current[evidence$gene])))
primary <- sort(unique(unname(to_current[evidence$gene[evidence$primary_selected]])))

# ---- Mapping check ------------------------------------------------------------
entrez <- suppressMessages(bitr(universe, fromType = "SYMBOL", toType = "ENTREZID",
                                OrgDb = org.Hs.eg.db))
write_csv(entrez, file.path(out_dir, "symbol_to_entrez.csv"))
check_value("Symbols updated to current official symbols", sum(symbol_map$rule %in% c("mitochondrial prefix", "unique alias")))
check_value("Universe genes mapped to Entrez", n_distinct(entrez$SYMBOL))
check_value("Candidate genes mapped to Entrez", sum(primary %in% entrez$SYMBOL), 21L)
unmapped <- setdiff(universe, entrez$SYMBOL)
if (length(unmapped)) message("  unmapped symbols: ", paste(unmapped, collapse = ", "))

# ---- Primary analysis: propagated annotations --------------------------------
run_enrichgo <- function(candidates, background, label) {
  ego <- enrichGO(gene = candidates, universe = background, OrgDb = org.Hs.eg.db,
                  keyType = "SYMBOL", ont = "BP", pAdjustMethod = "BH",
                  minGSSize = PARAMS$go_min_size, maxGSSize = PARAMS$go_max_size,
                  pvalueCutoff = 1, qvalueCutoff = 1, readable = FALSE)
  if (is.null(ego)) return(list(table = tibble(analysis = label), object = NULL))
  # @result holds every tested term (as.data.frame() would apply cutoffs)
  tab <- as_tibble(ego@result) |>
    mutate(analysis = label, n_candidates = length(intersect(candidates, background)),
           n_background = length(background), .before = 1) |>
    arrange(p.adjust, pvalue)
  list(table = tab, object = ego)
}

primary_go <- run_enrichgo(primary, universe, "all_primary_candidates")
write_csv(primary_go$table, file.path(out_dir, "GO_BP_propagated_primary.csv"))
n_sig <- sum(primary_go$table$p.adjust <= PARAMS$go_fdr, na.rm = TRUE)
check_value("GO terms tested (propagated)", nrow(primary_go$table))
check_value("GO terms with FDR <= 0.05 (propagated)", n_sig)
print(head(select(primary_go$table, Description, GeneRatio, BgRatio, pvalue, p.adjust, geneID), 10))

if (!is.null(primary_go$object) && n_sig > 1) {
  simp <- clusterProfiler::simplify(primary_go$object, cutoff = PARAMS$go_simplify_cutoff,
                                    by = "p.adjust", select_fun = min)
  write_csv(as_tibble(simp@result), file.path(out_dir, "GO_BP_propagated_primary_simplified.csv"))
}

# ---- Sensitivity (a): direct annotations only ---------------------------------
go_direct <- AnnotationDbi::select(org.Hs.eg.db, keys = universe, keytype = "SYMBOL",
                                   columns = c("GO", "ONTOLOGY")) |>
  as_tibble() |>
  filter(ONTOLOGY == "BP", !is.na(GO)) |>
  distinct(SYMBOL, GO)
go_names <- AnnotationDbi::select(GO.db, keys = unique(go_direct$GO), keytype = "GOID",
                                  columns = "TERM") |> as_tibble() |> distinct(GOID, .keep_all = TRUE)

direct_hypergeometric <- function(candidates, background, label) {
  bg <- unique(go_direct$SYMBOL[go_direct$SYMBOL %in% background])
  cand <- intersect(candidates, bg)
  go_direct |>
    filter(SYMBOL %in% bg) |>
    group_by(GO) |>
    summarise(m = n_distinct(SYMBOL), k = sum(SYMBOL %in% cand),
              genes = paste(sort(SYMBOL[SYMBOL %in% cand]), collapse = "/"), .groups = "drop") |>
    filter(m >= PARAMS$go_min_size, k >= 2) |>
    mutate(p_value = phyper(k - 1, m, length(bg) - m, length(cand), lower.tail = FALSE),
           p_adjust = p.adjust(p_value, "BH"), analysis = label,
           n_candidates = length(cand), n_background = length(bg)) |>
    left_join(go_names, by = c("GO" = "GOID")) |>
    arrange(p_adjust, p_value)
}
direct <- direct_hypergeometric(primary, universe, "all_primary_candidates")
write_csv(direct, file.path(out_dir, "GO_BP_direct_primary.csv"))
check_value("Background genes with a direct BP annotation", direct$n_background[1])
check_value("GO terms tested (direct annotations)", nrow(direct))
check_value("Top direct term", paste0(direct$TERM[1], " (P = ", signif(direct$p_value[1], 3), ")"))

# ---- Sensitivity (b): leave-region-out ----------------------------------------
# Gene coordinates from org.Hs.eg.db (CHRLOC = start, CHRLOCEND = end; negative
# values indicate the minus strand). Alternative-haplotype contigs are dropped.
coords <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys = universe, keytype = "SYMBOL",
                                                 columns = c("CHRLOC", "CHRLOCEND"))) |>
  as_tibble() |>
  filter(!is.na(CHRLOC), CHRLOCCHR %in% c(1:22, "X", "Y")) |>
  mutate(start = pmin(abs(CHRLOC), abs(CHRLOCEND)), end = pmax(abs(CHRLOC), abs(CHRLOCEND))) |>
  group_by(SYMBOL, chr = CHRLOCCHR) |>
  summarise(start = min(start), end = max(end), .groups = "drop") |>
  group_by(SYMBOL) |> slice_min(start, n = 1, with_ties = FALSE) |> ungroup() |>
  arrange(chr, start) |>
  group_by(chr) |>
  mutate(reach = cummax(end),
         new_region = row_number() == 1 | start - lag(reach) > PARAMS$region_gap_bp,
         region = paste0("chr", chr, "_", cumsum(new_region))) |>
  ungroup()
write_csv(coords, file.path(out_dir, "gene_coordinates_regions.csv"))

candidate_regions <- coords |>
  filter(SYMBOL %in% primary) |>
  group_by(region) |>
  summarise(candidates = paste(sort(SYMBOL), collapse = "; "), .groups = "drop")

removal_sets <- c(
  setNames(lapply(candidate_regions$region, function(r) coords$SYMBOL[coords$region == r]),
           paste0("region_", candidate_regions$region)),
  setNames(PARAMS$predefined_groups, paste0("group_", names(PARAMS$predefined_groups)))
)

lro <- imap(removal_sets, function(removed, nm) {
  res <- run_enrichgo(setdiff(primary, removed), setdiff(universe, removed), paste0("without_", nm))
  res$table |> mutate(removed_genes = paste(sort(removed), collapse = "; "))
})
lro_tbl <- bind_rows(lro)
write_csv(lro_tbl, file.path(out_dir, "GO_BP_propagated_leave_region_out.csv"))

lro_summary <- lro_tbl |>
  group_by(analysis, removed_genes) |>
  summarise(n_candidates = first(n_candidates), n_terms = n(),
            n_fdr05 = sum(p.adjust <= PARAMS$go_fdr, na.rm = TRUE),
            top_term = first(Description), top_p = first(pvalue), top_fdr = first(p.adjust),
            .groups = "drop")
write_csv(lro_summary, file.path(out_dir, "GO_BP_leave_region_out_summary.csv"))
print(lro_summary, n = Inf)

# direct-annotation version of the two predefined groups
direct_lro <- imap_dfr(PARAMS$predefined_groups, function(removed, nm)
  direct_hypergeometric(setdiff(primary, removed), setdiff(universe, removed), paste0("without_", nm)))
write_csv(direct_lro, file.path(out_dir, "GO_BP_direct_leave_region_out.csv"))
