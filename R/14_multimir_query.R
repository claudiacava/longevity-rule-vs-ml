# =============================================================================
# 14_multimir_query.R: validated miRNA-target interactions from multiMiR
#
# Queried genes = combined LongevityMap + Open Genes gene universe (2,483 genes
# in the 24 July 2026 data). This covers both candidate sets and provides the
# background from which coverage-matched null gene sets are drawn (15_mirna).
# Queries run in chunks and each chunk is cached, so an interrupted run resumes.
#
# NB: multiMiR is loaded only through `multiMiR::` calls because the package
# attaches AnnotationDbi::select(), which masks dplyr::select().
# =============================================================================

# Locate the project root even if the working directory was not set:
# when this file is run with source("<path>/R/<script>.R"), move to <path>.
if (!file.exists(file.path("R", "config.R"))) {
  this_file <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NULL)
  if (!is.null(this_file)) setwd(dirname(dirname(this_file)))
}
source(file.path("R", "config.R"))
log_step("multiMiR validated interactions")

out_file <- file.path(PATHS$raw_mm_dir, "multimir_validated_interactions.csv.gz")
chunk_dir <- file.path(PATHS$raw_mm_dir, "chunks")
dir.create(chunk_dir, recursive = TRUE, showWarnings = FALSE)

lm_genes <- read_csv(res_path("10_longevitymap", "gene_evidence_genetic.csv"), show_col_types = FALSE)$gene
og_genes <- read_csv(res_path("13_open_genes", "open_genes_records_reconciled.csv"),
                     col_types = cols(.default = col_character()))$gene
query_genes <- sort(unique(norm_gene(c(lm_genes, og_genes))))
write_lines(query_genes, file.path(PATHS$raw_mm_dir, "queried_genes.txt"))
check_value("Genes to query", length(query_genes), 2483L)

# The interaction table used in the article is included in data/raw/multimir.
# Delete it to query the current multiMiR database again.
if (file.exists(out_file)) {
  message("  Using the existing interaction table: ", out_file)
} else {
  chunks <- split(query_genes, ceiling(seq_along(query_genes) / PARAMS$multimir_chunk))
  for (i in seq_along(chunks)) {
    f <- file.path(chunk_dir, sprintf("chunk_%03d.rds", i))
    if (file.exists(f)) next
    message("  chunk ", i, "/", length(chunks))
    res <- tryCatch(
      multiMiR::get_multimir(org = "hsa", target = chunks[[i]], table = "validated",
                             summary = FALSE),
      error = function(e) { message("    failed: ", conditionMessage(e)); NULL })
    if (is.null(res)) next                       # retried on the next run
    dat <- if (methods::is(res, "mmquery_bioc")) res@data else as.data.frame(res)
    saveRDS(as_tibble(dat), f)
  }

  chunk_files <- list.files(chunk_dir, pattern = "^chunk_\\d+\\.rds$", full.names = TRUE)
  if (length(chunk_files) < length(chunks)) {
    stop("Some chunks failed (", length(chunks) - length(chunk_files), "). Re-run this script.")
  }
  edges <- map_dfr(chunk_files, ~ mutate(readRDS(.x), across(everything(), as.character)))
  write_csv(edges, out_file)
  check_value("Raw validated interaction rows", nrow(edges))
  check_value("Genes with >= 1 validated interaction", n_distinct(edges$target_symbol))

  # Database versions for the Methods section
  info <- tryCatch(multiMiR::multimir_dbInfoVersions(), error = function(e) NULL)
  if (!is.null(info)) write_csv(as_tibble(info), file.path(PATHS$raw_mm_dir, "multimir_db_versions.csv"))
  info2 <- tryCatch(multiMiR::multimir_dbInfo(), error = function(e) NULL)
  if (!is.null(info2)) write_csv(as_tibble(info2), file.path(PATHS$raw_mm_dir, "multimir_db_info.csv"))
  writeLines(c(paste("query_date:", Sys.Date()),
               paste("multiMiR package:", as.character(packageVersion("multiMiR")))),
             file.path(PATHS$raw_mm_dir, "query_info.txt"))
}
