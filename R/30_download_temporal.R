# =============================================================================
# 30_download_temporal.R: evidence published after LongevityMap Build 3
#
# Follow-up study (temporal validation). Retrieves, as dated snapshots:
#   1. the GWAS Catalog association file (all studies, with ontology
#      annotations) and the subset of longevity-related traits;
#   2. the current Open Genes datasets (the July 2026 files used in the
#      published article are NOT modified);
#   3. publication dates from Europe PMC for every PMID of LongevityMap,
#      Open Genes and the selected GWAS studies.
#
# Output
#   data/raw/temporal/<date>/            raw snapshots (+ download_manifest.csv)
#   data/derived/temporal_trait_selection.csv   GWAS traits: review the
#       'include' column by hand, then run the script again
#   data/derived/pmid_publication_dates.csv     cached Europe PMC dates
#   data/derived/doi_to_pmid_map_temporal.csv   cached DOI -> PMID mapping
#   results/30_temporal/download_summary.csv
#
# The first run creates a snapshot named after today's date; later runs reuse
# the most recent snapshot (reproducible). Set FORCE_DOWNLOAD <- TRUE for a
# new snapshot. The full GWAS Catalog file is large (several hundred MB) and
# should not be committed to GitHub.
# =============================================================================

# Locate the project root even if the working directory was not set:
# when this file is run with source("<path>/R/<script>.R"), move to <path>.
if (!file.exists(file.path("R", "config.R"))) {
  this_file <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NULL)
  if (!is.null(this_file)) setwd(dirname(dirname(this_file)))
}
source(file.path("R", "config.R"))
suppressPackageStartupMessages(library(jsonlite))
options(timeout = max(3600, getOption("timeout")))

FORCE_DOWNLOAD <- FALSE

TEMPORAL <- list(
  # LongevityMap Build 3 was released in June 2017
  lm_release_date = as.Date("2017-06-01"),
  # Full association file with ontology annotations (tried in this order)
  # (file names as listed on the GWAS Catalog download page, November 2026)
  gwas_urls = c(
    "https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/gwas-catalog-associations_ontology-annotated-full.zip",
    "https://www.ebi.ac.uk/gwas/api/search/downloads/alternative"
  ),
  # Ancestry of every study (sample sizes and ancestral groups), small file
  gwas_ancestry_url = "https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/gwas-catalog-ancestry.tsv",
  gwas_manual_page = "https://www.ebi.ac.uk/gwas/docs/file-downloads",
  # Traits shown for review (searched in MAPPED_TRAIT and DISEASE/TRAIT)
  trait_candidate_regex = "longevity|life ?span|centenarian|age at death|healthspan|extreme age|exceptional age|survival to",
  # Default decision for new traits (MAPPED_TRAIT only); edit the selection file
  trait_include_regex = "longevity|life ?span|centenarian",
  epmc_url = "https://www.ebi.ac.uk/europepmc/webservices/rest/search",
  # NCBI gene -> PubMed links, used to measure research attention per gene
  gene2pubmed_url = "https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2pubmed.gz",
  epmc_batch = 20   # 100 PMIDs per query were rejected by Europe PMC (2 Nov 2026)
)

# ---- Snapshot folder ------------------------------------------------------------
temporal_root <- file.path(ROOT, "data", "raw", "temporal")
dir.create(temporal_root, recursive = TRUE, showWarnings = FALSE)
snapshots <- list.dirs(temporal_root, recursive = FALSE, full.names = FALSE)
snapshots <- sort(snapshots[grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}$", snapshots)])
snapshot <- if (length(snapshots) && !FORCE_DOWNLOAD) tail(snapshots, 1) else as.character(Sys.Date())
snap_dir <- file.path(temporal_root, snapshot)
dir.create(file.path(snap_dir, "open_genes"), recursive = TRUE, showWarnings = FALSE)
out_dir <- res_path("30_temporal")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
log_step("Temporal validation: downloads (snapshot ", snapshot, ")")

manifest_file <- file.path(snap_dir, "download_manifest.csv")
manifest <- if (file.exists(manifest_file)) read_chr_csv(manifest_file) else
  tibble(source = character(), file_name = character(), source_url = character(),
         access_date = character(), file_size_bytes = character(), md5 = character(), note = character())
record_download <- function(source, path, url, note = "") {
  row <- tibble(source = source, file_name = basename(path), source_url = url,
                access_date = as.character(Sys.Date()), file_size_bytes = as.character(file.size(path)),
                md5 = unname(tools::md5sum(path)), note = note)
  manifest <<- bind_rows(filter(manifest, file_name != basename(path)), row)
  write_csv(manifest, manifest_file)
}

split_mapped_genes <- function(x) {
  # MAPPED_GENE: "A" (in gene), "A, B" (several genes), "A - B" (intergenic:
  # upstream - downstream); haplotypes are separated by ";" or " x "
  str_split(norm_text(x), "\\s+-\\s+|\\s*[,;]\\s*|\\s+x\\s+")
}

# ---- 1. GWAS Catalog ----------------------------------------------------------
log_step("GWAS Catalog")
gwas_full <- file.path(snap_dir, "gwas_catalog_associations_full.tsv")
gwas_cand_file <- file.path(snap_dir, "gwas_catalog_longevity_candidate_associations.tsv")

if (!file.exists(gwas_full) && !file.exists(gwas_cand_file)) {
  ok <- FALSE
  for (u in TEMPORAL$gwas_urls) {
    message("  trying ", u, " (large file, may take several minutes)")
    tmp <- tempfile(fileext = ".download")
    res <- tryCatch(utils::download.file(u, tmp, mode = "wb", quiet = TRUE), error = function(e) e)
    if (inherits(res, "error") || !file.exists(tmp) || file.size(tmp) < 1e6) {
      message("   failed: ", if (inherits(res, "error")) conditionMessage(res) else "empty or very small file")
      unlink(tmp); next
    }
    if (identical(readBin(tmp, "raw", 2), as.raw(c(0x50, 0x4b)))) {   # zip archive
      inside <- utils::unzip(tmp, list = TRUE)$Name
      tsv <- inside[grepl("\\.tsv$", inside)][1]
      utils::unzip(tmp, files = tsv, exdir = snap_dir, junkpaths = TRUE)
      file.rename(file.path(snap_dir, basename(tsv)), gwas_full)
    } else {
      file.copy(tmp, gwas_full, overwrite = TRUE)
    }
    unlink(tmp)
    header <- readLines(gwas_full, n = 1, warn = FALSE)
    if (!all(c("PUBMEDID", "MAPPED_TRAIT", "MAPPED_GENE", "DATE") %in% strsplit(header, "\t")[[1]])) {
      message("   unexpected file format; trying the next address")
      unlink(gwas_full); next
    }
    record_download("GWAS Catalog associations (ontology annotated)", gwas_full, u)
    ok <- TRUE
    break
  }
  if (!ok) stop("GWAS Catalog download failed. Download 'All associations v1.0.2 - with added ontology annotations' (full, zip) from ",
                TEMPORAL$gwas_manual_page, ", unzip it, save it as\n  ", gwas_full, "\nand run this script again.")
}
gwas_anc <- file.path(snap_dir, "gwas_catalog_ancestry.tsv")
if (!file.exists(gwas_anc)) {
  res <- tryCatch(utils::download.file(TEMPORAL$gwas_ancestry_url, gwas_anc, mode = "wb", quiet = TRUE), error = function(e) e)
  if (inherits(res, "error")) { message("   ancestry file not downloaded: ", conditionMessage(res)); unlink(gwas_anc) } else
    record_download("GWAS Catalog study ancestry", gwas_anc, TEMPORAL$gwas_ancestry_url)
}
if (file.exists(gwas_full) && !any(manifest$file_name == basename(gwas_full)))
  record_download("GWAS Catalog associations (ontology annotated)", gwas_full, "manual download", "downloaded by hand")

if (!file.exists(gwas_cand_file)) {
  message("  reading the full file and extracting longevity-related traits")
  full <- readr::read_tsv(gwas_full, col_types = cols(.default = col_character()),
                          quote = "", progress = FALSE, na = c("", "NA"))
  check_value("GWAS Catalog associations (all traits)", nrow(full))
  rx <- regex(TEMPORAL$trait_candidate_regex, ignore_case = TRUE)
  cand <- full |>
    filter(str_detect(coalesce(MAPPED_TRAIT, ""), rx) | str_detect(coalesce(`DISEASE/TRAIT`, ""), rx))
  write_tsv(cand, gwas_cand_file, na = "")
  record_download("GWAS Catalog longevity-related candidate associations", gwas_cand_file,
                  basename(gwas_full), paste("filter:", TEMPORAL$trait_candidate_regex))
  rm(full); invisible(gc())
}
cand <- readr::read_tsv(gwas_cand_file, col_types = cols(.default = col_character()),
                        quote = "", progress = FALSE, na = c("", "NA"))

# Trait selection: one row per (mapped trait, reported trait). New rows get a
# default decision; decisions already in the file are kept, so manual edits
# survive later runs.
sel_file <- file.path(PATHS$derived, "temporal_trait_selection.csv")
traits <- cand |>
  group_by(mapped_trait = MAPPED_TRAIT, mapped_trait_uri = MAPPED_TRAIT_URI, disease_trait = `DISEASE/TRAIT`) |>
  summarise(n_associations = n(), n_studies = n_distinct(`STUDY ACCESSION`), n_pmids = n_distinct(PUBMEDID),
            first_date = if (all(is.na(DATE))) NA_character_ else min(DATE, na.rm = TRUE),
            last_date = if (all(is.na(DATE))) NA_character_ else max(DATE, na.rm = TRUE), .groups = "drop") |>
  mutate(default_include = str_detect(coalesce(mapped_trait, ""), regex(TEMPORAL$trait_include_regex, ignore_case = TRUE)))
old_sel <- if (file.exists(sel_file)) read_chr_csv(sel_file) else NULL
key <- function(d) paste(d$mapped_trait, d$disease_trait, sep = "\u0001")
traits <- traits |>
  mutate(include = if (is.null(old_sel)) default_include else
           ifelse(key(traits) %in% key(old_sel), as.logical(old_sel$include[match(key(traits), key(old_sel))]), default_include),
         decision = if (is.null(old_sel)) "default" else
           ifelse(key(traits) %in% key(old_sel), coalesce(old_sel$decision[match(key(traits), key(old_sel))], "manual"), "default"),
         note = if (is.null(old_sel) || !"note" %in% names(old_sel)) "" else
           coalesce(old_sel$note[match(key(traits), key(old_sel))], "")) |>
  arrange(desc(include), desc(n_studies), mapped_trait)
write_csv(traits, sel_file, na = "")
check_value("Candidate trait rows (review temporal_trait_selection.csv)", nrow(traits))
check_value("Trait rows included", sum(traits$include))

gwas_sel <- cand |>
  semi_join(filter(traits, include), by = c("MAPPED_TRAIT" = "mapped_trait", "DISEASE/TRAIT" = "disease_trait"))

# ---- 1b. NCBI gene2pubmed (research attention) ---------------------------------------
log_step("NCBI gene2pubmed")
g2p_file <- file.path(snap_dir, "gene2pubmed_human.tsv")
if (!file.exists(g2p_file)) {
  tmp <- tempfile(fileext = ".gz")
  message("  downloading ", TEMPORAL$gene2pubmed_url, " (large file)")
  res <- tryCatch(utils::download.file(TEMPORAL$gene2pubmed_url, tmp, mode = "wb", quiet = TRUE), error = function(e) e)
  if (inherits(res, "error")) {
    message("   failed: ", conditionMessage(res), " - research attention will fall back to LongevityMap counts")
  } else {
    # ~40 million rows for all species: read in blocks and keep human rows only
    g2p <- readr::read_tsv_chunked(tmp, readr::DataFrameCallback$new(function(x, pos) x[x[[1]] == "9606", 2:3]),
                                   chunk_size = 1e6, col_types = cols(.default = col_character()), progress = FALSE)
    names(g2p) <- c("gene_id", "pmid")
    write_tsv(g2p, g2p_file)
    record_download("NCBI gene2pubmed (human rows only)", g2p_file, TEMPORAL$gene2pubmed_url)
    check_value("gene2pubmed human gene-publication links", nrow(g2p))
    rm(g2p)
  }
  unlink(tmp)
} else message("  using existing file: ", basename(g2p_file))

# ---- 2. Open Genes, current version ------------------------------------------------
log_step("Open Genes (current version)")
OG_API <- "https://open-genes.com/api/"
flatten_api_table <- function(x) {
  if (is.list(x) && !is.data.frame(x) && "items" %in% names(x)) x <- x$items
  df <- as.data.frame(x, stringsAsFactors = FALSE)
  for (nm in names(df)) {
    if (is.list(df[[nm]])) {
      df[[nm]] <- vapply(df[[nm]], function(v) {
        if (is.null(v) || length(v) == 0) NA_character_ else as.character(jsonlite::toJSON(v, auto_unbox = TRUE))
      }, character(1))
    }
  }
  df
}
for (i in seq_len(nrow(OG_FILES))) {
  dest <- file.path(snap_dir, "open_genes", OG_FILES$file[i])
  if (file.exists(dest)) { message("  using existing file: ", OG_FILES$file[i]); next }
  url <- paste0(OG_API, OG_FILES$endpoint[i], "?pageSize=0")
  message("  downloading: ", url)
  df <- tryCatch(flatten_api_table(jsonlite::fromJSON(url, flatten = TRUE)), error = function(e) e)
  if (inherits(df, "error")) { message("   failed: ", conditionMessage(df)); next }
  write_csv(df, dest, na = "")
  record_download(paste("Open Genes:", OG_FILES$evidence_type[i]), dest, url)
}
og_new_file <- file.path(snap_dir, "open_genes", "human_longevity_associations.csv")
og_new <- if (file.exists(og_new_file)) read_chr_csv(og_new_file) else NULL
og_old <- read_chr_csv(file.path(PATHS$raw_og_dir, "human_longevity_associations.csv"))

# DOI -> PMID for the current Open Genes human records (cache seeded with the
# mapping used in the article)
norm_doi <- function(x) norm_text(x) |> str_to_lower() |> str_remove("^https?://(dx\\.)?doi\\.org/")
doi_cache_file <- file.path(PATHS$derived, "doi_to_pmid_map_temporal.csv")
doi_cache <- if (file.exists(doi_cache_file)) read_chr_csv(doi_cache_file) else
  read_chr_csv(file.path(PATHS$derived, "doi_to_pmid_map.csv"))
doi_cache <- filter(doi_cache, status != "query_error")
if (!is.null(og_new)) {
  todo <- if ("doi" %in% names(og_new)) setdiff(unique(norm_doi(og_new$doi)), c("", doi_cache$doi)) else character()
  message("  DOIs to map to PMIDs: ", length(todo))
  for (i in seq_along(todo)) {
    url <- paste0(TEMPORAL$epmc_url, "?query=", utils::URLencode(paste0('DOI:"', todo[i], '"'), reserved = TRUE),
                  "&format=json&resultType=lite&pageSize=5")
    out <- tryCatch(jsonlite::fromJSON(url), error = function(e) NULL)
    res <- out$resultList$result
    pm <- if (!is.null(res) && NROW(res) && "pmid" %in% names(res)) unique(stats::na.omit(res$pmid)) else character()
    doi_cache <- bind_rows(doi_cache, tibble(doi = todo[i], pmid = if (length(pm)) pm[1] else NA_character_,
      status = if (is.null(out)) "query_error" else if (length(pm) == 1) "mapped" else if (length(pm) > 1) "mapped_first_of_several" else "not_found"))
    if (i %% 25 == 0 || i == length(todo)) { write_csv(doi_cache, doi_cache_file); message("  ", i, "/", length(todo)) }
    Sys.sleep(0.15)
  }
  write_csv(doi_cache, doi_cache_file)
}
og_pmids <- function(d) {
  if (is.null(d)) return(character())
  for (col in c("pmid", "doi")) if (!col %in% names(d)) d[[col]] <- NA_character_
  direct <- str_extract(norm_text(d$pmid), "[0-9]{6,9}")
  via_doi <- doi_cache$pmid[match(norm_doi(d$doi), doi_cache$doi)]
  coalesce(direct, via_doi)
}

# ---- 3. Publication dates (Europe PMC) ---------------------------------------------
log_step("Publication dates (Europe PMC)")
lm_raw <- read_chr_csv(PATHS$raw_lm)
pmid_sources <- bind_rows(
  tibble(pmid = str_extract(norm_text(lm_raw$PubMed), "[0-9]{6,9}"), source = "LongevityMap"),
  tibble(pmid = og_pmids(og_old), source = "Open Genes (Jul 2026)"),
  tibble(pmid = og_pmids(og_new), source = paste0("Open Genes (", snapshot, ")")),
  tibble(pmid = str_extract(norm_text(gwas_sel$PUBMEDID), "[0-9]{6,9}"), source = "GWAS Catalog (selected traits)")
) |> filter(!is.na(pmid)) |> distinct()

dates_file <- file.path(PATHS$derived, "pmid_publication_dates.csv")
dates <- if (file.exists(dates_file)) read_chr_csv(dates_file) else
  tibble(pmid = character(), first_publication_date = character(), pub_year = character(),
         journal = character(), title = character(), status = character())
dates <- filter(dates, status != "query_error")
todo <- setdiff(unique(pmid_sources$pmid), dates$pmid)
message("  PMIDs: ", n_distinct(pmid_sources$pmid), "; to query: ", length(todo))
batches <- split(todo, ceiling(seq_along(todo) / TEMPORAL$epmc_batch))
for (b in seq_along(batches)) {
  ids <- batches[[b]]
  q <- paste0("(", paste0("EXT_ID:", ids, collapse = " OR "), ") AND SRC:MED")
  url <- paste0(TEMPORAL$epmc_url, "?query=", utils::URLencode(q, reserved = TRUE),
                "&format=json&resultType=lite&pageSize=1000")
  out <- tryCatch(jsonlite::fromJSON(url), error = function(e) NULL)
  res <- out$resultList$result
  got <- if (!is.null(res) && NROW(res)) {
    as_tibble(res) |>
      transmute(pmid = as.character(pmid),
                first_publication_date = if ("firstPublicationDate" %in% names(res)) as.character(firstPublicationDate) else NA_character_,
                pub_year = if ("pubYear" %in% names(res)) as.character(pubYear) else NA_character_,
                journal = if ("journalTitle" %in% names(res)) as.character(journalTitle) else NA_character_,
                title = if ("title" %in% names(res)) as.character(title) else NA_character_,
                status = "found") |>
      distinct(pmid, .keep_all = TRUE)
  } else tibble(pmid = character())
  missing <- setdiff(ids, got$pmid)
  dates <- bind_rows(dates, got, tibble(pmid = missing, status = if (is.null(out)) "query_error" else "not_found"))
  write_csv(dates, dates_file, na = "")
  message("  batch ", b, "/", length(batches), ": ", nrow(got), " found")
  Sys.sleep(0.3)
}
check_value("PMIDs with a publication date", sum(dates$status == "found" & dates$pmid %in% pmid_sources$pmid))
check_value("PMIDs not found / query errors",
            paste(sum(dates$status == "not_found"), sum(dates$status == "query_error"), sep = " / "))

# ---- 4. Summary ---------------------------------------------------------------------
log_step("Summary")
pub_date <- function(p) {
  d <- dates[match(p, dates$pmid), ]
  as.Date(coalesce(d$first_publication_date, ifelse(is.na(d$pub_year), NA_character_, paste0(d$pub_year, "-01-01"))), optional = TRUE)
}
cutoff <- TEMPORAL$lm_release_date
universe <- read_csv(res_path("10_longevitymap", "gene_evidence_genetic.csv"), show_col_types = FALSE)$gene
primary <- read_lines(res_path("10_longevitymap", "primary_genes.txt"))

lm_dates <- pub_date(unique(str_extract(norm_text(lm_raw$PubMed), "[0-9]{6,9}")))
og_new_post <- if (is.null(og_new)) NA else sum(pub_date(og_pmids(og_new)) >= cutoff, na.rm = TRUE)
gwas_sel <- gwas_sel |> mutate(pub = coalesce(as.Date(DATE, optional = TRUE), pub_date(str_extract(PUBMEDID, "[0-9]{6,9}"))), post = pub >= cutoff)
gwas_genes_post <- unique(unlist(split_mapped_genes(gwas_sel$MAPPED_GENE[gwas_sel$post %in% TRUE])))
gwas_genes_post <- gwas_genes_post[gwas_genes_post != ""]

summary_tbl <- tribble(
  ~item, ~value,
  "snapshot", snapshot,
  "LongevityMap cut-off date", as.character(cutoff),
  "LongevityMap: latest publication date of its PMIDs", as.character(max(lm_dates, na.rm = TRUE)),
  "LongevityMap: PMIDs published after the cut-off", as.character(sum(lm_dates >= cutoff, na.rm = TRUE)),
  "Open Genes (Jul 2026): human records", as.character(nrow(og_old)),
  "Open Genes (current): human records", if (is.null(og_new)) "not downloaded" else as.character(nrow(og_new)),
  "Open Genes (current): human records published after the cut-off", as.character(og_new_post),
  "GWAS Catalog: selected trait rows", as.character(sum(traits$include)),
  "GWAS Catalog: associations (selected traits)", as.character(nrow(gwas_sel)),
  "GWAS Catalog: studies (selected traits)", as.character(n_distinct(gwas_sel$`STUDY ACCESSION`)),
  "GWAS Catalog: associations published after the cut-off", as.character(sum(gwas_sel$post, na.rm = TRUE)),
  "GWAS Catalog: genome-wide significant (P <= 5e-8) after the cut-off",
    as.character(sum(gwas_sel$post & suppressWarnings(as.numeric(gwas_sel$`P-VALUE`)) <= 5e-8, na.rm = TRUE)),
  "GWAS Catalog: mapped genes after the cut-off", as.character(length(gwas_genes_post)),
  "  of which in the 213-gene LongevityMap universe", as.character(sum(gwas_genes_post %in% universe)),
  "  of which among the 21 primary genes", paste0(sum(primary %in% gwas_genes_post), " (",
                                                 paste(sort(primary[primary %in% gwas_genes_post]), collapse = ", "), ")")
)
write_csv(summary_tbl, file.path(out_dir, "download_summary.csv"))
for (i in seq_len(nrow(summary_tbl))) check_value(summary_tbl$item[i], summary_tbl$value[i])
message("\nNext: open data/derived/temporal_trait_selection.csv, check the 'include' column ",
        "(TRUE/FALSE), set 'decision' to 'manual' for rows you change, and run this script again.")
