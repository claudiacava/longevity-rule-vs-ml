# =============================================================================
# 33_download_ml_features.R: biological features for the machine-learning
# models, in the version closest to (and not after) the cut-off 2017-06-01.
#
# Analysis plan (3 November 2026),
# section on information-leakage control.
#
# For every source the script tries a list of candidate URLs in order of
# preference, keeps the first that answers, and records what it used:
#   GO annotations   archived GO release, latest day <= 2017-06-01 (goa_human.gaf + go-basic.obo)
#   Pathways         Reactome release 60; fallback MSigDB 6.0 (April 2017) Reactome / KEGG gene sets
#   GO fallback      MSigDB 6.0 c5 Biological Process gene sets
#   PPI              STRING 10.5 (human, detailed scores per channel) + aliases
#                    BioGRID 3.4.149 (human physical interactions are filtered in script 34)
#   Expression       GTEx v8 gene median TPM (v6p has no public gene-median file)
#   Labels           GenAge human and models: current build (older builds only on request
#                    from HAGR); evaluation genes are never used as training labels
# Files already present are not downloaded again. Nothing is processed here.
#
# Output : data/raw/ml/<date>/  (+ url_check.csv, download_manifest.csv)
# =============================================================================

if (!file.exists(file.path("R", "config.R"))) {
  this_file <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NULL)
  if (!is.null(this_file)) setwd(dirname(dirname(this_file)))
}
source(file.path("R", "config.R"))
options(timeout = max(3600, getOption("timeout")))

CUTOFF <- as.Date("2017-06-01")
snap_dir <- file.path(ROOT, "data", "raw", "ml", format(Sys.Date()))
dir.create(snap_dir, recursive = TRUE, showWarnings = FALSE)
log_step("Machine-learning features (snapshot ", basename(snap_dir), ")")

go_days <- format(rev(seq(as.Date("2017-03-01"), CUTOFF, by = "day")))
msig <- "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/6.0/"
SOURCES <- list(
  go_gaf = list(file = "goa_human.gaf.gz", version = paste("GO release", go_days),
                urls = paste0("https://release.geneontology.org/", go_days, "/annotations/goa_human.gaf.gz")),
  go_obo = list(file = "go-basic.obo", version = paste("GO release", go_days),
                urls = paste0("https://release.geneontology.org/", go_days, "/ontology/go-basic.obo")),
  reactome_genes = list(file = "NCBI2Reactome_All_Levels.txt", version = "Reactome 60",
                        urls = c("https://download.reactome.org/60/NCBI2Reactome_All_Levels.txt",
                                 "https://reactome.org/download/60/NCBI2Reactome_All_Levels.txt")),
  reactome_pathways = list(file = "ReactomePathways.txt", version = "Reactome 60",
                           urls = c("https://download.reactome.org/60/ReactomePathways.txt",
                                    "https://reactome.org/download/60/ReactomePathways.txt")),
  msigdb_reactome = list(file = "c2.cp.reactome.v6.0.entrez.gmt", version = "MSigDB 6.0",
                         urls = paste0(msig, "c2.cp.reactome.v6.0.entrez.gmt")),
  msigdb_kegg = list(file = "c2.cp.kegg.v6.0.entrez.gmt", version = "MSigDB 6.0",
                     urls = paste0(msig, "c2.cp.kegg.v6.0.entrez.gmt")),
  msigdb_go_bp = list(file = "c5.bp.v6.0.entrez.gmt", version = "MSigDB 6.0",
                      urls = paste0(msig, "c5.bp.v6.0.entrez.gmt")),
  string_links = list(file = "9606.protein.links.detailed.v10.5.txt.gz", version = "STRING 10.5",
                      urls = c("https://version-10-5.string-db.org/download/protein.links.detailed.v10.5/9606.protein.links.detailed.v10.5.txt.gz",
                               "https://stringdb-static.org/download/protein.links.detailed.v10.5/9606.protein.links.detailed.v10.5.txt.gz")),
  string_aliases = list(file = "9606.protein.aliases.v10.5.txt.gz", version = "STRING 10.5",
                        urls = c("https://version-10-5.string-db.org/download/protein.aliases.v10.5/9606.protein.aliases.v10.5.txt.gz",
                                 "https://stringdb-static.org/download/protein.aliases.v10.5/9606.protein.aliases.v10.5.txt.gz")),
  biogrid = list(file = "BIOGRID-ORGANISM-3.4.149.tab2.zip", version = "BioGRID 3.4.149",
                 urls = c("https://downloads.thebiogrid.org/Download/BioGRID/Release-Archive/BIOGRID-3.4.149/BIOGRID-ORGANISM-3.4.149.tab2.zip",
                          "https://downloads.thebiogrid.org/File/BioGRID/Release-Archive/BIOGRID-3.4.149/BIOGRID-ORGANISM-3.4.149.tab2.zip")),
  # GTEx v6p has no gene-median file at its public address (404, 3 Nov 2026):
  # v8 is used and declared (expression is measurement, not curated knowledge)
  gtex_v8 = list(file = "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz", version = "GTEx v8",
                 urls = c("https://storage.googleapis.com/adult-gtex/bulk-gex/v8/rna-seq/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz",
                          "https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz")),
  genage_human = list(file = "human_genes.zip", version = paste("GenAge, downloaded", Sys.Date()),
                      urls = c("https://genomics.senescence.info/genes/human_genes.zip", "https://hagr.ageing-map.org/genes/human_genes.zip")),
  genage_models = list(file = "models_genes.zip", version = paste("GenAge, downloaded", Sys.Date()),
                       urls = c("https://genomics.senescence.info/genes/models_genes.zip", "https://hagr.ageing-map.org/genes/models_genes.zip"))
)

status_of <- function(u) {
  s <- tryCatch(attr(curlGetHeaders(u, redirect = TRUE, verify = TRUE), "status"), error = function(e) NA_integer_)
  if (is.null(s)) NA_integer_ else as.integer(s)
}

checks <- list(); manifest <- list()
for (nm in names(SOURCES)) {
  s <- SOURCES[[nm]]
  dest <- file.path(snap_dir, s$file)
  versions <- rep_len(s$version, length(s$urls))
  if (file.exists(dest) && file.size(dest) > 0) {
    message(sprintf("  %-18s already present", nm))
    manifest[[nm]] <- tibble(source = nm, file = s$file, url = NA_character_, version = "already present",
                             bytes = file.size(dest), status = "present")
    next
  }
  got <- FALSE
  for (i in seq_along(s$urls)) {
    # header check only to pick the GO release among many candidate days; some
    # servers refuse header requests, so other sources are downloaded directly
    many <- length(s$urls) > 3
    st <- if (many) status_of(s$urls[i]) else NA_integer_
    if (many && (is.na(st) || st != 200L)) {
      checks[[length(checks) + 1]] <- tibble(source = nm, url = s$urls[i], http_status = st, error = NA_character_)
      next
    }
    err <- NA_character_
    ok <- tryCatch({ download.file(s$urls[i], dest, mode = "wb", quiet = TRUE, method = "libcurl"); TRUE },
                   error = function(e) { err <<- conditionMessage(e); FALSE },
                   warning = function(w) { err <<- conditionMessage(w); FALSE })
    checks[[length(checks) + 1]] <- tibble(source = nm, url = s$urls[i], http_status = st, error = err)
    # a server can answer with an HTML page instead of the file: check the format
    if (ok && file.exists(dest)) {
      head_bytes <- readBin(dest, "raw", 4)
      is_html <- grepl("^\\s*<", rawToChar(readBin(dest, "raw", 200)[readBin(dest, "raw", 200) != as.raw(0)]))
      bad <- (grepl("\\.gz$", dest) && !identical(head_bytes[1:2], as.raw(c(0x1f, 0x8b)))) ||
             (grepl("\\.zip$", dest) && !identical(head_bytes[1:2], as.raw(c(0x50, 0x4b)))) ||
             (!grepl("\\.(gz|zip)$", dest) && is_html)
      if (bad) { ok <- FALSE; err <- "server returned a web page, not the file"; checks[[length(checks)]]$error <- err }
    }
    if (!ok) message("    ", basename(s$urls[i]), ": ", err)
    if (TRUE) {
      if (ok && file.exists(dest) && file.size(dest) > 0) {
        manifest[[nm]] <- tibble(source = nm, file = s$file, url = s$urls[i], version = versions[i],
                                 bytes = file.size(dest), status = "downloaded")
        message(sprintf("  %-18s %s (%.1f MB)", nm, versions[i], file.size(dest) / 1e6))
        got <- TRUE; break
      }
      if (file.exists(dest)) unlink(dest)
    }
  }
  if (!got) {
    manifest[[nm]] <- tibble(source = nm, file = s$file, url = NA_character_, version = NA_character_,
                             bytes = NA_real_, status = "not found")
    message(sprintf("  %-18s NOT FOUND at any candidate URL", nm))
  }
}
man <- bind_rows(manifest) |> mutate(downloaded_on = format(Sys.Date()))
write_csv(man, file.path(snap_dir, "download_manifest.csv"), na = "")
write_csv(bind_rows(checks), file.path(snap_dir, "url_check.csv"), na = "")

log_step("Summary")
print(as.data.frame(select(man, source, version, status, bytes)), row.names = FALSE)
message("\nWritten: ", snap_dir)
