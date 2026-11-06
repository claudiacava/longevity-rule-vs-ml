# =============================================================================
# 01_download_open_genes.R: retrieve the six Open Genes datasets
#
# The manuscript uses the files retrieved on 24 July 2026 (see
# data/raw/open_genes/open_genes_download_manifest.csv). If those files are
# present they are used as they are, so results are reproducible. Set
# FORCE_DOWNLOAD <- TRUE to retrieve the current database version instead
# (results may then differ from the manuscript).
# =============================================================================

# Locate the project root even if the working directory was not set:
# when this file is run with source("<path>/R/<script>.R"), move to <path>.
if (!file.exists(file.path("R", "config.R"))) {
  this_file <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NULL)
  if (!is.null(this_file)) setwd(dirname(dirname(this_file)))
}
source(file.path("R", "config.R"))
suppressPackageStartupMessages(library(jsonlite))

FORCE_DOWNLOAD <- FALSE
OG_API <- "https://open-genes.com/api/"

flatten_api_table <- function(x) {
  # The API returns either a bare array or an object with an "items" array.
  if (is.list(x) && !is.data.frame(x) && "items" %in% names(x)) x <- x$items
  df <- as.data.frame(x, stringsAsFactors = FALSE)
  # Nested list columns are stored as JSON strings so the CSV stays flat.
  for (nm in names(df)) {
    if (is.list(df[[nm]])) {
      df[[nm]] <- vapply(df[[nm]], function(v) {
        if (is.null(v) || length(v) == 0) NA_character_ else
          as.character(jsonlite::toJSON(v, auto_unbox = TRUE))
      }, character(1))
    }
  }
  df
}

log_step("Open Genes datasets")
manifest_rows <- list()
for (i in seq_len(nrow(OG_FILES))) {
  f <- OG_FILES$file[i]
  dest <- file.path(PATHS$raw_og_dir, f)
  url <- paste0(OG_API, OG_FILES$endpoint[i], "?pageSize=0")
  if (file.exists(dest) && !FORCE_DOWNLOAD) {
    message("  using existing file: ", f)
    next
  }
  message("  downloading: ", url)
  x <- jsonlite::fromJSON(url, flatten = TRUE)
  df <- flatten_api_table(x)
  write_csv(df, dest, na = "")
  manifest_rows[[f]] <- tibble(dataset = OG_FILES$evidence_type[i], file_name = f,
                               source_url = url, access_date = as.character(Sys.Date()),
                               file_size_bytes = file.size(dest))
}

manifest_file <- file.path(PATHS$raw_og_dir, "open_genes_download_manifest.csv")
if (length(manifest_rows)) {
  old <- if (file.exists(manifest_file)) read_chr_csv(manifest_file) else tibble()
  new <- bind_rows(manifest_rows) |> mutate(across(everything(), as.character))
  manifest <- bind_rows(old |> filter(!file_name %in% new$file_name), new)
  write_csv(manifest, manifest_file)
}
stopifnot(all(file.exists(file.path(PATHS$raw_og_dir, OG_FILES$file))))
if (file.exists(manifest_file)) print(read_chr_csv(manifest_file))
