# =============================================================================
# 02_map_doi_to_pmid.R: convert Open Genes DOIs to PubMed identifiers
#
# LongevityMap identifies publications by PMID, whereas almost all Open Genes
# human-longevity records report only a DOI (1,412 of 1,460 in the 24 July 2026
# download). Without this mapping, publication-level overlap between the two
# databases cannot be detected. DOIs are resolved through the Europe PMC REST
# API; results are cached so the query runs only once.
# =============================================================================

# Locate the project root even if the working directory was not set:
# when this file is run with source("<path>/R/<script>.R"), move to <path>.
if (!file.exists(file.path("R", "config.R"))) {
  this_file <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NULL)
  if (!is.null(this_file)) setwd(dirname(dirname(this_file)))
}
source(file.path("R", "config.R"))
suppressPackageStartupMessages(library(jsonlite))

cache_file <- file.path(PATHS$derived, "doi_to_pmid_map.csv")

human <- read_chr_csv(file.path(PATHS$raw_og_dir, "human_longevity_associations.csv"))
dois <- human$doi |> norm_text() |> str_to_lower() |> str_remove("^https?://(dx\\.)?doi\\.org/")
dois <- sort(unique(dois[dois != ""]))

cache <- if (file.exists(cache_file)) read_chr_csv(cache_file) else
  tibble(doi = character(), pmid = character(), status = character())
# DOIs not yet queried, plus those whose query failed (retried on each run)
cache <- filter(cache, status != "query_error")
todo <- setdiff(dois, cache$doi)

lookup_doi <- function(doi) {
  url <- paste0("https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=",
                utils::URLencode(paste0('DOI:"', doi, '"'), reserved = TRUE),
                "&format=json&resultType=lite&pageSize=5")
  out <- tryCatch(jsonlite::fromJSON(url), error = function(e) NULL)
  if (is.null(out)) return(tibble(doi = doi, pmid = NA_character_, status = "query_error"))
  res <- out$resultList$result
  if (is.null(res) || !NROW(res) || !"pmid" %in% names(res)) {
    return(tibble(doi = doi, pmid = NA_character_, status = "not_found"))
  }
  pm <- unique(stats::na.omit(res$pmid))
  tibble(doi = doi, pmid = if (length(pm)) pm[1] else NA_character_,
         status = if (length(pm) == 1) "mapped" else if (length(pm) > 1) "mapped_first_of_several" else "no_pmid")
}

log_step("DOI -> PMID mapping")
message("  distinct DOIs: ", length(dois), "; to query: ", length(todo))
for (i in seq_along(todo)) {
  cache <- bind_rows(cache, lookup_doi(todo[i]))
  if (i %% 25 == 0 || i == length(todo)) {
    write_csv(cache, cache_file)
    message("  ", i, "/", length(todo))
  }
  Sys.sleep(0.15)  # be polite to the API
}
write_csv(cache, cache_file)
print(count(cache, status))
