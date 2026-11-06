# =============================================================================
# 30b_fix_doi_map.R: correct three Open Genes DOIs that could not be mapped to
# PubMed (found by the independent recomputation of 5 November 2026)
#   10.7554/eLife.39856.001  -> 10.7554/eLife.39856  (".001" suffix: Timmers et al. 2019)
#   0.1111/acel.12471        -> 10.1111/acel.12471   (leading "1" missing in Open Genes)
#   10.1016%2fj.neurobiolaging.2006.09.010 -> URL-encoded "/" decoded
# The corrected DOIs are queried in Europe PMC; the PMID is written into
# data/derived/doi_to_pmid_map_temporal.csv under the ORIGINAL (uncorrected) DOI,
# so that scripts 31-32 find it, and the publication date is added to
# data/derived/pmid_publication_dates.csv.
# =============================================================================
if (!file.exists(file.path("R", "config.R"))) {
  this_file <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NULL)
  if (!is.null(this_file)) setwd(dirname(dirname(this_file)))
}
source(file.path("R", "config.R"))
suppressPackageStartupMessages(library(jsonlite))
log_step("Correcting unmapped Open Genes DOIs")
epmc <- "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
fix <- tibble(doi = c("10.7554/elife.39856.001", "0.1111/acel.12471", "10.1016%2fj.neurobiolaging.2006.09.010"),
              doi_query = c("10.7554/eLife.39856", "10.1111/acel.12471", "10.1016/j.neurobiolaging.2006.09.010"))
map_file <- file.path(PATHS$derived, "doi_to_pmid_map_temporal.csv")
dates_file <- file.path(PATHS$derived, "pmid_publication_dates.csv")
map <- read_chr_csv(map_file); dates <- read_chr_csv(dates_file)
for (i in seq_len(nrow(fix))) {
  url <- paste0(epmc, "?query=", utils::URLencode(paste0('DOI:"', fix$doi_query[i], '"'), reserved = TRUE), "&format=json&resultType=lite&pageSize=5")
  res <- tryCatch(jsonlite::fromJSON(url)$resultList$result, error = function(e) NULL)
  if (is.null(res) || !NROW(res) || !"pmid" %in% names(res)) { message("  NOT FOUND: ", fix$doi_query[i]); next }
  r1 <- res[!is.na(res$pmid), ][1, ]
  pm <- as.character(r1$pmid)
  map <- filter(map, doi != fix$doi[i]) |> bind_rows(tibble(doi = fix$doi[i], pmid = pm, status = "mapped_after_doi_correction"))
  if (!pm %in% dates$pmid[dates$status == "found"]) {
    dates <- filter(dates, pmid != pm) |> bind_rows(tibble(pmid = pm, first_publication_date = as.character(r1$firstPublicationDate),
      pub_year = as.character(r1$pubYear), journal = as.character(r1$journalTitle), title = as.character(r1$title), status = "found"))
  }
  message(sprintf("  %-40s -> PMID %s (%s) %s", fix$doi_query[i], pm, r1$firstPublicationDate, substr(r1$title, 1, 60)))
  Sys.sleep(0.3)
}
write_csv(map, map_file, na = ""); write_csv(dates, dates_file, na = "")
message("Updated: ", map_file, " and ", dates_file)
