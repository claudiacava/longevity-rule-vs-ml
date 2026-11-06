# =============================================================================
# 00_setup.R: install required packages and record the software environment
# =============================================================================

cran_pkgs <- c("readr", "dplyr", "tidyr", "stringr", "purrr", "tibble",
               "ggplot2", "scales", "patchwork", "broom", "jsonlite",
               "openxlsx", "BiocManager")
bioc_pkgs <- c("AnnotationDbi", "org.Hs.eg.db", "GO.db", "clusterProfiler",
               "enrichplot", "multiMiR")

missing_cran <- setdiff(cran_pkgs, rownames(installed.packages()))
if (length(missing_cran)) {
  install.packages(missing_cran, repos = "https://cloud.r-project.org")
}
missing_bioc <- setdiff(bioc_pkgs, rownames(installed.packages()))
if (length(missing_bioc)) {
  BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)
}

# Locate the project root even if the working directory was not set:
# when this file is run with source("<path>/R/<script>.R"), move to <path>.
if (!file.exists(file.path("R", "config.R"))) {
  this_file <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NULL)
  if (!is.null(this_file)) setwd(dirname(dirname(this_file)))
}
source(file.path("R", "config.R"))

log_step("Software environment")
writeLines(capture.output(sessionInfo()), res_path("sessionInfo.txt"))

versions <- tibble(
  package = c(cran_pkgs, bioc_pkgs),
  version = vapply(c(cran_pkgs, bioc_pkgs),
                   function(p) as.character(utils::packageVersion(p)), character(1))
)
write_csv(versions, res_path("package_versions.csv"))
check_value("R version", R.version.string)
check_value("org.Hs.eg.db", as.character(packageVersion("org.Hs.eg.db")))
check_value("multiMiR", as.character(packageVersion("multiMiR")))
