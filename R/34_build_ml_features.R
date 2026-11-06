# =============================================================================
# 34_build_ml_features.R: feature tables for the rule-versus-machine-learning
# comparison. Nothing is modelled here; the models (script 35) read these files.
#
# Analysis plan: 3 November 2026.
#
# Gene identifier everywhere: NCBI Entrez Gene ID (org.Hs.eg.db).
# Universe for the genome-wide models: human protein-coding genes.
#
#   A. Same-evidence features (model M1), from LongevityMap gene-study units
#      published before a reference date:
#        train_2013  reference 2013-01-01, label = independent confirmation in
#                    [2013-01-01, 2017-06-01) (script 32)
#        train_2011  reference 2011-01-01 (sensitivity)
#        apply       reference 2017-06-01 = all LongevityMap evidence; the
#                    outcome (after 2017-06-01) is NOT attached here
#   B. Annotation-dependent features (model M2-A)
#        GO Biological Process, release 2017-05-30 (goa_human.gaf), excluding
#        NOT annotations, propagated to ancestors (is_a, part_of)
#        Reactome and KEGG gene sets, MSigDB 6.0 (April 2017)
#        BioGRID 3.4.149 human physical interactions
#   C. Less annotation-dependent features (model M2-B)
#        GTEx v8 median TPM by tissue, log2(TPM + 1)
#   D. Labels for the positive-unlabelled models: GenAge Human (Build 20);
#      the 213 evaluation genes are never labelled (always held out)
#   E. Summary: coverage of each family and its correlation with research
#      attention (PubMed publications before the cut-off) among the 213 genes
#
# Inputs : data/raw/ml/<snapshot>/ (script 33); results of 10, 12, 31, 32
# Output : data/derived/ml/*.csv.gz ; results/34_ml_features/*.csv
# =============================================================================

if (!file.exists(file.path("R", "config.R"))) {
  this_file <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NULL)
  if (!is.null(this_file)) setwd(dirname(dirname(this_file)))
}
source(file.path("R", "config.R"))
suppressPackageStartupMessages(library(org.Hs.eg.db))

FB <- list(
  refs        = c(train_2013 = "2013-01-01", train_2011 = "2011-01-01", apply = "2017-06-01"),
  min_genes_term = 10,     # GO terms / gene sets kept if annotated to >= 10 universe genes
  max_genes_term = 1000    # ... and to <= 1000 (very general terms carry little information)
)
ml_snaps <- sort(list.dirs(file.path(ROOT, "data", "raw", "ml"), recursive = FALSE, full.names = TRUE))
raw_dir <- tail(ml_snaps[grepl("[0-9]{4}-[0-9]{2}-[0-9]{2}$", ml_snaps)], 1)
stopifnot(length(raw_dir) == 1)
der_dir <- file.path(PATHS$derived, "ml"); dir.create(der_dir, recursive = TRUE, showWarnings = FALSE)
out_dir <- res_path("34_ml_features"); dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
log_step("ML features (raw snapshot ", basename(raw_dir), ")")
raw <- function(f) file.path(raw_dir, f)
gz  <- function(x, f) write_csv(x, file.path(der_dir, f), na = "")

# ---- Universe and evaluation genes -----------------------------------------------------
pc <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys = AnnotationDbi::keys(org.Hs.eg.db, "ENTREZID"),
                                             keytype = "ENTREZID", columns = c("SYMBOL", "GENETYPE"))) |>
  filter(GENETYPE == "protein-coding") |> distinct(ENTREZID, .keep_all = TRUE)
upd <- read_chr_csv(res_path("12_go", "symbol_updates.csv"))
s2e <- read_chr_csv(res_path("12_go", "symbol_to_entrez.csv"))
official <- function(g) { o <- upd$symbol[match(g, upd$original)]; str_to_upper(coalesce(o, g)) }
rec <- read_chr_csv(res_path("10_longevitymap", "gene_level_records_all_modalities.csv")) |>
  filter(evidence_modality == "genetic association")
dates <- read_chr_csv(file.path(PATHS$derived, "pmid_publication_dates.csv"))
rec$date <- as.Date(coalesce(dates$first_publication_date[match(rec$pub_med, dates$pmid)],
                             paste0(dates$pub_year[match(rec$pub_med, dates$pmid)], "-01-01")), optional = TRUE)
eval_genes <- tibble(gene = sort(unique(rec$gene))) |>
  # LOC symbols carry their Entrez Gene ID (LOC100631378 -> 100631378)
  mutate(entrez = coalesce(s2e$ENTREZID[match(official(gene), str_to_upper(s2e$SYMBOL))],
                           str_match(gene, "^LOC([0-9]+)$")[, 2]))
check_value("Evaluation genes / with Entrez ID", paste(nrow(eval_genes), sum(!is.na(eval_genes$entrez)), sep = " / "), "213 / 213")
universe <- tibble(entrez = union(pc$ENTREZID, na.omit(eval_genes$entrez))) |>
  mutate(symbol = unname(suppressMessages(AnnotationDbi::mapIds(org.Hs.eg.db, entrez, "SYMBOL", "ENTREZID"))),
         evaluation_gene = entrez %in% eval_genes$entrez)
check_value("Universe (protein-coding + evaluation genes)", nrow(universe))
gz(universe, "universe.csv.gz")

# ---- A. Same-evidence features (M1) ----------------------------------------------------
lm_features <- function(ref) {
  ref <- as.Date(ref)
  r <- filter(rec, date < ref)
  units <- r |> group_by(gene, study_unit) |>
    summarise(sig = !any(association == "non-significant"), date = min(date), .groups = "drop")
  g <- units |> group_by(gene) |>
    summarise(n_studies = n(), n_sig = sum(sig), n_nonsig = sum(!sig),
              years_since_first = as.numeric(ref - min(date)) / 365.25,
              years_since_last = as.numeric(ref - max(date)) / 365.25,
              years_since_last_sig = if (any(sig)) as.numeric(ref - max(date[sig])) / 365.25 else NA_real_,
              n_studies_last5y = sum(date >= ref - 5 * 365.25), .groups = "drop")
  p <- r |> left_join(select(units, gene, study_unit, sig), by = c("gene", "study_unit")) |>
    group_by(gene) |> summarise(n_populations = n_distinct(population), n_populations_sig = n_distinct(population[sig]),
                                n_variants = n_distinct(variant_s), .groups = "drop")
  tibble(gene = eval_genes$gene, entrez = eval_genes$entrez, reference_date = format(ref)) |>
    left_join(g, by = "gene") |> left_join(p, by = "gene") |>
    mutate(across(c(n_studies, n_sig, n_nonsig, n_studies_last5y, n_populations, n_populations_sig, n_variants), ~ coalesce(.x, 0L)),
           any_evidence = n_studies > 0,
           proportion_sig = if_else(n_studies > 0, n_sig / n_studies, 0),
           replication_score = log2(1 + n_studies) * (2 * proportion_sig - 1),
           rule_selected = n_sig >= PARAMS$min_sig_studies & n_populations >= PARAMS$min_populations &
             proportion_sig >= PARAMS$min_prop_sig - 1e-9)
}
fz <- res_path("32_ml_feasibility", "feasibility_gene_level.csv")
fz_tab <- if (file.exists(fz)) read_csv(fz, show_col_types = FALSE) else NULL
if (is.null(fz_tab) || !"independent" %in% names(fz_tab)) {
  message("  script 32 output without the 'independent' label: running script 32 first")
  source(file.path("R", "32_ml_feasibility.R"), local = new.env())
  fz_tab <- read_csv(fz, show_col_types = FALSE)
}
m1 <- bind_rows(lapply(names(FB$refs), function(s) {
  f <- lm_features(FB$refs[[s]]) |> mutate(set = s, .before = 1)
  if (s != "apply") {
    lab <- fz_tab |> filter(pseudo_cutoff == as.Date(FB$refs[[s]])) |> select(gene, label_independent = independent)
    f <- left_join(f, lab, by = "gene")
  }
  f
}))
check_value("M1 apply set: rule selects", sum(m1$rule_selected[m1$set == "apply"]), 21L)
check_value("M1 train_2013: genes with evidence / independent labels among them",
            paste(sum(m1$any_evidence[m1$set == "train_2013"]),
                  sum(m1$label_independent[m1$set == "train_2013" & m1$any_evidence], na.rm = TRUE), sep = " / "), "201 / 17")
gz(m1, "m1_same_evidence_features.csv.gz")

# ---- B1. GO Biological Process, release 2017 --------------------------------------------
log_step("GO Biological Process (archived release)")
gaf <- readr::read_tsv(raw("goa_human.gaf.gz"), comment = "!", col_names = FALSE, quote = "",
                       col_types = cols(.default = col_character()), progress = FALSE) |>
  filter(X9 == "P", !str_detect(coalesce(X4, ""), "NOT")) |>
  transmute(uniprot = X2, symbol = X3, go = X5, evidence = X7)
map_uni <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys = unique(gaf$uniprot), keytype = "UNIPROT", columns = "ENTREZID"))
map_sym <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys = unique(gaf$symbol), keytype = "SYMBOL", columns = "ENTREZID"))
map_ali <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys = unique(gaf$symbol), keytype = "ALIAS", columns = "ENTREZID")) |>
  group_by(ALIAS) |> filter(n_distinct(ENTREZID) == 1) |> ungroup()
gaf <- gaf |>
  left_join(filter(map_uni, !is.na(ENTREZID)) |> distinct(UNIPROT, .keep_all = TRUE), by = c("uniprot" = "UNIPROT")) |>
  mutate(entrez = coalesce(ENTREZID, map_sym$ENTREZID[match(symbol, map_sym$SYMBOL)],
                           map_ali$ENTREZID[match(symbol, map_ali$ALIAS)])) |>
  select(-ENTREZID)
check_value("GO BP annotations mapped to Entrez (%)", round(100 * mean(!is.na(gaf$entrez)), 1))
direct <- gaf |> filter(!is.na(entrez), entrez %in% universe$entrez) |> distinct(entrez, go)

# ontology: is_a and part_of parents, biological_process only, obsolete terms dropped
obo <- readLines(raw("go-basic.obo"), warn = FALSE)
stanza <- cumsum(startsWith(obo, "["))
is_term <- stanza %in% stanza[obo == "[Term]"]
ol <- tibble(stanza = stanza[is_term], line = obo[is_term])
terms <- ol |> group_by(stanza) |>
  summarise(go = sub("^id: ", "", line[startsWith(line, "id: ")][1]),
            name = sub("^name: ", "", line[startsWith(line, "name: ")][1]),
            namespace = sub("^namespace: ", "", line[startsWith(line, "namespace: ")][1]),
            obsolete = any(line == "is_obsolete: true"), .groups = "drop")
rels <- ol |>
  filter(startsWith(line, "is_a: ") | startsWith(line, "relationship: part_of ")) |>
  transmute(stanza, parent = str_extract(line, "GO:[0-9]{7}")) |>
  left_join(select(terms, stanza, go), by = "stanza")
bp <- filter(terms, namespace == "biological_process", !obsolete)
edges <- rels |> filter(go %in% bp$go, parent %in% bp$go) |> distinct(go, parent)
check_value("GO BP terms / parent links (ontology 2017-05-30)", paste(nrow(bp), nrow(edges), sep = " / "))
prop <- filter(direct, go %in% bp$go)
repeat {
  nxt <- prop |> inner_join(edges, by = "go", relationship = "many-to-many") |> transmute(entrez, go = parent) |>
    anti_join(prop, by = c("entrez", "go")) |> distinct()
  if (nrow(nxt) == 0) break
  prop <- bind_rows(prop, nxt)
}
go_size <- count(prop, go, name = "n_genes")
keep <- go_size$go[go_size$n_genes >= FB$min_genes_term & go_size$n_genes <= FB$max_genes_term]
go_long <- filter(prop, go %in% keep)
check_value("GO BP: genes annotated / terms kept", paste(n_distinct(go_long$entrez), length(keep), sep = " / "))
gz(go_long, "go_bp_2017_propagated_long.csv.gz")
gz(transmute(bp, go, name) |> filter(go %in% keep), "go_bp_2017_terms.csv.gz")

# ---- B2. Pathway gene sets, MSigDB 6.0 --------------------------------------------------
read_gmt <- function(f) {
  l <- strsplit(readLines(raw(f), warn = FALSE), "\t")
  bind_rows(lapply(l, function(x) tibble(set = x[1], entrez = x[-(1:2)])))
}
pw <- bind_rows(read_gmt("c2.cp.reactome.v6.0.entrez.gmt") |> mutate(source = "Reactome"),
                read_gmt("c2.cp.kegg.v6.0.entrez.gmt") |> mutate(source = "KEGG")) |>
  filter(entrez %in% universe$entrez) |> distinct()
pw_size <- count(pw, set, name = "n_genes")
pw <- filter(pw, set %in% pw_size$set[pw_size$n_genes >= FB$min_genes_term & pw_size$n_genes <= FB$max_genes_term])
check_value("Pathways: genes / sets (Reactome + KEGG)", paste(n_distinct(pw$entrez), n_distinct(pw$set), sep = " / "))
gz(pw, "pathways_msigdb6_long.csv.gz")

# ---- B3. BioGRID 3.4.149 physical interactions --------------------------------------------
zf <- raw("BIOGRID-ORGANISM-3.4.149.tab2.zip")
inner <- grep("Homo_sapiens", unzip(zf, list = TRUE)$Name, value = TRUE)[1]
bg <- readr::read_tsv(unz(zf, inner), col_types = cols(.default = col_character()), quote = "", progress = FALSE)
names(bg) <- make.names(sub("^#", "", names(bg)))
ppi <- bg |>
  filter(Experimental.System.Type == "physical", Organism.Interactor.A == "9606", Organism.Interactor.B == "9606",
         Entrez.Gene.Interactor.A != Entrez.Gene.Interactor.B) |>
  transmute(a = pmin(Entrez.Gene.Interactor.A, Entrez.Gene.Interactor.B), b = pmax(Entrez.Gene.Interactor.A, Entrez.Gene.Interactor.B),
            pmid = Pubmed.ID, throughput = Throughput) |>
  filter(a %in% universe$entrez, b %in% universe$entrez)
edges_ppi <- ppi |> group_by(a, b) |>
  summarise(n_publications = n_distinct(pmid), low_throughput = any(str_detect(throughput, "Low")), .groups = "drop")
check_value("BioGRID: human physical interactions (unique pairs) / genes",
            paste(nrow(edges_ppi), n_distinct(c(edges_ppi$a, edges_ppi$b)), sep = " / "))
gz(edges_ppi, "biogrid_3.4.149_physical_edges.csv.gz")

# ---- C. GTEx v8 median TPM -------------------------------------------------------------
gt <- readr::read_tsv(raw("GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz"), skip = 2,
                      show_col_types = FALSE, progress = FALSE) |>
  mutate(ensembl = sub("\\.[0-9]+$", "", Name))
ens <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys = unique(gt$ensembl), keytype = "ENSEMBL", columns = "ENTREZID")) |>
  filter(!is.na(ENTREZID)) |> group_by(ENSEMBL) |> filter(n() == 1) |> ungroup()
gtex <- gt |> inner_join(ens, by = c("ensembl" = "ENSEMBL")) |> filter(ENTREZID %in% universe$entrez) |>
  select(-Name, -Description, -ensembl) |>
  group_by(ENTREZID) |> summarise(across(everything(), mean), .groups = "drop") |>
  mutate(across(-ENTREZID, ~ log2(.x + 1))) |> rename(entrez = ENTREZID)
check_value("GTEx v8: genes / tissues", paste(nrow(gtex), ncol(gtex) - 1, sep = " / "))
gz(gtex, "gtex_v8_log2tpm_wide.csv.gz")

# ---- D. Positive labels (GenAge Human), evaluation genes held out ------------------------
zh <- raw("human_genes.zip")
ga <- read_csv(unz(zh, grep("\\.csv$", unzip(zh, list = TRUE)$Name, value = TRUE)[1]),
               col_types = cols(.default = col_character()), locale = locale(encoding = "latin1"))
names(ga) <- make.names(names(ga))
labels <- universe |>
  mutate(genage_human = entrez %in% ga$entrez.gene.id,
         label = case_when(evaluation_gene ~ NA_integer_, genage_human ~ 1L, TRUE ~ 0L),
         label_role = case_when(evaluation_gene ~ "held out (evaluation gene)", genage_human ~ "positive", TRUE ~ "unlabelled"))
check_value("GenAge Human genes in universe / of which evaluation genes (held out)",
            paste(sum(labels$genage_human), sum(labels$genage_human & labels$evaluation_gene), sep = " / "))
check_value("Training positives / unlabelled", paste(sum(labels$label == 1, na.rm = TRUE), sum(labels$label == 0, na.rm = TRUE), sep = " / "))
gz(labels, "pu_labels_genage_human.csv.gz")

# ---- E. Coverage and research attention among the evaluation genes -------------------------
att <- read_csv(res_path("30_temporal", "gene_outcomes.csv"), show_col_types = FALSE) |> select(gene, pubs_pre)
per_gene <- eval_genes |> left_join(att, by = "gene") |>
  mutate(go_terms = vapply(entrez, function(e) sum(go_long$entrez == e), 0L),
         pathway_sets = vapply(entrez, function(e) sum(pw$entrez == e), 0L),
         ppi_degree = vapply(entrez, function(e) sum(edges_ppi$a == e | edges_ppi$b == e), 0L),
         gtex_present = entrez %in% gtex$entrez,
         gtex_max_log2tpm = vapply(entrez, function(e) { r <- gtex[gtex$entrez == e, -1]; if (nrow(r)) max(unlist(r)) else NA_real_ }, 0))
write_csv(per_gene, file.path(out_dir, "evaluation_gene_feature_counts.csv"), na = "")
sp <- function(x) suppressWarnings(cor(x, per_gene$pubs_pre, method = "spearman", use = "complete.obs"))
fam <- tibble(
  family = c("GO BP 2017", "Pathways MSigDB 6.0", "BioGRID 3.4.149", "GTEx v8"),
  model = c("M2-A", "M2-A", "M2-A", "M2-B"),
  universe_genes_covered = c(n_distinct(go_long$entrez), n_distinct(pw$entrez), n_distinct(c(edges_ppi$a, edges_ppi$b)), nrow(gtex)),
  n_features = c(length(keep), n_distinct(pw$set), NA, ncol(gtex) - 1),
  evaluation_genes_without_data = c(sum(per_gene$go_terms == 0), sum(per_gene$pathway_sets == 0), sum(per_gene$ppi_degree == 0), sum(!per_gene$gtex_present)),
  spearman_with_publications = c(sp(per_gene$go_terms), sp(per_gene$pathway_sets), sp(per_gene$ppi_degree), sp(per_gene$gtex_max_log2tpm)),
  measure_correlated = c("number of GO terms", "number of gene sets", "interaction degree", "maximum expression across tissues"))
write_csv(fam, file.path(out_dir, "feature_family_summary.csv"), na = "")

log_step("Feature families (evaluation genes: 213)")
print(as.data.frame(fam), row.names = FALSE)
message("\nWritten: ", der_dir, " and ", out_dir)
