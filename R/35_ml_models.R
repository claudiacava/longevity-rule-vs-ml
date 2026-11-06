# =============================================================================
# 35_ml_models.R: train the machine-learning competitors of the published
# rule and score the 213 evaluation genes. No future outcome is read here:
# the comparison with the outcomes after 2017-06-01 is script 36.
#
# Analysis plan: 3 November 2026.
#
#   M1  same LongevityMap evidence as the rule (script 34, set train_2013)
#       label = independent confirmation in [2013-01-01, 2017-06-01)
#       M1-LR  ridge logistic regression (glmnet), lambda by repeated CV
#       M1-RF  balanced random forest (ranger, equal class sampling per tree)
#       then applied to all LongevityMap evidence (set apply)
#       sensitivity: pseudo-cut-off 2011 (set train_2011)
#   M2  positive-unlabelled learning on the protein-coding genome, two-step
#       (Paz-Ruza et al., Comput Biol Med 2024):
#       1. reliable negatives: an unlabelled gene is a reliable negative if its
#          nearest neighbour is unlabelled AND the fraction of unlabelled genes
#          among its k nearest neighbours is >= t; similarity = Jaccard on the
#          most discriminative binary features (Euclidean on scaled values for
#          expression); k, t from the grid of Paz-Ruza (Table A.1)
#       2. balanced random forest, positives vs reliable negatives
#       k and t chosen by 5-fold cross-validation (F1 of the positive class);
#       reliable negatives re-selected inside every training fold
#       families: M2-GO, M2-PATH, M2-PPI, M2-ANN (GO + pathways + PPI) = M2-A;
#                 M2-EXPR (GTEx v8) = M2-B
#       positives = GenAge Human; the 213 evaluation genes are never in training
#   Deviations from Paz-Ruza, declared: genome-wide universe instead of ageing
#   genes only; random forest only (no CatBoost); one 5-fold CV for k and t.
#
# Inputs : data/derived/ml/*.csv.gz (script 34)
# Output : results/35_ml_models/*.csv
# =============================================================================

if (!file.exists(file.path("R", "config.R"))) {
  this_file <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NULL)
  if (!is.null(this_file)) setwd(dirname(dirname(this_file)))
}
source(file.path("R", "config.R"))
for (p in c("glmnet", "ranger", "Matrix")) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
suppressPackageStartupMessages({ library(Matrix); library(glmnet); library(ranger) })

ML <- list(
  seed = 20261103,
  m1_features = c("n_studies", "n_sig", "n_nonsig", "n_populations", "n_populations_sig", "n_variants",
                  "n_studies_last5y", "proportion_sig", "replication_score", "years_since_first",
                  "years_since_last", "years_since_last_sig", "any_sig"),
  m1_cv_repeats = 20, m1_oof_repeats = 10,
  pu_grid = tibble(k = c(3, 3, 5, 5, 8, 8, 8), t = c(2/3, 1, 4/5, 1, 6/8, 7/8, 1)),
  pu_folds = 5,
  n_features = 300,      # most discriminative features kept for the classifier (binary families)
  n_features_knn = 100,  # ... and for the Jaccard similarity
  ppi_min_partner_degree = 100,
  trees = 300,
  chunk = 500,
  threads = max(1, parallel::detectCores() - 1)
)
der_dir <- file.path(PATHS$derived, "ml")
# Negative control (leakage check): with ML_PERMUTE <- TRUE set before sourcing, the
# training labels are shuffled among the training genes, so the models cannot learn
# anything real; their predictions should then be no better than chance.
PERMUTE <- exists("ML_PERMUTE") && isTRUE(ML_PERMUTE)
out_dir <- res_path(if (PERMUTE) "35_ml_models_permuted" else "35_ml_models"); dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
if (PERMUTE) ML$seed <- ML$seed + 777
rd <- function(f) read_csv(file.path(der_dir, f), col_types = cols(.default = col_character()), progress = FALSE)
set.seed(ML$seed)
log_step("Machine-learning models (", ML$threads, " threads)", if (PERMUTE) " - LABELS PERMUTED (negative control)" else "")

# ---- Helpers -------------------------------------------------------------------------------
avg_precision <- function(score, y) {
  ok <- !is.na(score); score <- score[ok]; y <- y[ok]
  if (sum(y) == 0) return(NA_real_)
  o <- order(-score); y <- y[o]
  sum((cumsum(y) / seq_along(y))[y == 1]) / sum(y)
}
f1_pos <- function(pred, y) {
  tp <- sum(pred & y == 1); fp <- sum(pred & y == 0); fn <- sum(!pred & y == 1)
  if (tp == 0) 0 else 2 * tp / (2 * tp + fp + fn)
}
strat_folds <- function(y, k) {
  f <- integer(length(y))
  for (cl in unique(y)) { i <- which(y == cl); f[i] <- sample(rep_len(seq_len(k), length(i))) }
  f
}
fit_brf <- function(x, y, seed) {
  y <- factor(y, levels = c(0, 1))
  f <- sum(y == "1") / length(y)
  ranger(x = x, y = y, probability = TRUE, num.trees = ML$trees, sample.fraction = c(f, f), replace = TRUE,
         min.node.size = 5, num.threads = ML$threads, seed = seed)
}
prob1 <- function(m, x) predict(m, data = x, num.threads = ML$threads)$predictions[, "1"]

# ---- M1: same evidence as the rule -------------------------------------------------------
log_step("M1: LongevityMap evidence")
m1 <- rd("m1_same_evidence_features.csv.gz") |>
  mutate(across(c(n_studies, n_sig, n_nonsig, n_populations, n_populations_sig, n_variants, n_studies_last5y,
                  proportion_sig, replication_score, years_since_first, years_since_last, years_since_last_sig), as.numeric),
         across(c(any_evidence, rule_selected), as.logical),
         label_independent = as.logical(label_independent))
m1_matrix <- function(d) {
  d |> mutate(any_sig = as.numeric(n_sig > 0),
              years_since_last_sig = coalesce(years_since_last_sig, 40),  # no significant study yet
              across(c(n_studies, n_sig, n_nonsig, n_populations, n_populations_sig, n_variants, n_studies_last5y), log1p)) |>
    select(all_of(ML$m1_features)) |> as.matrix()
}
choose_lambda <- function(x, y) {
  base <- glmnet(x, y, family = "binomial", alpha = 0)
  cvm <- sapply(seq_len(ML$m1_cv_repeats), function(r) {
    cv <- cv.glmnet(x, y, family = "binomial", alpha = 0, lambda = base$lambda, foldid = strat_folds(y, 5),
                    type.measure = "deviance")
    cv$cvm[match(base$lambda, cv$lambda)]
  })
  base$lambda[which.min(rowMeans(cvm, na.rm = TRUE))]
}
run_m1 <- function(train_set) {
  tr <- filter(m1, set == train_set, any_evidence)
  ap <- filter(m1, set == "apply")
  x <- m1_matrix(tr); y <- as.numeric(tr$label_independent)
  if (PERMUTE) y <- sample(y)
  mu <- colMeans(x); sdv <- apply(x, 2, sd); sdv[sdv == 0] <- 1
  sc <- function(z) sweep(sweep(z, 2, mu, "-"), 2, sdv, "/")
  lam <- choose_lambda(sc(x), y)
  lr <- glmnet(sc(x), y, family = "binomial", alpha = 0, lambda = lam)
  rf <- fit_brf(x, y, ML$seed)
  # out-of-fold average precision on the training window (internal check only)
  oof <- sapply(seq_len(ML$m1_oof_repeats), function(r) {
    fo <- strat_folds(y, 5); p_lr <- p_rf <- numeric(length(y))
    for (k in 1:5) {
      a <- fo != k; b <- fo == k
      m_k <- colMeans(x[a, ]); s_k <- apply(x[a, ], 2, sd); s_k[s_k == 0] <- 1
      sc_k <- function(z) sweep(sweep(z, 2, m_k, "-"), 2, s_k, "/")
      g <- glmnet(sc_k(x[a, ]), y[a], family = "binomial", alpha = 0, lambda = lam)
      p_lr[b] <- predict(g, sc_k(x[b, , drop = FALSE]), type = "response")[, 1]
      p_rf[b] <- prob1(fit_brf(x[a, ], y[a], ML$seed + r * 10 + k), x[b, , drop = FALSE])
    }
    c(lr = avg_precision(p_lr, y), rf = avg_precision(p_rf, y))
  })
  coefs <- tibble(training = train_set, feature = ML$m1_features, ridge_coefficient_standardised = as.numeric(coef(lr))[-1],
                  rf_importance = NA_real_)
  scores <- tibble(gene = ap$gene, entrez = ap$entrez,
                   lr = as.numeric(predict(lr, sc(m1_matrix(ap)), type = "response")),
                   rf = prob1(rf, m1_matrix(ap)))
  names(scores)[3:4] <- paste0("M1_", c("LR", "RF"), "_", sub("train_", "", train_set))
  perf <- tibble(model = paste0("M1_", c("LR", "RF"), "_", sub("train_", "", train_set)), training_genes = length(y),
                 training_positives = sum(y), lambda = c(lam, NA), k = NA, t = NA,
                 cv_avg_precision = rowMeans(oof), cv_f1 = NA_real_, reliable_negatives = NA_integer_)
  check_value(paste("M1", train_set, ": training genes / positives / OOF avg precision LR, RF"),
              paste(length(y), sum(y), paste(round(rowMeans(oof), 3), collapse = ", "), sep = " / "))
  list(scores = scores, perf = perf, coefs = coefs)
}
# trained models are cached: delete results/35_ml_models/models_cache.rds to retrain
cache_file <- file.path(out_dir, "models_cache.rds")
cache <- if (file.exists(cache_file)) readRDS(cache_file) else NULL
if (!is.null(cache)) message("  using trained models from ", cache_file)
m1_runs <- if (!is.null(cache)) cache$m1_runs else lapply(c("train_2013", "train_2011"), run_m1)

# ---- M2: feature matrices --------------------------------------------------------------------
log_step("M2: feature matrices")
uni <- rd("universe.csv.gz") |> mutate(evaluation_gene = as.logical(evaluation_gene))
lab <- rd("pu_labels_genage_human.csv.gz") |> mutate(label = as.integer(label))
uni$label <- lab$label[match(uni$entrez, lab$entrez)]
ix <- setNames(seq_len(nrow(uni)), uni$entrez)
sparse_from <- function(long, col) {
  long <- distinct(filter(long, entrez %in% uni$entrez), entrez, .data[[col]])
  feats <- sort(unique(long[[col]]))
  sparseMatrix(i = ix[long$entrez], j = match(long[[col]], feats), x = 1, dims = c(nrow(uni), length(feats)),
               dimnames = list(uni$entrez, feats))
}
X_go <- sparse_from(rd("go_bp_2017_propagated_long.csv.gz"), "go")
X_pw <- sparse_from(rd("pathways_msigdb6_long.csv.gz"), "set")
ed <- rd("biogrid_3.4.149_physical_edges.csv.gz")
deg <- table(c(ed$a, ed$b))
hubs <- names(deg)[deg >= ML$ppi_min_partner_degree]
ppi_long <- bind_rows(transmute(ed, entrez = a, partner = b), transmute(ed, entrez = b, partner = a)) |> filter(partner %in% hubs)
X_ppi <- sparse_from(ppi_long, "partner"); colnames(X_ppi) <- paste0("ppi_", colnames(X_ppi))
gt <- rd("gtex_v8_log2tpm_wide.csv.gz")
X_ex <- matrix(NA_real_, nrow(uni), ncol(gt) - 1, dimnames = list(uni$entrez, names(gt)[-1]))
X_ex[ix[gt$entrez], ] <- as.matrix(mutate(gt[, -1], across(everything(), as.numeric)))
families <- list(
  M2_GO   = list(X = X_go, binary = TRUE),
  M2_PATH = list(X = X_pw, binary = TRUE),
  M2_PPI  = list(X = X_ppi, binary = TRUE),
  M2_ANN  = list(X = cbind(X_go, X_pw, X_ppi), binary = TRUE),
  M2_EXPR = list(X = X_ex, binary = FALSE))
for (nm in names(families)) {
  X <- families[[nm]]$X
  has <- if (families[[nm]]$binary) Matrix::rowSums(X) > 0 else rowSums(is.na(X)) == 0
  families[[nm]]$has <- has
  check_value(paste(nm, ": features / genes with data / positives with data / evaluation genes with data"),
              paste(ncol(X), sum(has), sum(has & uni$label %in% 1), sum(has & uni$evaluation_gene), sep = " / "))
}

# ---- M2: two-step PU learning -------------------------------------------------------------------
disc_score <- function(X, y) {
  p1 <- Matrix::colMeans(X[y == 1, , drop = FALSE]); p0 <- Matrix::colMeans(X[y == 0, , drop = FALSE])
  pb <- Matrix::colMeans(X)
  abs(p1 - p0) / sqrt(pb * (1 - pb) + 1e-9)
}
# for each unlabelled row: is each of its kmax nearest neighbours (among all rows) unlabelled?
knn_unlabelled <- function(S, is_u, kmax, binary) {
  u <- which(is_u); out <- matrix(NA, length(u), kmax)
  if (binary) { S <- as(S, "CsparseMatrix"); sz <- Matrix::rowSums(S) } else { S <- as.matrix(S); sq <- rowSums(S^2) }
  for (ch in split(seq_along(u), ceiling(seq_along(u) / ML$chunk))) {
    r <- u[ch]
    if (binary) {
      inter <- as.matrix(tcrossprod(S[r, , drop = FALSE], S))
      un <- outer(sz[r], sz, "+") - inter
      sim <- ifelse(un > 0, inter / un, 0)
    } else {
      sim <- -(outer(sq[r], sq, "+") - 2 * tcrossprod(S[r, , drop = FALSE], S))
    }
    sim <- sim + runif(length(sim), 0, 1e-9)      # random tie-breaking
    sim[cbind(seq_along(r), r)] <- -Inf           # a gene is not its own neighbour
    out[ch, ] <- t(apply(sim, 1, function(s) {
      thr <- -sort(-s, partial = kmax)[kmax]
      cand <- which(s >= thr)
      is_u[cand[order(-s[cand])][seq_len(kmax)]]
    }))
  }
  out
}
reliable <- function(nn, k, t) nn[, 1] & rowMeans(nn[, seq_len(k), drop = FALSE]) >= t - 1e-9
prep <- function(X, y, binary) {
  if (!binary) return(list(cls = colnames(X), knn = colnames(X)))
  s <- disc_score(X, y)
  o <- order(-s)
  list(cls = colnames(X)[o[seq_len(min(ML$n_features, length(o)))]],
       knn = colnames(X)[o[seq_len(min(ML$n_features_knn, length(o)))]])
}
scale_with <- function(X, ref) { mu <- colMeans(ref); s <- apply(ref, 2, sd); s[s == 0] <- 1; sweep(sweep(X, 2, mu, "-"), 2, s, "/") }
# dense block of the selected columns, with plain column names (GO:..., pathway names)
mx <- function(X, rows, cols) { m <- as.matrix(X[rows, cols, drop = FALSE]); colnames(m) <- paste0("v", seq_along(cols)); m }
pu_fit <- function(X, y, binary, grid_row, nn, f, seed) {
  rn <- reliable(nn, grid_row$k, grid_row$t)
  keep <- y == 1; keep[which(y == 0)[rn]] <- TRUE
  fit_brf(mx(X, keep, f$cls), y[keep], seed)
}
run_pu <- function(nm) {
  fam <- families[[nm]]; X <- fam$X
  tr <- which(fam$has & !is.na(uni$label))           # training genes (evaluation genes excluded)
  ev <- which(fam$has & uni$evaluation_gene)
  y <- uni$label[tr]
  if (PERMUTE) y <- sample(y)
  kmax <- max(ML$pu_grid$k)
  folds <- strat_folds(y, ML$pu_folds)
  cvres <- list()
  for (fo in seq_len(ML$pu_folds)) {
    a <- tr[folds != fo]; b <- tr[folds == fo]; ya <- y[folds != fo]; yb <- y[folds == fo]
    Xa <- X[a, , drop = FALSE]
    f <- prep(Xa, ya, fam$binary)
    S <- if (fam$binary) Xa[, f$knn, drop = FALSE] else scale_with(Xa, Xa)
    nn <- knn_unlabelled(S, ya == 0, kmax, fam$binary)
    for (g in seq_len(nrow(ML$pu_grid))) {
      m <- pu_fit(Xa, ya, fam$binary, ML$pu_grid[g, ], nn, f, ML$seed + fo * 100 + g)
      p <- prob1(m, mx(X, b, f$cls))
      cvres[[length(cvres) + 1]] <- tibble(fold = fo, k = ML$pu_grid$k[g], t = ML$pu_grid$t[g],
                                           f1 = f1_pos(p >= 0.5, yb), avg_precision = avg_precision(p, yb),
                                           reliable_negatives = sum(reliable(nn, ML$pu_grid$k[g], ML$pu_grid$t[g])))
    }
    message(sprintf("  %-8s fold %d of %d done", nm, fo, ML$pu_folds))
  }
  cv <- bind_rows(cvres) |> group_by(k, t) |>
    summarise(cv_f1 = mean(f1), cv_f1_sd = sd(f1), cv_avg_precision = mean(avg_precision),
              reliable_negatives = round(mean(reliable_negatives)), .groups = "drop") |> arrange(desc(cv_f1))
  best <- cv[1, ]
  # final model on all training genes
  Xt <- X[tr, , drop = FALSE]
  f <- prep(Xt, y, fam$binary)
  S <- if (fam$binary) Xt[, f$knn, drop = FALSE] else scale_with(Xt, Xt)
  nn <- knn_unlabelled(S, y == 0, kmax, fam$binary)
  m <- pu_fit(Xt, y, fam$binary, best, nn, f, ML$seed)
  score <- rep(NA_real_, nrow(uni)); score[ev] <- prob1(m, mx(X, ev, f$cls))
  check_value(paste(nm, ": best k, t / CV F1 / reliable negatives (final)"),
              paste(paste(best$k, round(best$t, 3), sep = ", "), round(best$cv_f1, 3),
                    sum(reliable(nn, best$k, best$t)), sep = " / "))
  list(score = score, cv = mutate(cv, model = nm, .before = 1),
       perf = tibble(model = nm, training_genes = length(y), training_positives = sum(y), lambda = NA,
                     k = best$k, t = best$t, cv_avg_precision = best$cv_avg_precision, cv_f1 = best$cv_f1,
                     reliable_negatives = sum(reliable(nn, best$k, best$t))),
       features = tibble(model = nm, feature = f$cls))
}
log_step("M2: positive-unlabelled models")
pu_runs <- if (!is.null(cache)) cache$pu_runs else lapply(setNames(names(families), names(families)), run_pu)
if (is.null(cache)) saveRDS(list(m1_runs = m1_runs, pu_runs = pu_runs), cache_file)
# bind tables from the runs (unnamed list, so bind_rows stacks them as rows)
stack <- function(runs, what) bind_rows(unname(lapply(runs, `[[`, what)))

# ---- Output ------------------------------------------------------------------------------------
ev_genes <- filter(m1, set == "apply") |> select(gene, entrez, rule_selected, n_studies, n_sig, proportion_sig, replication_score)
scores <- ev_genes |>
  left_join(m1_runs[[1]]$scores, by = c("gene", "entrez")) |>
  left_join(m1_runs[[2]]$scores, by = c("gene", "entrez"))
for (nm in names(pu_runs)) scores[[nm]] <- pu_runs[[nm]]$score[match(scores$entrez, uni$entrez)]
write_csv(scores, file.path(out_dir, "scores_evaluation_genes.csv"), na = "")
perf <- bind_rows(stack(m1_runs, "perf") |> mutate(across(c(k, t, lambda), as.numeric), reliable_negatives = as.integer(reliable_negatives)),
                  stack(pu_runs, "perf") |> mutate(across(c(k, t, lambda), as.numeric), reliable_negatives = as.integer(reliable_negatives)))
write_csv(perf, file.path(out_dir, "model_training_summary.csv"), na = "")
write_csv(stack(pu_runs, "cv"), file.path(out_dir, "pu_cv_grid.csv"), na = "")
write_csv(stack(m1_runs, "coefs"), file.path(out_dir, "m1_coefficients.csv"), na = "")
write_csv(stack(pu_runs, "features"), file.path(out_dir, "pu_selected_features.csv"), na = "")

log_step("Summary")
print(as.data.frame(read_csv(file.path(out_dir, "model_training_summary.csv"), show_col_types = FALSE)), row.names = FALSE)
message("\nScores for the 213 genes written to ", file.path(out_dir, "scores_evaluation_genes.csv"),
        "\nNo outcome after 2017-06-01 has been used. Next: script 36 (comparison).")
