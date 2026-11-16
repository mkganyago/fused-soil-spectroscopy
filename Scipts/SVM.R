#SVM implementation for modelling soil properties

# 1: Load required libraries
suppressPackageStartupMessages({
  library(pls)
  library(caret)
  library(dplyr)
  library(ggplot2)
  library(corrplot)
  library(prospectr)
  library(readr)
  library(moments)
  library(kernlab) 
})

# Clear workspace Remove functions and libraries in work space. # Memory management
rm(list = ls())
set.seed(42)

# Helpers to avoid NA Coerces inputs to numeric, keeps only finite pairs usually NA is used by we use zero
.safe_cor <- function(x, y) {
  x <- as.numeric(x); y <- as.numeric(y)
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 2L) return(NA_real_)
  if (sd(x[ok]) == 0 || sd(y[ok]) == 0) return(0)
  suppressWarnings(stats::cor(x[ok], y[ok], use = "complete.obs"))
}
.ccc_fallback <- function(x, y) {
  x <- as.numeric(x); y <- as.numeric(y)
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]; y <- y[ok]
  if (length(x) < 2L) return(NA_real_)
  mx <- mean(x); my <- mean(y)
  vx <- stats::var(x); vy <- stats::var(y)
  sxy <- stats::cov(x, y)
  (2 * sxy) / (vx + vy + (mx - my)^2)
}
.compute_ccc <- function(obs, pred) {
  if (requireNamespace("epiR", quietly = TRUE)) {
    out <- tryCatch(epiR::epi.ccc(as.numeric(obs), as.numeric(pred))$rho.c["est"],
                    error = function(e) NA_real_)
    as.numeric(out)
  } else .ccc_fallback(obs, pred)
}
.compute_rpiq <- function(obs, rmse) {
  if (!is.finite(rmse) || rmse == 0) return(NA_real_)
  stats::IQR(as.numeric(obs), na.rm = TRUE) / rmse
}

#Checks skewness and then applies the appropriate transformation next time will try Box-Cox
check_normality <- function(y, alpha = 0.05) {
  sw_test <- shapiro.test(y)
  skew <- moments::skewness(y)
  kurt <- moments::kurtosis(y) - 3
  cat("Normality Assessment:\n")
  cat("Shapiro-Wilk p-value:", signif(sw_test$p.value, 4), "\n")
  cat("Skewness:", round(skew, 4), "\n")
  cat("Excess Kurtosis:", round(kurt, 4), "\n")
  is_normal <- sw_test$p.value > alpha && abs(skew) < 1
  if (!is_normal) {
    cat("Data is not normally distributed\n")
    if (abs(skew) > 1) cat("High skewness detected (|skew| > 1)\n")
  } else {
    cat("Data appears normally distributed\n")
  }
  list(
    is_normal = is_normal,
    skewness = skew,
    kurtosis = kurt,
    shapiro_p = sw_test$p.value,
    suggested_transform = suggest_transformation(skew)
  )
}
suggest_transformation <- function(skewness) {
  if (abs(skewness) < 0.5) "none"
  else if (skewness > 0.5) if (skewness > 1) "log" else "sqrt"
  else if (skewness < -1) "square" else "cube_root"
}
apply_transformation <- function(y, transform_type) {
  switch(transform_type,
         "none" = y,
         "log" = if (any(y <= 0)) { cat("Warning: Non-positive values found. Using log(y + 1)\n"); log(y + 1) } else log(y),
         "sqrt" = if (any(y < 0)) { cat("Warning: Negative values found. Using sqrt(y + min(y) + 1)\n"); sqrt(y + abs(min(y)) + 1) } else sqrt(y),
         "square" = y^2,
         "cube_root" = sign(y) * abs(y)^(1/3),
         y
  )
}
inverse_transformation <- function(y_pred, transform_type, original_y = NULL) {
  switch(transform_type,
         "none" = y_pred,
         "log" = if (!is.null(original_y) && any(original_y <= 0)) exp(y_pred) - 1 else exp(y_pred),
         "sqrt" = if (!is.null(original_y) && any(original_y < 0)) y_pred^2 - abs(min(original_y)) - 1 else y_pred^2,
         "square" = sqrt(y_pred),
         "cube_root" = sign(y_pred) * abs(y_pred)^3,
         y_pred
  )
}

#  Data loading
cat("Loading combined data...\n")
# sand_Wet_chem <- read_csv("C:\\Users\\aplenip\\Documents\\Noma\\FINAL DATA FOR MODELLING\\Clay Wet Chem.csv")
# MIR_averaged_noma <- read_csv("C:\\Users\\aplenip\\Documents\\Noma\\FINAL DATA FOR MODELLING\\VNIR_SWIR.csv")
# combined_data <- merge(sand_Wet_chem, MIR_averaged_noma, by = "SNM")
combined_data <- read_csv("C:\\Users\\aplenip\\Documents\\Gift\\final\\final_lead\\lead_full_duplicates.csv")
cat("Combined data dimensions:", dim(combined_data), "\n")
cat("Missing values in combined data:", sum(!is.finite(as.matrix(combined_data))), "\n")

id_col_name <- names(combined_data)[1]
cu_col_name <- names(combined_data)[2]
spectra_cols <- 3:ncol(combined_data)

spectral_data <- as.matrix(combined_data[, spectra_cols]); storage.mode(spectral_data) <- "numeric"
cu_data <- as.numeric(combined_data[[cu_col_name]])
ssn_data <- combined_data[[id_col_name]]

keep_rows <- is.finite(cu_data)
spectral_data <- spectral_data[keep_rows, , drop = FALSE]
cu_data <- cu_data[keep_rows]
ssn_data <- ssn_data[keep_rows]

cu_data_original <- cu_data

# Normality & transformation 
cat("\n", strrep("=", 60), "\n")
cat("NORMALITY ASSESSMENT AND TRANSFORMATION\n")
cat(strrep("=", 60), "\n")
normality_results <- check_normality(cu_data)
if (normality_results$suggested_transform != "none") {
  cat("\nApplying", normality_results$suggested_transform, "transformation...\n")
  cu_data_transformed <- apply_transformation(cu_data, normality_results$suggested_transform)
  cat("\nPost-transformation normality check:\n")
  normality_post <- check_normality(cu_data_transformed)
  if (abs(normality_post$skewness) < abs(normality_results$skewness)) {
    cu_data <- cu_data_transformed
    transformation_applied <- normality_results$suggested_transform
    cat("Transformation improved normality. Using transformed data.\n")
  } else {
    transformation_applied <- "none"
    cat("Transformation did not improve normality. Using original data.\n")
  }
} else {
  transformation_applied <- "none"
  cat("No transformation needed.\n")
}

# Wavelength vector (numeric if in names also removes X if data is imported using read.csv)
wl_all_chr <- colnames(spectral_data)
wl_all_num <- suppressWarnings(as.numeric(sub("^X", "", wl_all_chr)))

cat("\nTarget variable (Cu) statistics - Combined dataset:\n")
print(summary(cu_data))

# 80/20 representative split by lab data
cat("\nApplying representative sampling for 70/20 split based on Cu concentrations...\n")
N <- nrow(spectral_data)
train_size <- round(0.7 * N); val_size <- N - train_size
cat("Total samples:", N, "\n")
cat("Training samples:", train_size, "\n")
cat("Validation samples:", val_size, "\n")

ord <- order(cu_data)
n_val_per_interval <- max(1L, round(N / val_size))
val_idx <- integer(0)
for (i in seq(1, N, by = n_val_per_interval)) {
  j <- min(i + n_val_per_interval - 1L, N)
  pick <- if (i == j) i else sample.int(j - i + 1L, 1L) + (i - 1L)
  val_idx <- c(val_idx, ord[pick])
}
val_idx <- unique(val_idx)[1:val_size]
train_idx <- setdiff(seq_len(N), val_idx)

train_X <- spectral_data[train_idx, , drop = FALSE]
train_y <- cu_data[train_idx]
train_y_original <- cu_data_original[train_idx]
val_X <- spectral_data[val_idx, , drop = FALSE]
val_y <- cu_data[val_idx]
val_y_original <- cu_data_original[val_idx]

ks_test <- ks.test(train_y, val_y)

#Sminrnov test to check distibutions at 0.05
cat("Kolmogorov-Smirnov p:", signif(ks_test$p.value, 4), ifelse(ks_test$p.value > 0.05, "(similar)\n", "(different)\n"))


# Preprocessing helpers to avoid errors such as dividing by zero and finite values
.sanitize_pair <- function(Xtr, Xva, wl_num, nzv_thresh = 1e-12) {
  Xtr <- as.matrix(Xtr); Xva <- as.matrix(Xva)
  storage.mode(Xtr) <- "numeric"; storage.mode(Xva) <- "numeric"
  finite_cols <- colSums(is.finite(Xtr)) == nrow(Xtr) & colSums(is.finite(Xva)) == nrow(Xva)
  if (!any(finite_cols)) stop("All columns non-finite after initial check.")
  Xtr <- Xtr[, finite_cols, drop = FALSE]
  Xva <- Xva[, finite_cols, drop = FALSE]
  wl <- wl_num[finite_cols]
  sds <- apply(Xtr, 2, sd)
  keep_sd <- is.finite(sds) & (sds > nzv_thresh)
  if (!any(keep_sd)) stop("All columns had near-zero variance in training.")
  Xtr <- Xtr[, keep_sd, drop = FALSE]
  Xva <- Xva[, keep_sd, drop = FALSE]
  wl <- wl[keep_sd]
  nzv_idx <- tryCatch(caret::nearZeroVar(as.data.frame(Xtr), saveMetrics = FALSE), error = function(e) integer(0))
  if (length(nzv_idx)) {
    keep <- setdiff(seq_len(ncol(Xtr)), nzv_idx)
    if (length(keep) < 2L) stop("Too few columns after NZV removal.")
    Xtr <- Xtr[, keep, drop = FALSE]
    Xva <- Xva[, keep, drop = FALSE]
    wl <- wl[keep]
  }
  list(train = Xtr, val = Xva, wl = wl)
}
.snv_prospectr <- function(Xtr, Xva, wl_num) {
  Ztr <- prospectr::standardNormalVariate(as.matrix(Xtr))
  Zva <- prospectr::standardNormalVariate(as.matrix(Xva))
  .sanitize_pair(Ztr, Zva, wl_num)
}
.sg_prospectr <- function(Xtr, Xva, wl_num, w, p, deriv = 1) {
  Ztr <- prospectr::savitzkyGolay(as.matrix(Xtr), m = deriv, p = p, w = w)
  Zva <- prospectr::savitzkyGolay(as.matrix(Xva), m = deriv, p = p, w = w)
  .sanitize_pair(Ztr, Zva, wl_num)
}
.msc_prospectr <- function(Xtr, Xva, wl_num) {
  ref <- colMeans(Xtr, na.rm = TRUE)
  Ztr <- prospectr::msc(as.matrix(Xtr), ref_spectrum = ref)
  Zva <- prospectr::msc(as.matrix(Xva), ref_spectrum = ref)
  .sanitize_pair(Ztr, Zva, wl_num)
}

#Outlier detection & removal (data-level 3x rmse using PLS) 
detect_outliers <- function(observed, predicted, threshold_multiplier = 3) {
  residuals <- observed - predicted
  rmse <- sqrt(mean(residuals^2, na.rm = TRUE))
  threshold <- threshold_multiplier * rmse 
  outliers <- abs(residuals) > threshold
  outliers[is.na(outliers)] <- FALSE
  list(outliers = outliers, threshold = threshold, rmse = rmse,
       n_outliers = sum(outliers, na.rm = TRUE))
}
remove_outliers_iterative <- function(X, y, max_iterations = 3) {
  original_n <- nrow(X); idx <- seq_len(original_n); iteration <- 1; total_removed <- 0
  while (iteration <= max_iterations) {
    cat("Outlier detection iteration", iteration, "- Current samples:", nrow(X), "\n")
    temp_df <- data.frame(y = y, X)
    temp_model <- tryCatch(plsr(y ~ ., data = temp_df,
                                ncomp = min(10, nrow(X) - 1, ncol(X)),
                                center = TRUE, scale = TRUE, validation = "none"),
                           error = function(e) NULL)
    if (is.null(temp_model)) { cat("Could not fit model for outlier detection. Stopping.\n"); break }
    predictions <- as.numeric(predict(temp_model, newdata = temp_df,
                                      ncomp = temp_model$ncomp)[,,1])
    outlier_results <- detect_outliers(y, predictions)
    if (outlier_results$n_outliers == 0) { cat("No outliers detected. Stopping.\n"); break }
    cat("Detected", outlier_results$n_outliers, "outliers (threshold:",
        round(outlier_results$threshold, 4), ")\n")
    keep_idx <- !outlier_results$outliers
    X <- X[keep_idx, , drop = FALSE]; y <- y[keep_idx]; idx <- idx[keep_idx]
    total_removed <- total_removed + outlier_results$n_outliers
    iteration <- iteration + 1
  }
  cat("Total outliers removed:", total_removed, "out of", original_n,
      "(", round(100 * total_removed / original_n, 2), "%)\n")
  list(X = X, y = y, kept_idx = idx, removed_count = total_removed)
}

cat("\n", strrep("=", 60), "\n")
cat("OUTLIER DETECTION AND REMOVAL\n")
cat(strrep("=", 60), "\n")
outlier_results <- remove_outliers_iterative(train_X, train_y)
train_X_clean <- outlier_results$X
train_y_clean <- outlier_results$y
kept_idx <- outlier_results$kept_idx
train_y_original_clean <- train_y_original[kept_idx]

#SVM-RBF tuning
kfold_rmse <- function(X, y, sigma, C, k = 3, eps = 0.1) {
  folds <- caret::createFolds(y, k = k, list = TRUE, returnTrain = FALSE)
  se <- 0; n <- 0
  for (fold in folds) {
    X_te <- X[fold, , drop = FALSE]; y_te <- y[fold]
    X_tr <- X[-fold, , drop = FALSE]; y_tr <- y[-fold]
    fit <- tryCatch(
      kernlab::ksvm(x = X_tr, y = y_tr,
                    kernel = "rbfdot", kpar = list(sigma = sigma),
                    C = C, epsilon = eps, type = "eps-svr",
                    scaled = TRUE, prob.model = FALSE),
      error = function(e) NULL
    )
    if (is.null(fit)) next
    pred <- as.numeric(predict(fit, X_te))
    se <- se + sum((y_te - pred)^2, na.rm = TRUE)
    n  <- n  + length(y_te)
  }
  if (n == 0) return(Inf)
  sqrt(se / n)
}
tune_svm_rbf_fast <- function(X_train, y_train,
                              C_grid = c(0.5, 2, 8),
                              sigma_expand = c(0.5, 1, 2),
                              k = 3, eps = 0.1) {
  
  sig <- as.numeric(kernlab::sigest(X_train, frac = 0.5))[2]
  sigma_grid <- sig * sigma_expand
  grid <- expand.grid(C = C_grid, sigma = sigma_grid, KEEP.OUT.ATTRS = FALSE)
  grid$cv_rmse <- NA_real_
  best <- list(cv_rmse = Inf, C = NA_real_, sigma = NA_real_)
  for (i in seq_len(nrow(grid))) {
    C_i <- grid$C[i]; s_i <- grid$sigma[i]
    rm <- kfold_rmse(X_train, y_train, sigma = s_i, C = C_i, k = k, eps = eps)
    grid$cv_rmse[i] <- rm
    cat(sprintf("  SVM tune: C=%.3f, sigma=%.6g -> CV RMSE=%.4f\n", C_i, s_i, rm))
    if (is.finite(rm) && rm < best$cv_rmse) best <- list(cv_rmse = rm, C = C_i, sigma = s_i)
  }
  # final fit on full training
  final_model <- kernlab::ksvm(x = X_train, y = y_train,
                               kernel = "rbfdot", kpar = list(sigma = best$sigma),
                               C = best$C, epsilon = eps, type = "eps-svr",
                               scaled = TRUE, prob.model = FALSE)
  list(model = final_model, best = best, grid = grid)
}

calculate_metrics <- function(observed, predicted) {
  rmse <- sqrt(mean((observed - predicted)^2, na.rm = TRUE))
  r2 <- .safe_cor(observed, predicted)^2
  mae <- mean(abs(observed - predicted), na.rm = TRUE)
  ccc <- .compute_ccc(observed, predicted)
  rpiq <- .compute_rpiq(observed, rmse)
  list(RMSE = rmse, R2 = r2, MAE = mae, CCC = ccc, RPIQ = rpiq)
}

all_models <- list()
all_predictions <- list()
all_metrics <- list()

# MODEL 1: Base model
cat("\n", strrep("=", 60), "\n")
cat("TRAINING SVM-RBF MODELS\n")
cat(strrep("=", 60), "\n")

cat("Training Base SVM-RBF...\n")
base_pair <- .sanitize_pair(train_X_clean, val_X, wl_all_num)
base_model <- tune_svm_rbf_fast(base_pair$train, train_y_clean)
base_train_pred <- as.numeric(predict(base_model$model, base_pair$train))
base_val_pred   <- as.numeric(predict(base_model$model, base_pair$val))

if (transformation_applied != "none") {
  base_train_pred_orig <- inverse_transformation(base_train_pred, transformation_applied, train_y_original_clean)
  base_val_pred_orig   <- inverse_transformation(base_val_pred,   transformation_applied, val_y_original)
} else {
  base_train_pred_orig <- base_train_pred
  base_val_pred_orig   <- base_val_pred
}

base_train_metrics <- calculate_metrics(train_y_original_clean, base_train_pred_orig)
base_val_metrics   <- calculate_metrics(val_y_original, base_val_pred_orig)

all_models$base <- list(model = base_model$model, best = base_model$best, grid = base_model$grid, wl = base_pair$wl)
all_predictions$base <- list(train = base_train_pred_orig, val = base_val_pred_orig)
all_metrics$base <- list(train = base_train_metrics, val = base_val_metrics)

cat(sprintf("Base SVM-RBF - Best C=%.3f, sigma=%.6g; Val RMSE=%.4f, R²=%.4f\n",
            all_models$base$best$C, all_models$base$best$sigma,
            base_val_metrics$RMSE, base_val_metrics$R2))

# MODEL 2: SG (prospectr)
cat("\nTraining SVM-RBF + SG (prospectr)...\n")
sg_grid <- expand.grid(w = c(11, 15, 21), p = c(1, 2, 3))
best_sg <- NULL; best_sg_rmse <- Inf

for (g in seq_len(nrow(sg_grid))) {
  w <- sg_grid$w[g]; p <- sg_grid$p[g]
  sg_pair <- tryCatch(.sg_prospectr(train_X_clean, val_X, wl_all_num, w = w, p = p, deriv = 1),
                      error = function(e) NULL)
  if (is.null(sg_pair) || ncol(sg_pair$train) < 2) next
  sg_fit <- tune_svm_rbf_fast(sg_pair$train, train_y_clean)
  sg_val_pred <- as.numeric(predict(sg_fit$model, sg_pair$val))
  if (transformation_applied != "none") {
    sg_val_pred_orig <- inverse_transformation(sg_val_pred, transformation_applied, val_y_original)
    val_rmse <- sqrt(mean((val_y_original - sg_val_pred_orig)^2, na.rm = TRUE))
  } else {
    val_rmse <- sqrt(mean((val_y - sg_val_pred)^2, na.rm = TRUE))
  }
  if (is.finite(val_rmse) && val_rmse < best_sg_rmse) {
    best_sg_rmse <- val_rmse
    best_sg <- list(fit = sg_fit, w = w, p = p, pair = sg_pair)
  }
  cat(sprintf("  SG w=%d p=%d -> best C=%.3f sigma=%.6g | Val RMSE=%.4f\n",
              w, p, sg_fit$best$C, sg_fit$best$sigma, val_rmse))
}
if (is.null(best_sg)) stop("No valid SG configuration after sanitization.")
sg_train_pred <- as.numeric(predict(best_sg$fit$model, best_sg$pair$train))
sg_val_pred   <- as.numeric(predict(best_sg$fit$model, best_sg$pair$val))
if (transformation_applied != "none") {
  sg_train_pred_orig <- inverse_transformation(sg_train_pred, transformation_applied, train_y_original_clean)
  sg_val_pred_orig   <- inverse_transformation(sg_val_pred,   transformation_applied, val_y_original)
} else {
  sg_train_pred_orig <- sg_train_pred
  sg_val_pred_orig   <- sg_val_pred
}
sg_train_metrics <- calculate_metrics(train_y_original_clean, sg_train_pred_orig)
sg_val_metrics   <- calculate_metrics(val_y_original, sg_val_pred_orig)

all_models$sg <- list(model = best_sg$fit$model, best = best_sg$fit$best,
                      w = best_sg$w, p = best_sg$p, wl = best_sg$pair$wl)
all_predictions$sg <- list(train = sg_train_pred_orig, val = sg_val_pred_orig)
all_metrics$sg <- list(train = sg_train_metrics, val = sg_val_metrics)

cat(sprintf("Best SG config: w=%d, p=%d; Best C=%.3f, sigma=%.6g; Val RMSE=%.4f, R²=%.4f\n",
            all_models$sg$w, all_models$sg$p, all_models$sg$best$C, all_models$sg$best$sigma,
            sg_val_metrics$RMSE, sg_val_metrics$R2))

# MODEL 3: SNV (prospectr) 
cat("\nTraining SVM-RBF + SNV (prospectr)...\n")
snv_pair <- .snv_prospectr(train_X_clean, val_X, wl_all_num)
snv_fit <- tune_svm_rbf_fast(snv_pair$train, train_y_clean)
snv_train_pred <- as.numeric(predict(snv_fit$model, snv_pair$train))
snv_val_pred   <- as.numeric(predict(snv_fit$model, snv_pair$val))
if (transformation_applied != "none") {
  snv_train_pred_orig <- inverse_transformation(snv_train_pred, transformation_applied, train_y_original_clean)
  snv_val_pred_orig   <- inverse_transformation(snv_val_pred,   transformation_applied, val_y_original)
} else {
  snv_train_pred_orig <- snv_train_pred
  snv_val_pred_orig   <- snv_val_pred
}
snv_train_metrics <- calculate_metrics(train_y_original_clean, snv_train_pred_orig)
snv_val_metrics   <- calculate_metrics(val_y_original, snv_val_pred_orig)

all_models$snv <- list(model = snv_fit$model, best = snv_fit$best, wl = snv_pair$wl)
all_predictions$snv <- list(train = snv_train_pred_orig, val = snv_val_pred_orig)
all_metrics$snv <- list(train = snv_train_metrics, val = snv_val_metrics)

cat(sprintf("SNV SVM-RBF - Best C=%.3f, sigma=%.6g; Val RMSE=%.4f, R²=%.4f\n",
            all_models$snv$best$C, all_models$snv$best$sigma,
            snv_val_metrics$RMSE, snv_val_metrics$R2))

# MODEL 4: MSC (prospectr)
cat("\nTraining SVM-RBF + MSC (prospectr)...\n")
msc_pair <- .msc_prospectr(train_X_clean, val_X, wl_all_num)
msc_fit <- tune_svm_rbf_fast(msc_pair$train, train_y_clean)
msc_train_pred <- as.numeric(predict(msc_fit$model, msc_pair$train))
msc_val_pred   <- as.numeric(predict(msc_fit$model, msc_pair$val))
if (transformation_applied != "none") {
  msc_train_pred_orig <- inverse_transformation(msc_train_pred, transformation_applied, train_y_original_clean)
  msc_val_pred_orig   <- inverse_transformation(msc_val_pred,   transformation_applied, val_y_original)
} else {
  msc_train_pred_orig <- msc_train_pred
  msc_val_pred_orig   <- msc_val_pred
}
msc_train_metrics <- calculate_metrics(train_y_original_clean, msc_train_pred_orig)
msc_val_metrics   <- calculate_metrics(val_y_original, msc_val_pred_orig)

all_models$msc <- list(model = msc_fit$model, best = msc_fit$best, wl = msc_pair$wl)
all_predictions$msc <- list(train = msc_train_pred_orig, val = msc_val_pred_orig)
all_metrics$msc <- list(train = msc_train_metrics, val = msc_val_metrics)

cat(sprintf("MSC SVM-RBF - Best C=%.3f, sigma=%.6g; Val RMSE=%.4f, R²=%.4f\n",
            all_models$msc$best$C, all_models$msc$best$sigma,
            msc_val_metrics$RMSE, msc_val_metrics$R2))

# Model Summary 
model_summary <- data.frame(
  Model = c("Base SVM-RBF", "SVM-RBF + SG", "SVM-RBF + SNV", "SVM-RBF + MSC"),
  Best_C = c(all_models$base$best$C, all_models$sg$best$C, all_models$snv$best$C, all_models$msc$best$C),
  Best_sigma = c(all_models$base$best$sigma, all_models$sg$best$sigma, all_models$snv$best$sigma, all_models$msc$best$sigma),
  Train_RMSE = c(all_metrics$base$train$RMSE, all_metrics$sg$train$RMSE, all_metrics$snv$train$RMSE, all_metrics$msc$train$RMSE),
  Train_R2   = c(all_metrics$base$train$R2,   all_metrics$sg$train$R2,   all_metrics$snv$train$R2,   all_metrics$msc$train$R2),
  Val_RMSE   = c(all_metrics$base$val$RMSE,   all_metrics$sg$val$RMSE,   all_metrics$snv$val$RMSE,   all_metrics$msc$val$RMSE),
  Val_R2     = c(all_metrics$base$val$R2,     all_metrics$sg$val$R2,     all_metrics$snv$val$R2,     all_metrics$msc$val$R2),
  Val_CCC    = c(all_metrics$base$val$CCC,    all_metrics$sg$val$CCC,    all_metrics$snv$val$CCC,    all_metrics$msc$val$CCC),
  Val_RPIQ   = c(all_metrics$base$val$RPIQ,   all_metrics$sg$val$RPIQ,   all_metrics$snv$val$RPIQ,   all_metrics$msc$val$RPIQ),
  check.names = FALSE
)
model_summary[, c("Best_C","Best_sigma","Train_RMSE","Train_R2","Val_RMSE","Val_R2","Val_CCC","Val_RPIQ")] <-
  lapply(model_summary[, c("Best_C","Best_sigma","Train_RMSE","Train_R2","Val_RMSE","Val_R2","Val_CCC","Val_RPIQ")], function(v) round(as.numeric(v), 6))
print(model_summary)

results <-write_csv(model_summary, "Model summarry SVM.csv")

# Export model summary CSV
results_dir <- "results"; if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)
ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
summary_path <- file.path(results_dir, paste0("model_summary_SVM_RBF_", ts, ".csv"))
write.csv(model_summary, summary_path, row.names = FALSE)
cat("Model summary CSV saved to:", normalizePath(summary_path, winslash = "/"), "\n")

# Identify best model by Val_RMSE
best_model_idx <- which.min(model_summary$Val_RMSE)
cat("\nBest performing model:", model_summary$Model[best_model_idx], "\n")
cat("Best validation RMSE:", model_summary$Val_RMSE[best_model_idx], "\n")
cat("Best validation R²:", model_summary$Val_R2[best_model_idx], "\n")

# Combined plot (with model-outlier filterin 3x rmse) 
plot_colors <- list(
  train = list(base="#2ca02c", sg="#2ca02c", snv="#2ca02c", msc="#2ca02c"),
  val   = list(base="#ff7f0e", sg="#ff7f0e", snv="#ff7f0e", msc="#ff7f0e")
)
par_settings <- list(mar = c(4.5, 4.5, 3, 1), mgp = c(2.5, 1, 0), cex.main = 1.2, cex.lab = 1.1, cex.axis = 1.0)
model_outlier_k <- 3

make_combined_plot <- function(model_key, model_name, col_train, col_val,
                               pch_train = 16, pch_val = 17, legend_cex = 0.75) {
  train_obs  <- train_y_original_clean
  val_obs    <- val_y_original
  train_pred <- all_predictions[[model_key]]$train
  val_pred   <- all_predictions[[model_key]]$val
  tr_out <- detect_outliers(train_obs, train_pred, threshold_multiplier = model_outlier_k)
  va_out <- detect_outliers(val_obs,   val_pred,   threshold_multiplier = model_outlier_k)
  tr_mask <- !tr_out$outliers; va_mask <- !va_out$outliers
  met_tr <- calculate_metrics(train_obs[tr_mask], train_pred[tr_mask])
  met_va <- calculate_metrics(val_obs[va_mask],   val_pred[va_mask])
  rng <- range(c(train_obs[tr_mask], val_obs[va_mask], train_pred[tr_mask], val_pred[va_mask]), na.rm = TRUE)
  plot(train_obs[tr_mask], train_pred[tr_mask],
       xlim = rng, ylim = rng,
       pch = pch_train, cex = 0.80, col = col_train,
       main = model_name,
       xlab = paste("Actual", cu_col_name), ylab = paste("Predicted", cu_col_name),
       cex.main = par_settings$cex.main, cex.lab = par_settings$cex.lab)
  points(val_obs[va_mask], val_pred[va_mask], pch = pch_val, cex = 0.85, col = col_val)
  abline(0, 1, lwd = 2); grid()
  leg_metrics <- c(
    sprintf("Calibration: RMSE=%.3f  R²=%.3f  CCC=%.3f  RPIQ=%.3f", met_tr$RMSE, met_tr$R2, met_tr$CCC, met_tr$RPIQ),
    sprintf("Validation:  RMSE=%.3f  R²=%.3f  CCC=%.3f  RPIQ=%.3f", met_va$RMSE, met_va$R2, met_va$CCC, met_va$RPIQ)
  )
  legend("topleft", legend = leg_metrics, bty = "o", bg = rgb(1,1,1,0.92), cex = legend_cex)
  legend("bottomright", legend = c("Calibration", "Validation"),
         pch = c(pch_train, pch_val), pt.cex = c(0.90, 0.95),
         col = c(col_train, col_val), bty = "o", bg = rgb(1,1,1,0.92), cex = legend_cex)
}

cat("\nGenerating combined plots with model-outlier filtering (k = ", model_outlier_k, ")...\n", sep = "")
par(mfrow = c(2, 2), mar = par_settings$mar, mgp = par_settings$mgp)
make_combined_plot("base", "Base SVM-RBF", plot_colors$train$base, plot_colors$val$base)
make_combined_plot("sg",   "SVM-RBF + SG",  plot_colors$train$sg,   plot_colors$val$sg)
make_combined_plot("snv",  "SVM-RBF + SNV", plot_colors$train$snv,  plot_colors$val$snv)
make_combined_plot("msc",  "SVM-RBF + MSC", plot_colors$train$msc,  plot_colors$val$msc)

# SAVE COMBINED  
out_dir <- "plots_highres"; if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
save_png <- function(path, width = 3508, height = 2480, res = 300, plot_fun) {
  png(filename = path, width = width, height = height, res = res)
  old_par <- par(no.readonly = TRUE); on.exit({ par(old_par); dev.off() }, add = TRUE); plot_fun()
}
save_tiff <- function(path, width = 3508, height = 2480, res = 300, compression = "lzw", plot_fun) {
  tiff(filename = path, width = width, height = height, res = res, compression = compression)
  old_par <- par(no.readonly = TRUE); on.exit({ par(old_par); dev.off() }, add = TRUE); plot_fun()
}
plot_combined_grid_fun <- function() {
  par(mfrow = c(2, 2), mar = par_settings$mar, mgp = par_settings$mgp)
  make_combined_plot("base", "Base SVM-RBF", plot_colors$train$base, plot_colors$val$base)
  make_combined_plot("sg",   "SVM-RBF + SG",  plot_colors$train$sg,   plot_colors$val$sg)
  make_combined_plot("snv",  "SVM-RBF + SNV", plot_colors$train$snv,  plot_colors$val$snv)
  make_combined_plot("msc",  "SVM-RBF + MSC", plot_colors$train$msc,  plot_colors$val$msc)
}
grid_png_path  <- file.path(out_dir, paste0("combined_train_val_grid_SVM_RBF_", ts, ".png"))
grid_tiff_path <- file.path(out_dir, paste0("combined_train_val_grid_SVM_RBF_", ts, ".tiff"))
save_png (grid_png_path,  plot_fun = plot_combined_grid_fun)
save_tiff(grid_tiff_path, plot_fun = plot_combined_grid_fun)
cat("\nHigh-res combined plot saved to:\n  ",
    normalizePath(grid_png_path, winslash = "/"), "\n  ",
    normalizePath(grid_tiff_path, winslash = "/"), "\n", sep = "")

# Export predictions  
get_model_filtered_df <- function(model_key, model_name) {
  train_obs  <- train_y_original_clean
  val_obs    <- val_y_original
  train_pred <- all_predictions[[model_key]]$train
  val_pred   <- all_predictions[[model_key]]$val
  tr_out <- detect_outliers(train_obs, train_pred, threshold_multiplier = model_outlier_k)
  va_out <- detect_outliers(val_obs,   val_pred,   threshold_multiplier = model_outlier_k)
  tr_mask <- !tr_out$outliers; va_mask <- !va_out$outliers
  train_ids_clean <- ssn_data[train_idx][kept_idx]
  val_ids <- ssn_data[val_idx]
  df_tr <- data.frame(Model = model_name, Split = "Calibration",
                      ID = train_ids_clean[tr_mask], Actual = train_obs[tr_mask],
                      Predicted = train_pred[tr_mask], stringsAsFactors = FALSE)
  df_va <- data.frame(Model = model_name, Split = "Validation",
                      ID = val_ids[va_mask], Actual = val_obs[va_mask],
                      Predicted = val_pred[va_mask], stringsAsFactors = FALSE)
  list(all = rbind(df_tr, df_va), val = df_va)
}
df_base <- get_model_filtered_df("base", "Base SVM-RBF")
df_sg   <- get_model_filtered_df("sg",   "SVM-RBF + SG")
df_snv  <- get_model_filtered_df("snv",  "SVM-RBF + SNV")
df_msc  <- get_model_filtered_df("msc",  "SVM-RBF + MSC")

combined_filtered   <- dplyr::bind_rows(df_base$all, df_sg$all, df_snv$all, df_msc$all)
validation_filtered <- dplyr::bind_rows(df_base$val, df_sg$val, df_snv$val, df_msc$val)

if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)
combined_csv_path   <- file.path(results_dir, paste0("combined_predictions_filtered_SVM_RBF_", ts, ".csv"))
validation_csv_path <- file.path(results_dir, paste0("validation_actual_vs_predicted_filtered_SVM_RBF_", ts, ".csv"))
write.csv(combined_filtered, combined_csv_path, row.names = FALSE)
write.csv(validation_filtered, validation_csv_path, row.names = FALSE)
cat("\nFiltered combined predictions CSV saved to:\n  ",
    normalizePath(combined_csv_path, winslash = "/"), "\n", sep = "")
cat("Filtered Validation Actual vs Predicted CSV saved to:\n  ",
    normalizePath(validation_csv_path, winslash = "/"), "\n", sep = "")

# -------------------- Final Summary --------------------
cat("\n", strrep("=", 60), "\n")
cat("FINAL ANALYSIS SUMMARY (SVM-RBF)\n")
cat(strrep("=", 60), "\n")
cat("Data transformation applied:", transformation_applied, "\n")
cat("Total samples removed as outliers (data-level):", outlier_results$removed_count, "\n")
cat("Final training samples:", nrow(train_X_clean), "\n")
cat("Validation samples:", length(val_y_original), "\n")
cat("\nPreprocessing methods evaluated:\n")
cat("1. Base SVM-RBF\n2. SVM-RBF + Savitzky-Golay (1st deriv)\n3. SVM-RBF + SNV\n4. SVM-RBF + MSC\n")
cat("\nBest model by Val RMSE:", model_summary$Model[best_model_idx], "\n")

results_list <- list(
  models = all_models,
  predictions = all_predictions,
  metrics = all_metrics,
  summary_table = model_summary,
  transformation = transformation_applied,
  model_outlier_k = model_outlier_k,
  files = list(
    model_summary_csv = summary_path,
    combined_predictions_filtered_csv = combined_csv_path,
    validation_actual_vs_predicted_filtered_csv = validation_csv_path,
    combined_plot_png = grid_png_path,
    combined_plot_tiff = grid_tiff_path
  )
)

cat("\nAnalysis completed successfully!\n")
# =====================================================================
