#Monolithic version for running Soil spectroscopy models
#This is the plug and play version with all functions consolidated into this single script.
#This runs one property a time and this is Cu


# By Nomaxabiso Ngogodo :NgogodoN@arc.agric.za


# Load required libraries
suppressPackageStartupMessages({
  library(pls) 
  library(caret)
  library(dplyr)
  library(ggplot2)
  library(corrplot)
  library(prospectr)  
  library(readr)
  library(moments)   
})



# Clear workspace
rm(list = ls())
set.seed(42)

# -------------------- utilities --------------------

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

# -------------------- Normality and transformation functions --------------------
#Checks if data is Normal or not then applies the best transformation for that target property
#Here the bestNormalize function can also be used // 	https://petersonr.github.io/bestNormalize/,


check_normality <- function(y, alpha = 0.05) {
  sw_test <- shapiro.test(y)
  skew <- moments::skewness(y)
  kurt <- moments::kurtosis(y) - 3  # Threshold for extreme kurtosis
  
  cat("Normality Assessment:\n")
  cat("Shapiro-Wilk p-value:", signif(sw_test$p.value, 4), "\n")
  cat("Skewness:", round(skew, 4), "\n")
  cat("Excess Kurtosis:", round(kurt, 4), "\n")
  
  is_normal <- sw_test$p.value > alpha && abs(skew) < 1
  
  if (!is_normal) {
    cat("Data is not normally distributed\n")
    if (abs(skew) > 1) {
      cat("High skewness detected (|skew| > 1)\n")
    }
  } else {
    cat("Data appears normally distributed\n")
  }
  
  return(list(
    is_normal = is_normal,
    skewness = skew,
    kurtosis = kurt,
    shapiro_p = sw_test$p.value,
    suggested_transform = suggest_transformation(skew)
  ))
}



suggest_transformation <- function(skewness) {
  if (abs(skewness) < 0.5) {
    return("none")
  } else if (skewness > 0.5) {
    if (skewness > 1) {
      return("log")
    } else {
      return("sqrt")
    }
  } else {  # skewness < -0.5
    if (skewness < -1) {
      return("square")
    } else {
      return("cube_root")
    }
  }
}

apply_transformation <- function(y, transform_type) {
  switch(transform_type,
         "none" = y,
         "log" = {
           if (any(y <= 0)) {
             cat("Warning: Non-positive values found. Using log(y + 1)\n")
             log(y + 1)
           } else {
             log(y)
           }
         },
         "sqrt" = {
           if (any(y < 0)) {
             cat("Warning: Negative values found. Using sqrt(y + min(y) + 1)\n")
             sqrt(y + abs(min(y)) + 1)
           } else {
             sqrt(y)
           }
         },
         "square" = y^2,
         "cube_root" = sign(y) * abs(y)^(1/3),
         y  # default to no transformation
  )
}

inverse_transformation <- function(y_pred, transform_type, original_y = NULL) {
  switch(transform_type,
         "none" = y_pred,
         "log" = {
           if (!is.null(original_y) && any(original_y <= 0)) {
             exp(y_pred) - 1
           } else {
             exp(y_pred)
           }
         },
         "sqrt" = {
           if (!is.null(original_y) && any(original_y < 0)) {
             y_pred^2 - abs(min(original_y)) - 1
           } else {
             y_pred^2
           }
         },
         "square" = sqrt(y_pred),
         "cube_root" = sign(y_pred) * abs(y_pred)^3,
         y_pred  # default
  )
}


# -------------------- Data loading ----------------------------------------------
#Load data, make sure the data is in the working directory
cat("Loading combined data...\n")

cat("Loading combined data...\n")
sand_Wet_chem <- read_csv("Silt Wet Chem.csv")
MIR_averaged_noma <- read_csv("New_MIR_averaged_noma.csv")
AVAIL <- read_csv("available samples .csv")
MIR_averaged_noma <- merge(AVAIL, MIR_averaged_noma, by = "SNM")
combined_data <- merge(sand_Wet_chem, MIR_averaged_noma, by = "SNM")

cat("Combined data dimensions:", dim(combined_data), "\n")

cat("Missing values in combined data:", sum(!is.finite(as.matrix(combined_data))), "\n")

combined_data <- na.omit(combined_data)   
 
 cat("Missing values in combined data:", sum(!is.finite(as.matrix(combined_data))), "\n")
# -------------------- Columns & vectors --------------------
id_col_name <- names(combined_data)[1]
cu_col_name <- names(combined_data)[2]
spectra_cols <- 3:ncol(combined_data)

spectral_data <- as.matrix(combined_data[, spectra_cols])
storage.mode(spectral_data) <- "numeric"
cu_data <- as.numeric(combined_data[[cu_col_name]])
ssn_data <- combined_data[[id_col_name]]

# Remove rows with missing Cu
keep_rows <- is.finite(cu_data)
spectral_data <- spectral_data[keep_rows, , drop = FALSE]
cu_data <- cu_data[keep_rows]
ssn_data <- ssn_data[keep_rows]

# Store original Cu data
cu_data_original <- cu_data

# -------------------- Normality and transform (if helpful) --------------------
cat("\n", strrep("=", 60), "\n")
cat("NORMALITY ASSESSMENT AND TRANSFORMATION\n")
cat(strrep("=", 60), "\n")

normality_results <- check_normality(cu_data)

if (normality_results$suggested_transform != "none") {
  cat("\nApplying", normality_results$suggested_transform, "transformation...\n")
  cu_data_transformed <- apply_transformation(cu_data, normality_results$suggested_transform)
  
  # Check normality after transformation
  cat("\nPost-transformation normality check:\n")
  normality_post <- check_normality(cu_data_transformed)
  
  # Use transformed data if it improves normality
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

# Wavelength vector from column names
wl_all_chr <- colnames(spectral_data)
wl_all_num <- suppressWarnings(as.numeric(sub("^X", "", wl_all_chr)))

cat("\nTarget variable (Cu) statistics - Combined dataset:\n")
print(summary(cu_data))

# -------------------- 80/20 representative split by Cu --------------------
cat("\nApplying representative sampling for 80/20 split based on Cu concentrations...\n")
N <- nrow(spectral_data)
train_size <- round(0.8 * N)
val_size <- N - train_size
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
cat("Kolmogorov-Smirnov p:", signif(ks_test$p.value, 4), ifelse(ks_test$p.value > 0.05, "(similar)\n", "(different)\n"))

# -------------------- Preprocessing helpers --------------------
.sanitize_pair <- function(Xtr, Xva, wl_num, nzv_thresh = 1e-12) {
  Xtr <- as.matrix(Xtr); Xva <- as.matrix(Xva)
  storage.mode(Xtr) <- "numeric"; storage.mode(Xva) <- "numeric"
  
  finite_cols <- colSums(is.finite(Xtr)) == nrow(Xtr) & colSums(is.finite(Xva)) == nrow(Xva)
  if (!any(finite_cols)) stop("All columns non-finite after initial check.")
  Xtr <- Xtr[, finite_cols, drop = FALSE]
  Xva <- Xva[, finite_cols, drop = FALSE]
  wl <- wl_num[finite_cols]
  
  # Near-zero variance removal on TRAIN only
  sds <- apply(Xtr, 2, sd)
  keep_sd <- is.finite(sds) & (sds > nzv_thresh)
  if (!any(keep_sd)) stop("All columns had near-zero variance in training.")
  Xtr <- Xtr[, keep_sd, drop = FALSE]
  Xva <- Xva[, keep_sd, drop = FALSE]
  wl <- wl[keep_sd]
  
  # caret::nearZeroVar as extra guard
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

# Prospectr preprocessing functions
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

# Fourier self-deconvolution (for "deconv")
.msc_prospectr <- function(Xtr, Xva, wl_num) {
  ref <- colMeans(Xtr, na.rm = TRUE)
  Ztr <- prospectr::msc(as.matrix(Xtr), ref_spectrum = ref)
  Zva <- prospectr::msc(as.matrix(Xva), ref_spectrum = ref)
  .sanitize_pair(Ztr, Zva, wl_num)
}

# -------------------- Outlier detection and removal (data-level) --------------------
#adapted from Dangal et al., 2019
detect_outliers <- function(observed, predicted, threshold_multiplier = 5) {
  residuals <- observed - predicted
  rmse <- sqrt(mean(residuals^2, na.rm = TRUE))
  threshold <- threshold_multiplier * rmse
  outliers <- abs(residuals) > threshold                                                                                                   
  outliers[is.na(outliers)] <- FALSE
  return(list(
    outliers = outliers,
    threshold = threshold,
    rmse = rmse,
    n_outliers = sum(outliers, na.rm = TRUE)
  ))
}

remove_outliers_iterative <- function(X, y, max_iterations = 3) {
  original_n <- nrow(X)
  idx <- seq_len(original_n)  # Track mapping to original training rows
  iteration <- 1
  total_removed <- 0
  
  while (iteration <= max_iterations) {
    cat("Outlier detection iteration", iteration, "- Current samples:", nrow(X), "\n")
    
    # Fit simple PLS model for outlier detection
    temp_df <- data.frame(y = y, X)
    temp_model <- tryCatch({
      plsr(y ~ ., data = temp_df, ncomp = min(10, nrow(X) - 1, ncol(X)), 
           center = TRUE, scale = TRUE, validation = "none")
    }, error = function(e) NULL)
    
    if (is.null(temp_model)) {
      cat("Could not fit model for outlier detection. Stopping.\n")
      break
    }
    
    predictions <- as.numeric(predict(temp_model, newdata = temp_df, 
                                      ncomp = temp_model$ncomp)[,,1])
    
    outlier_results <- detect_outliers(y, predictions)
    
    if (outlier_results$n_outliers == 0) {
      cat("No outliers detected. Stopping.\n")
      break
    }
    
    cat("Detected", outlier_results$n_outliers, "outliers (threshold:", 
        round(outlier_results$threshold, 4), ")\n")
    
    # Remove outliers and keep index mapping
    keep_idx <- !outlier_results$outliers
    X <- X[keep_idx, , drop = FALSE]
    y <- y[keep_idx]
    idx <- idx[keep_idx]
    
    total_removed <- total_removed + outlier_results$n_outliers
    iteration <- iteration + 1
  }
  
  cat("Total outliers removed:", total_removed, "out of", original_n, 
      "(", round(100 * total_removed / original_n, 2), "%)\n")
  
  return(list(X = X, y = y, kept_idx = idx, removed_count = total_removed))
}

# Apply outlier removal to training data
cat("\n", strrep("=", 60), "\n")
cat("OUTLIER DETECTION AND REMOVAL\n")
cat(strrep("=", 60), "\n")

outlier_results <- remove_outliers_iterative(train_X, train_y)
train_X_clean <- outlier_results$X
train_y_clean <- outlier_results$y
kept_idx <- outlier_results$kept_idx
# Map back to original-scale target for the kept training rows
train_y_original_clean <- train_y_original[kept_idx]

# -------------------- Model optimization functions --------------------
optimize_pls_components <- function(X_train, y_train, X_val, y_val, max_comp = 20) {
  if (!is.matrix(X_train) || !is.matrix(X_val)) stop("optimize_pls_components expects matrices.")
  if (ncol(X_train) < 2L) stop("Too few predictors after sanitization (need >= 2).")
  
  max_components <- max(1L, min(max_comp, nrow(X_train) - 1L, ncol(X_train)))
  colnames(X_train) <- paste0("wl_", seq_len(ncol(X_train)))
  colnames(X_val) <- colnames(X_train)
  
  train_df <- data.frame(Cu = as.numeric(y_train), X_train, check.names = FALSE)
  comp_results <- data.frame(
    ncomp = 1:max_components,
    train_rmse = NA_real_, train_r2 = NA_real_,
    val_rmse = NA_real_, val_r2 = NA_real_
  )
  
  for (i in 1:max_components) {
    fit <- plsr(Cu ~ ., data = train_df, ncomp = i, center = TRUE, scale = TRUE,
                validation = "none", method = "oscorespls")
    pred_train <- as.numeric(predict(fit, newdata = as.data.frame(X_train), ncomp = i)[,,1])
    pred_val <- as.numeric(predict(fit, newdata = as.data.frame(X_val), ncomp = i)[,,1])
    
    comp_results$train_rmse[i] <- sqrt(mean((y_train - pred_train)^2, na.rm = TRUE))
    comp_results$train_r2[i] <- .safe_cor(y_train, pred_train)^2
    comp_results$val_rmse[i] <- sqrt(mean((y_val - pred_val)^2, na.rm = TRUE))
    comp_results$val_r2[i] <- .safe_cor(y_val, pred_val)^2
  }
  
  best_i <- which.min(comp_results$val_rmse)
  final_model <- plsr(Cu ~ ., data = train_df, ncomp = best_i, center = TRUE, scale = TRUE,
                      validation = "none", method = "oscorespls")
  list(model = final_model, optimal_comp = best_i, comp_results = comp_results,
       train_data = X_train, val_data = X_val)
}

calculate_metrics <- function(observed, predicted) {
  rmse <- sqrt(mean((observed - predicted)^2, na.rm = TRUE))
  r2 <- .safe_cor(observed, predicted)^2
  mae <- mean(abs(observed - predicted), na.rm = TRUE)
  ccc <- .compute_ccc(observed, predicted)
  rpiq <- .compute_rpiq(observed, rmse)
  list(RMSE = rmse, R2 = r2, MAE = mae, CCC = ccc, RPIQ = rpiq)
}

# -------------------- Storage --------------------
all_models <- list()
all_predictions <- list()
all_metrics <- list()

# -------------------- MODEL 1: Base (sanitize) --------------------
cat("\n", strrep("=", 60), "\n")
cat("TRAINING MODELS\n")
cat(strrep("=", 60), "\n")

cat("Training Base PLS model...\n")
base_pair <- .sanitize_pair(train_X_clean, val_X, wl_all_num)
base_model <- optimize_pls_components(base_pair$train, train_y_clean, base_pair$val, val_y)

base_train_pred <- as.numeric(predict(base_model$model, newdata = as.data.frame(base_model$train_data),
                                      ncomp = base_model$optimal_comp)[,,1])
base_val_pred <- as.numeric(predict(base_model$model, newdata = as.data.frame(base_model$val_data),
                                    ncomp = base_model$optimal_comp)[,,1])

# Transform predictions back if transformation was applied
if (transformation_applied != "none") {
  base_train_pred_orig <- inverse_transformation(base_train_pred, transformation_applied, train_y_original_clean)
  base_val_pred_orig <- inverse_transformation(base_val_pred, transformation_applied, val_y_original)
} else {
  base_train_pred_orig <- base_train_pred
  base_val_pred_orig <- base_val_pred
}

base_train_metrics <- calculate_metrics(train_y_original_clean, base_train_pred_orig)
base_val_metrics <- calculate_metrics(val_y_original, base_val_pred_orig)

all_models$base <- c(base_model, list(wl = base_pair$wl))
all_predictions$base <- list(train = base_train_pred_orig, val = base_val_pred_orig)
all_metrics$base <- list(train = base_train_metrics, val = base_val_metrics)

cat("Base PLS - Predictors:", ncol(base_model$train_data), "Optimal components:", base_model$optimal_comp, "\n")
cat("Base PLS - Validation RMSE:", round(base_val_metrics$RMSE, 4), "R²:", round(base_val_metrics$R2, 4), "\n")

# -------------------- MODEL 2: SG (prospectr) --------------------
cat("\nTraining PLS+SG models (prospectr)...\n")
sg_grid <- expand.grid(w = c(11, 15, 21), p = c(1, 2, 3))
best_sg <- NULL; best_sg_rmse <- Inf

for (g in seq_len(nrow(sg_grid))) {
  w <- sg_grid$w[g]; p <- sg_grid$p[g]
  sg_pair <- tryCatch(.sg_prospectr(train_X_clean, val_X, wl_all_num, w = w, p = p, deriv = 1),
                      error = function(e) NULL)
  if (is.null(sg_pair) || ncol(sg_pair$train) < 2) next
  
  sg_model <- optimize_pls_components(sg_pair$train, train_y_clean, sg_pair$val, val_y)
  sg_val_pred <- as.numeric(predict(sg_model$model, newdata = as.data.frame(sg_model$val_data),
                                    ncomp = sg_model$optimal_comp)[,,1])
  
  if (transformation_applied != "none") {
    sg_val_pred_orig <- inverse_transformation(sg_val_pred, transformation_applied, val_y_original)
    val_rmse <- sqrt(mean((val_y_original - sg_val_pred_orig)^2, na.rm = TRUE))
  } else {
    val_rmse <- sqrt(mean((val_y - sg_val_pred)^2, na.rm = TRUE))
  }
  
  if (is.finite(val_rmse) && val_rmse < best_sg_rmse) {
    best_sg_rmse <- val_rmse
    best_sg <- list(model = sg_model, w = w, p = p, wl = sg_pair$wl)
  }
  cat(sprintf("SG w=%d, p=%d: ncomp=%d, Val RMSE=%.4f, pcols=%d\n",
              w, p, sg_model$optimal_comp, val_rmse, ncol(sg_pair$train)))
}

if (is.null(best_sg)) stop("No valid SG configuration after sanitization.")

sg_train_pred <- as.numeric(predict(best_sg$model$model, newdata = as.data.frame(best_sg$model$train_data),
                                    ncomp = best_sg$model$optimal_comp)[,,1])
sg_val_pred <- as.numeric(predict(best_sg$model$model, newdata = as.data.frame(best_sg$model$val_data),
                                  ncomp = best_sg$model$optimal_comp)[,,1])

if (transformation_applied != "none") {
  sg_train_pred_orig <- inverse_transformation(sg_train_pred, transformation_applied, train_y_original_clean)
  sg_val_pred_orig <- inverse_transformation(sg_val_pred, transformation_applied, val_y_original)
} else {
  sg_train_pred_orig <- sg_train_pred
  sg_val_pred_orig <- sg_val_pred
}

sg_train_metrics <- calculate_metrics(train_y_original_clean, sg_train_pred_orig)
sg_val_metrics <- calculate_metrics(val_y_original, sg_val_pred_orig)

all_models$sg <- best_sg
all_predictions$sg <- list(train = sg_train_pred_orig, val = sg_val_pred_orig)
all_metrics$sg <- list(train = sg_train_metrics, val = sg_val_metrics)

cat("Best SG config: w =", best_sg$w, ", p =", best_sg$p, ", ncomp =", best_sg$model$optimal_comp, "\n")
cat("SG PLS - Validation RMSE:", round(sg_val_metrics$RMSE, 4), "R²:", round(sg_val_metrics$R2, 4), "\n")

# -------------------- MODEL 3: SNV (prospectr) --------------------
cat("\nTraining PLS+SNV model (prospectr)...\n")
snv_pair <- .snv_prospectr(train_X_clean, val_X, wl_all_num)
snv_model <- optimize_pls_components(snv_pair$train, train_y_clean, snv_pair$val, val_y)

snv_train_pred <- as.numeric(predict(snv_model$model, newdata = as.data.frame(snv_model$train_data),
                                     ncomp = snv_model$optimal_comp)[,,1])
snv_val_pred <- as.numeric(predict(snv_model$model, newdata = as.data.frame(snv_model$val_data),
                                   ncomp = snv_model$optimal_comp)[,,1])

if (transformation_applied != "none") {
  snv_train_pred_orig <- inverse_transformation(snv_train_pred, transformation_applied, train_y_original_clean)
  snv_val_pred_orig <- inverse_transformation(snv_val_pred, transformation_applied, val_y_original)
} else {
  snv_train_pred_orig <- snv_train_pred
  snv_val_pred_orig <- snv_val_pred
}

snv_train_metrics <- calculate_metrics(train_y_original_clean, snv_train_pred_orig)
snv_val_metrics <- calculate_metrics(val_y_original, snv_val_pred_orig)

all_models$snv <- c(snv_model, list(wl = snv_pair$wl))
all_predictions$snv <- list(train = snv_train_pred_orig, val = snv_val_pred_orig)
all_metrics$snv <- list(train = snv_train_metrics, val = snv_val_metrics)

cat("SNV PLS - Optimal components:", snv_model$optimal_comp, "\n")
cat("SNV PLS - Validation RMSE:", round(snv_val_metrics$RMSE, 4), "R²:", round(snv_val_metrics$R2, 4), "\n")

# -------------------- MODEL 4: MSC (prospectr) --------------------
cat("\nTraining PLS+MSC model (prospectr)...\n")
msc_pair <- .msc_prospectr(train_X_clean, val_X, wl_all_num)
msc_model <- optimize_pls_components(msc_pair$train, train_y_clean, msc_pair$val, val_y)

msc_train_pred <- as.numeric(predict(msc_model$model, newdata = as.data.frame(msc_model$train_data),
                                     ncomp = msc_model$optimal_comp)[,,1])
msc_val_pred <- as.numeric(predict(msc_model$model, newdata = as.data.frame(msc_model$val_data),
                                   ncomp = msc_model$optimal_comp)[,,1])

if (transformation_applied != "none") {
  msc_train_pred_orig <- inverse_transformation(msc_train_pred, transformation_applied, train_y_original_clean)
  msc_val_pred_orig <- inverse_transformation(msc_val_pred, transformation_applied, val_y_original)
} else {
  msc_train_pred_orig <- msc_train_pred
  msc_val_pred_orig <- msc_val_pred
}

msc_train_metrics <- calculate_metrics(train_y_original_clean, msc_train_pred_orig)
msc_val_metrics <- calculate_metrics(val_y_original, msc_val_pred_orig)

all_models$msc <- list(model = msc_model, wl = msc_pair$wl)
all_predictions$msc <- list(train = msc_train_pred_orig, val = msc_val_pred_orig)
all_metrics$msc <- list(train = msc_train_metrics, val = msc_val_metrics)

cat("MSC PLS - Optimal components:", msc_model$optimal_comp, "\n")
cat("MSC PLS - Validation RMSE:", round(msc_val_metrics$RMSE, 4), "R²:", round(msc_val_metrics$R2, 4), "\n")

# -------------------- Model Summary --------------------
print_header <- function(title) {
  cat("\n", strrep("=", 60), "\n")
  cat(title, "\n")
  cat(strrep("=", 60), "\n")
}

model_summary <- data.frame(
  Model = c("Base PLS", "PLS + SG", "PLS + SNV", "PLS + MSC"),
  Components = c(all_models$base$optimal_comp, all_models$sg$model$optimal_comp,
                 all_models$snv$optimal_comp, all_models$msc$model$optimal_comp),
  Train_RMSE = c(all_metrics$base$train$RMSE, all_metrics$sg$train$RMSE,
                 all_metrics$snv$train$RMSE, all_metrics$msc$train$RMSE),
  Train_R2 = c(all_metrics$base$train$R2, all_metrics$sg$train$R2,
               all_metrics$snv$train$R2, all_metrics$msc$train$R2),
  Val_RMSE = c(all_metrics$base$val$RMSE, all_metrics$sg$val$RMSE,
               all_metrics$snv$val$RMSE, all_metrics$msc$val$RMSE),
  Val_R2 = c(all_metrics$base$val$R2, all_metrics$sg$val$R2,
             all_metrics$snv$val$R2, all_metrics$msc$val$R2),
  Val_CCC = c(all_metrics$base$val$CCC, all_metrics$sg$val$CCC,
              all_metrics$snv$val$CCC, all_metrics$msc$val$CCC),
  Val_RPIQ = c(all_metrics$base$val$RPIQ, all_metrics$sg$val$RPIQ,
               all_metrics$snv$val$RPIQ, all_metrics$msc$val$RPIQ)
)

# Round numerical columns
model_summary[, 3:8] <- round(model_summary[, 3:8], 4)
print(model_summary)

# -------------------- Export model summary CSV --------------------
results_dir <- "results"
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)
ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
summary_path <- file.path(results_dir, paste0("model_summary_", ts, ".csv"))
write.csv(model_summary, summary_path, row.names = FALSE)
cat("Model summary CSV saved to:", normalizePath(summary_path, winslash = "/"), "\n")

# Identify best model
best_model_idx <- which.min(model_summary$Val_RMSE)
cat("\nBest performing model:", model_summary$Model[best_model_idx], "\n")
cat("Best validation RMSE:", model_summary$Val_RMSE[best_model_idx], "\n")
cat("Best validation R²:", model_summary$Val_R2[best_model_idx], "\n")

# -------------------- Plotting (existing curves) --------------------
plot_colors <- list(
  train = list(
    base = "#2ca02c",
    sg   = "#2ca02c", 
    snv  = "#2ca02c",
    msc  = "#2ca02c"
  ),
  val = list(
    base = "#ff7f0e",
    sg   = "#ff7f0e",
    snv  = "#ff7f0e", 
    msc  = "#ff7f0e"
  )
)
par_settings <- list(mar = c(4.5, 4.5, 3, 1), mgp = c(2.5, 1, 0), cex.main = 1.2, cex.lab = 1.1, cex.axis = 1.0)

# Component optimization curves (display only; not exported individually)
par(mfrow = c(2, 3), mar = par_settings$mar, mgp = par_settings$mgp)

# Training RMSE vs Components
y_range_train_rmse <- range(c(all_models$base$comp_results$train_rmse,
                              all_models$sg$model$comp_results$train_rmse,
                              all_models$snv$comp_results$train_rmse,
                              all_models$msc$model$comp_results$train_rmse), na.rm = TRUE)

plot(all_models$base$comp_results$ncomp, all_models$base$comp_results$train_rmse,
     type="l", col=plot_colors$train$base, lwd=2, main="Training RMSE vs Components",
     xlab="Components", ylab="RMSE", ylim=y_range_train_rmse,
     cex.main=par_settings$cex.main, cex.lab=par_settings$cex.lab)
lines(all_models$sg$model$comp_results$ncomp, all_models$sg$model$comp_results$train_rmse, 
      col=plot_colors$train$sg, lwd=2)
lines(all_models$snv$comp_results$ncomp, all_models$snv$comp_results$train_rmse, 
      col=plot_colors$train$snv, lwd=2)
lines(all_models$msc$model$comp_results$ncomp, all_models$msc$model$comp_results$train_rmse, 
      col=plot_colors$train$msc, lwd=2)
legend("topright", c("Base","SG","SNV","MSC"), col=unlist(plot_colors$train), lwd=2, cex=0.8)
grid()

# Validation RMSE vs Components
y_range_val_rmse <- range(c(all_models$base$comp_results$val_rmse,
                            all_models$sg$model$comp_results$val_rmse,
                            all_models$snv$comp_results$val_rmse,
                            all_models$msc$model$comp_results$val_rmse), na.rm = TRUE)

plot(all_models$base$comp_results$ncomp, all_models$base$comp_results$val_rmse,
     type="l", col=plot_colors$val$base, lwd=2, main="Validation RMSE vs Components",
     xlab="Components", ylab="RMSE", ylim=y_range_val_rmse,
     cex.main=par_settings$cex.main, cex.lab=par_settings$cex.lab)
lines(all_models$sg$model$comp_results$ncomp, all_models$sg$model$comp_results$val_rmse, 
      col=plot_colors$val$sg, lwd=2)
lines(all_models$snv$comp_results$ncomp, all_models$snv$comp_results$val_rmse, 
      col=plot_colors$val$snv, lwd=2)
lines(all_models$msc$model$comp_results$ncomp, all_models$msc$model$comp_results$val_rmse, 
      col=plot_colors$val$msc, lwd=2)
legend("topright", c("Base","SG","SNV","MSC"), col=unlist(plot_colors$val), lwd=2, cex=0.8)
grid()

# Training R² vs Components
y_range_train_r2 <- range(c(all_models$base$comp_results$train_r2,
                            all_models$sg$model$comp_results$train_r2,
                            all_models$snv$comp_results$train_r2,
                            all_models$msc$model$comp_results$train_r2), na.rm = TRUE)

plot(all_models$base$comp_results$ncomp, all_models$base$comp_results$train_r2,
     type="l", col=plot_colors$train$base, lwd=2, main="Training R² vs Components",
     xlab="Components", ylab="R²", ylim=y_range_train_r2,
     cex.main=par_settings$cex.main, cex.lab=par_settings$cex.lab)
lines(all_models$sg$model$comp_results$ncomp, all_models$sg$model$comp_results$train_r2, 
      col=plot_colors$train$sg, lwd=2)
lines(all_models$snv$comp_results$ncomp, all_models$snv$comp_results$train_r2, 
      col=plot_colors$train$snv, lwd=2)
lines(all_models$msc$model$comp_results$ncomp, all_models$msc$model$comp_results$train_r2, 
      col=plot_colors$train$msc, lwd=2)
grid()

# Validation R² vs Components
y_range_val_r2 <- range(c(all_models$base$comp_results$val_r2,
                          all_models$sg$model$comp_results$val_r2,
                          all_models$snv$comp_results$val_r2,
                          all_models$msc$model$comp_results$val_r2), na.rm = TRUE)

plot(all_models$base$comp_results$ncomp, all_models$base$comp_results$val_r2,
     type="l", col=plot_colors$val$base, lwd=2, main="Validation R² vs Components",
     xlab="Components", ylab="R²", ylim=y_range_val_r2,
     cex.main=par_settings$cex.main, cex.lab=par_settings$cex.lab)
lines(all_models$sg$model$comp_results$ncomp, all_models$sg$model$comp_results$val_r2, 
      col=plot_colors$val$sg, lwd=2)
lines(all_models$snv$comp_results$ncomp, all_models$snv$comp_results$val_r2, 
      col=plot_colors$val$snv, lwd=2)
lines(all_models$msc$model$comp_results$ncomp, all_models$msc$model$comp_results$val_r2, 
      col=plot_colors$val$msc, lwd=2)
grid()

# -------------------- Combined plot (with model-outlier filtering) --------------------
model_outlier_k <- 3  # threshold multiplier; set 2.0–2.5 for stricter

make_combined_plot <- function(model_key, model_name, col_train, col_val,
                               pch_train = 16,  # filled circle = Calibration/Train
                               pch_val   = 17,  # triangle = Validation
                               legend_cex = 0.75) {
  # Gather observed & predicted (original scale)
  train_obs  <- train_y_original_clean
  val_obs    <- val_y_original
  train_pred <- all_predictions[[model_key]]$train
  val_pred   <- all_predictions[[model_key]]$val
  
  # Model-outlier detection per split
  tr_out <- detect_outliers(train_obs, train_pred, threshold_multiplier = model_outlier_k)
  va_out <- detect_outliers(val_obs,   val_pred,   threshold_multiplier = model_outlier_k)
  tr_mask <- !tr_out$outliers
  va_mask <- !va_out$outliers
  
  # Metrics AFTER model-outlier removal
  met_tr <- calculate_metrics(train_obs[tr_mask], train_pred[tr_mask])
  met_va <- calculate_metrics(val_obs[va_mask],   val_pred[va_mask])
  
  # Plot both sets together with distinct symbols
  rng <- range(c(train_obs[tr_mask], val_obs[va_mask],
                 train_pred[tr_mask], val_pred[va_mask]), na.rm = TRUE)
  
  plot(train_obs[tr_mask], train_pred[tr_mask],
       xlim = rng, ylim = rng,
       pch = pch_train, cex = 0.80, col = col_train,
       main = model_name,  # simplified title as requested
       xlab = paste("Actual", cu_col_name), ylab = paste("Predicted", cu_col_name),
       cex.main = par_settings$cex.main, cex.lab = par_settings$cex.lab)
  points(val_obs[va_mask], val_pred[va_mask],
         pch = pch_val, cex = 0.85, col = col_val)
  
  abline(0, 1, lwd = 2)
  grid()
  
  # Metrics legend
  leg_metrics <- c(
    sprintf("Calibration: RMSE=%.3f  R²=%.3f  CCC=%.3f  RPIQ=%.3f",
            met_tr$RMSE, met_tr$R2, met_tr$CCC, met_tr$RPIQ),
    sprintf("Validation:  RMSE=%.3f  R²=%.3f  CCC=%.3f  RPIQ=%.3f",
            met_va$RMSE, met_va$R2, met_va$CCC, met_va$RPIQ)
  )
  legend("topleft", legend = leg_metrics, bty = "o",
         bg = rgb(1,1,1,0.92), cex = legend_cex)
  
  # Symbol legend: BOTH labels, no counts
  legend("bottomright",
         legend = c("Calibration", "Validation"),
         pch    = c(pch_train, pch_val),
         pt.cex = c(0.90, 0.95),
         col    = c(col_train, col_val),
         bty    = "o",
         bg     = rgb(1, 1, 1, 0.92),
         cex    = legend_cex)
}

cat("\nGenerating combined plots with model-outlier filtering (k = ",
    model_outlier_k, ")...\n", sep = "")
par(mfrow = c(2, 2), mar = par_settings$mar, mgp = par_settings$mgp)
make_combined_plot("base", " Base PLS", plot_colors$train$base, plot_colors$val$base)
make_combined_plot("sg",   "PLS + SG",  plot_colors$train$sg,   plot_colors$val$sg)
make_combined_plot("snv",  " PLS + SNV", plot_colors$train$snv,  plot_colors$val$snv)
make_combined_plot("msc",  "PLS + MSC", plot_colors$train$msc,  plot_colors$val$msc)

# -------------------- SAVE COMBINED GRID ONLY (HIGH-RES IMAGE) --------------------
out_dir <- "plots_highres"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

save_png <- function(path, width = 3508, height = 2480, res = 300, plot_fun) {
  png(filename = path, width = width, height = height, res = res)
  old_par <- par(no.readonly = TRUE)
  on.exit({ par(old_par); dev.off() }, add = TRUE)
  plot_fun()
}
save_tiff <- function(path, width = 3508, height = 2480, res = 300, compression = "lzw", plot_fun) {
  tiff(filename = path, width = width, height = height, res = res, compression = compression)
  old_par <- par(no.readonly = TRUE)
  on.exit({ par(old_par); dev.off() }, add = TRUE)
  plot_fun()
}

plot_combined_grid_fun <- function() {
  par(mfrow = c(2, 2), mar = par_settings$mar, mgp = par_settings$mgp)
  make_combined_plot("base", "Base PLS", plot_colors$train$base, plot_colors$val$base)
  make_combined_plot("sg",   "PLS + SG",  plot_colors$train$sg,   plot_colors$val$sg)
  make_combined_plot("snv",  "PLS + SNV", plot_colors$train$snv,  plot_colors$val$snv)
  make_combined_plot("msc",  "PLS + MSC", plot_colors$train$msc,  plot_colors$val$msc)
}

# Save the 2x2 combined grid (A4 landscape ~ 3508x2480 px @ 300 dpi)
grid_png_path  <- file.path(out_dir, paste0("combined_train_val_grid_", ts, ".png"))
grid_tiff_path <- file.path(out_dir, paste0("combined_train_val_grid_", ts, ".tiff"))
save_png (grid_png_path,  plot_fun = plot_combined_grid_fun)
save_tiff(grid_tiff_path, plot_fun = plot_combined_grid_fun)
cat("\nHigh-res combined plot saved to:\n  ",
    normalizePath(grid_png_path, winslash = "/"), "\n  ",
    normalizePath(grid_tiff_path, winslash = "/"), "\n", sep = "")

# -------------------- NEW: Export predictions with MODEL-OUTLIERS REMOVED --------------------
# Helper to build filtered data frames for each model
get_model_filtered_df <- function(model_key, model_name) {
  train_obs  <- train_y_original_clean
  val_obs    <- val_y_original
  train_pred <- all_predictions[[model_key]]$train
  val_pred   <- all_predictions[[model_key]]$val
  
  tr_out <- detect_outliers(train_obs, train_pred, threshold_multiplier = model_outlier_k)
  va_out <- detect_outliers(val_obs,   val_pred,   threshold_multiplier = model_outlier_k)
  tr_mask <- !tr_out$outliers
  va_mask <- !va_out$outliers
  
  # IDs for rows after data-level outlier removal (training) and original val
  train_ids_clean <- ssn_data[train_idx][kept_idx]
  val_ids <- ssn_data[val_idx]
  
  df_tr <- data.frame(
    Model = model_name,
    Split = "Calibration",
    ID = train_ids_clean[tr_mask],
    Actual = train_obs[tr_mask],
    Predicted = train_pred[tr_mask],
    stringsAsFactors = FALSE
  )
  df_va <- data.frame(
    Model = model_name,
    Split = "Validation",
    ID = val_ids[va_mask],
    Actual = val_obs[va_mask],
    Predicted = val_pred[va_mask],
    stringsAsFactors = FALSE
  )
  list(all = rbind(df_tr, df_va), val = df_va)
}

# Build combined and validation-only (filtered) data frames across all models
df_base <- get_model_filtered_df("base", "Base PLS")
df_sg   <- get_model_filtered_df("sg",   "PLS + SG")
df_snv  <- get_model_filtered_df("snv",  "PLS + SNV")
df_msc  <- get_model_filtered_df("msc",  "PLS + MSC")

combined_filtered <- dplyr::bind_rows(df_base$all, df_sg$all, df_snv$all, df_msc$all)
validation_filtered <- dplyr::bind_rows(df_base$val, df_sg$val, df_snv$val, df_msc$val)

# Export CSVs
combined_csv_path <- file.path(results_dir, paste0("combined_predictions_filtered_", ts, ".csv"))
validation_csv_path <- file.path(results_dir, paste0("validation_actual_vs_predicted_filtered_", ts, ".csv"))
write.csv(combined_filtered, combined_csv_path, row.names = FALSE)
write.csv(validation_filtered, validation_csv_path, row.names = FALSE)

cat("\nFiltered combined predictions CSV (Calibration + Validation) saved to:\n  ",
    normalizePath(combined_csv_path, winslash = "/"), "\n", sep = "")
cat("Filtered Validation Actual vs Predicted CSV saved to:\n  ",
    normalizePath(validation_csv_path, winslash = "/"), "\n", sep = "")

# -------------------- Final Summary --------------------
cat("\n", strrep("=", 60), "\n")
cat("FINAL ANALYSIS SUMMARY\n")
cat(strrep("=", 60), "\n")

cat("Data transformation applied:", transformation_applied, "\n")
cat("Total samples removed as outliers (data-level):", outlier_results$removed_count, "\n")
cat("Final training samples:", nrow(train_X_clean), "\n")
cat("Validation samples:", length(val_y_original), "\n")

cat("\nPreprocessing methods evaluated:\n")
cat("1. Base PLS (no preprocessing)\n")
cat("2. PLS + Savitzky-Golay smoothing + 1st derivative\n")
cat("3. PLS + Standard Normal Variate (SNV)\n") 
cat("4. PLS + Multiplicative Scatter Correction (MSC)\n")

cat("\nBest SG parameters: w =", all_models$sg$w, ", p =", all_models$sg$p, "\n")

# Optional container object for downstream use
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
