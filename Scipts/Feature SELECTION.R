library(tidyverse)
library(prospectr)
library(caret)
library(bestNormalize)
library(DescTools)
library(pls)
library(gridExtra)

set.seed(123)

# -----------------------
# 0. PARAMETERS
# -----------------------
cv_folds <- 10
plsr_max_comps <- 30
sg_grid <- expand.grid(p = c(1,2,3), w = seq(11, 25, 2), m = c(0,1,2))
max_outlier_pct <- 0.01

# CARS parameters
cars_sampling_num <- 100
cars_cv_folds <- 10
cars_max_iter <- 500

# Properties to model
properties <- c("Clay", "Cu", "Mn", "pH", "Sand", "Se", "Silt", "Zn")

# -----------------------
# HELPER: SAFE DIMENSION & DIMNAMES RESTORATION
# -----------------------
restore_spectral_dimnames <- function(mat, ref_mat) {
  if (!is.matrix(mat)) mat <- as.matrix(mat)
  if (length(dim(mat)) != 2) stop("Input must be 2D matrix")
  
  target_nrow <- nrow(ref_mat)
  target_ncol <- ncol(ref_mat)
  
  # Adjust rows
  if (nrow(mat) < target_nrow) {
    repeat_rows <- target_nrow - nrow(mat)
    mat <- rbind(mat, mat[rep(nrow(mat), repeat_rows), , drop = FALSE])
  } else if (nrow(mat) > target_nrow) {
    mat <- mat[1:target_nrow, , drop = FALSE]
  }
  
  # Adjust columns
  if (ncol(mat) < target_ncol) {
    repeat_cols <- target_ncol - ncol(mat)
    mat <- cbind(mat, mat[, rep(ncol(mat), repeat_cols), drop = FALSE])
  } else if (ncol(mat) > target_ncol) {
    mat <- mat[, 1:target_ncol, drop = FALSE]
  }
  
  rownames(mat) <- rownames(ref_mat)
  colnames(mat) <- colnames(ref_mat)
  return(mat)
}

# -----------------------
# 1. LOAD SPECTRA (COMMON FOR ALL)
# -----------------------
cat("\n=== LOADING SPECTRAL DATA ===\n")
spectra <- read.csv("combined_data1.csv", stringsAsFactors = FALSE)

spectral_cols <- names(spectra)[names(spectra) != "SNM"]
wavenumbers <- suppressWarnings(as.numeric(gsub("[^0-9.]", "", spectral_cols)))
if (all(is.na(wavenumbers))) wavenumbers <- NULL

if (!is.null(wavenumbers)) {
  co2_region <- which(wavenumbers >= 2268 & wavenumbers <= 2389)
  if (length(co2_region) > 0) {
    keep_cols <- c("SNM", spectral_cols[-co2_region])
    spectra <- spectra[, keep_cols, drop = FALSE]
    wavenumbers <- wavenumbers[-co2_region]
    spectral_cols <- names(spectra)[names(spectra) != "SNM"]
  }
}

cat("Spectral data loaded: n =", nrow(spectra), "features =", length(spectral_cols), "\n")

# -----------------------
# 2. SPECTRAL QC (ONE-TIME, COMMON FOR ALL)
# -----------------------
cat("\n=== SPECTRAL QC: PCA MAHALANOBIS ===\n")
spectral_mat <- as.matrix(spectra[, spectral_cols])
rownames(spectral_mat) <- spectra$SNM

if (!is.null(wavenumbers)) {
  target_wav <- 500
  wav_idx <- which.min(abs(wavenumbers - target_wav))
  bad_samples <- which(spectral_mat[, wav_idx] < 1)
  if (length(bad_samples) > 0) {
    cat("Removing", length(bad_samples), "low-reflectance samples at ~", target_wav, "cm⁻¹\n")
    spectral_mat <- spectral_mat[-bad_samples, , drop = FALSE]
    rownames(spectral_mat) <- spectra$SNM[-bad_samples]
  }
}

X <- scale(spectral_mat, center = TRUE, scale = TRUE)
pca <- prcomp(X, center = FALSE, scale. = FALSE)
var_explained <- cumsum(pca$sdev^2) / sum(pca$sdev^2)
n_pcs <- min(which(var_explained >= 0.95))
n_pcs <- max(3, min(n_pcs, 20))
scores <- pca$x[, 1:n_pcs, drop = FALSE]

alpha_md <- 0.995
md_cutoff <- qchisq(alpha_md, df = n_pcs)
md <- mahalanobis(scores, center = colMeans(scores), cov = cov(scores))
spec_outliers <- which(md > md_cutoff)

cat(sprintf("PCA MD: %d PCs, p=%.3f -> flagged %d samples\n", n_pcs, alpha_md, length(spec_outliers)))

if (length(spec_outliers) > 0) {
  spectral_mat <- spectral_mat[-spec_outliers, , drop = FALSE]
  rownames(spectral_mat) <- spectra$SNM[-spec_outliers]
}

cat("After spectral QC: n =", nrow(spectral_mat), "\n")

spectra_clean <- as_tibble(spectral_mat) %>%
  mutate(SNM = rownames(spectral_mat)) %>%
  relocate(SNM)

# -----------------------
# 3. FUNCTION: DANGAL SINGLE-PASS OUTLIER REMOVAL
# -----------------------
dangal_outlier_removal <- function(spectral_data, lab_data, property_name, max_pct = 0.01) {
  cat("\n--- Dangal outlier removal for", property_name, "---\n")
  
  lab_data_clean <- lab_data %>% select(SNM, Lab = all_of(property_name))
  spec_data_clean <- spectral_data
  
  lab_data_unique <- lab_data_clean %>% distinct(SNM, .keep_all = TRUE)
  spectra_unique <- spec_data_clean %>% distinct(SNM, .keep_all = TRUE)
  
  data_merged <- inner_join(lab_data_unique, spectra_unique, by = "SNM") %>% na.omit()
  n_initial <- nrow(data_merged)
  cat("Initial samples:", n_initial, "\n")
  
  if (n_initial < 50) {
    cat("WARNING: Too few samples (", n_initial, "). Skipping outlier removal.\n")
    return(list(data = data_merged, removed = tibble(), n_removed = 0))
  }
  
  X_spec <- as.matrix(data_merged %>% select(-SNM, -Lab))
  y_lab <- data_merged$Lab
  
  bn_obj <- bestNormalize(y_lab, allow_orderNorm = FALSE, standardize = FALSE)
  y_transformed <- predict(bn_obj)
  
  transform_name <- if (!is.null(bn_obj$chosen_transform)) class(bn_obj$chosen_transform)[1] else class(bn_obj)[1]
  
  calc_skew <- function(x) {
    n <- length(x)
    x_centered <- x - mean(x)
    m3 <- sum(x_centered^3) / n
    m2 <- sum(x_centered^2) / n
    return(m3 / (m2^(3/2)))
  }
  
  orig_skew <- calc_skew(y_lab)
  trans_skew <- calc_skew(y_transformed)
  cat(sprintf("Best transformation selected: %s | Skew: %.3f → %.3f\n", 
              transform_name, orig_skew, trans_skew))
  
  png(paste0(property_name, "_transformation_diagnostic.png"), width = 1200, height = 500, res = 100)
  par(mfrow = c(1, 2))
  hist(y_lab, breaks = 30, main = paste(property_name, "- Original"), 
       xlab = property_name, col = "lightblue", border = "white")
  text(x = quantile(y_lab, 0.95), y = par("usr")[4] * 0.9, 
       labels = sprintf("Skew: %.2f", orig_skew), pos = 2)
  hist(y_transformed, breaks = 30, main = paste(property_name, "- Transformed"), 
       xlab = paste("Transformed", property_name), col = "lightgreen", border = "white")
  text(x = quantile(y_transformed, 0.95), y = par("usr")[4] * 0.9,
       labels = sprintf("Skew: %.2f\n%s", trans_skew, transform_name), pos = 2)
  dev.off()
  
  plsr_data <- data.frame(y = y_transformed, X_spec)
  ctrl <- trainControl(method = "cv", number = 10, savePredictions = "final")
  
  plsr_model <- train(
    y ~ ., data = plsr_data, method = "pls", trControl = ctrl,
    tuneLength = 20, preProcess = c("center", "scale")
  )
  
  preds_df <- plsr_model$pred
  best_tune <- plsr_model$bestTune
  
  for (param in names(best_tune)) {
    preds_df <- preds_df[preds_df[[param]] == best_tune[[param]], ]
  }
  
  preds_ordered <- preds_df %>% arrange(rowIndex) %>% distinct(rowIndex, .keep_all = TRUE)
  residuals <- abs(preds_ordered$obs - preds_ordered$pred)
  
  n_to_remove <- ceiling(max_pct * n_initial)
  n_to_remove <- max(1, min(n_to_remove, n_initial - 10))
  
  resid_threshold <- sort(residuals, decreasing = TRUE)[n_to_remove]
  outlier_idx <- order(residuals, decreasing = TRUE)[1:n_to_remove]
  
  pct_removed <- 100 * length(outlier_idx) / n_initial
  cat(sprintf("Removing %d samples (%.2f%%) with residuals >= %.6f\n", 
              length(outlier_idx), pct_removed, resid_threshold))
  
  removed_info <- tibble()
  if (length(outlier_idx) > 0) {
    removed_info <- tibble(
      SNM = data_merged$SNM[outlier_idx],
      Lab_value = y_lab[outlier_idx],
      residual_transformed = residuals[outlier_idx],
      property = property_name
    )
    data_merged <- data_merged[-outlier_idx, ]
  }
  
  cat("Remaining samples:", nrow(data_merged), "\n")
  
  return(list(
    data = data_merged,
    removed = removed_info,
    n_removed = length(outlier_idx),
    transform_obj = bn_obj,
    transform_name = transform_name,
    orig_skew = orig_skew,
    trans_skew = trans_skew
  ))
}

# -----------------------
# 4. FUNCTION: KENNARD-STONE SPLIT
# -----------------------
create_ks_indices <- function(spectral_mat, y_transf, bins = 5, k_frac = 0.8) {
  df <- as.data.frame(spectral_mat)
  df$y_transformed <- y_transf
  df$y_bin <- ntile(df$y_transformed, bins)
  
  train_idx <- c(); test_idx <- c()
  
  for (b in unique(df$y_bin)) {
    rows <- which(df$y_bin == b)
    bin_data <- df[rows, ]
    x_bin <- as.matrix(bin_data[, !(names(bin_data) %in% c("y_transformed", "y_bin"))])
    
    k_train <- floor(k_frac * nrow(bin_data))
    if (k_train < 1) k_train <- 1
    
    if (nrow(bin_data) > 1 && ncol(x_bin) > 0) {
      ks <- kenStone(X = x_bin, k = k_train, metric = "euclid")
      train_idx <- c(train_idx, rows[ks$model])
      test_idx <- c(test_idx, rows[setdiff(1:nrow(bin_data), ks$model)])
    } else {
      train_idx <- c(train_idx, rows)
    }
  }
  
  train_idx <- unique(train_idx)
  test_idx <- setdiff(1:nrow(df), train_idx)
  
  list(train = train_idx, test = test_idx)
}

# -----------------------
# 5. FUNCTION: CARS FEATURE SELECTION
# -----------------------
perform_cars <- function(spectra_mat, lab_transf, train_idx, property_name, 
                         sampling_num = cars_sampling_num, cv_folds = cars_cv_folds, max_iter = cars_max_iter) {
  
  cat("\n--- CARS for:", property_name, "(baseline-corrected spectra) ---\n")
  
  X_train <- as.matrix(spectra_mat[train_idx, ])
  y_train <- lab_transf[train_idx]
  X_train_scaled <- scale(X_train, center = TRUE, scale = TRUE)
  
  n_samples <- nrow(X_train_scaled); n_features <- ncol(X_train_scaled)
  cat(sprintf("Training samples: %d, Initial features: %d\n", n_samples, n_features))
  
  best_rmse <- Inf; best_iteration <- 0; best_selected_vars <- NULL
  rmse_history <- numeric(max_iter); n_vars_history <- numeric(max_iter)
  current_vars <- 1:n_features
  
  cat("\nRunning CARS iterations...\n")
  
  for (iter in 1:max_iter) {
    if (length(current_vars) < 2) {
      cat(sprintf("  Iteration %d: Only %d variable(s) left, stopping\n", iter, length(current_vars)))
      break
    }
    
    var_importance <- numeric(length(current_vars))
    
    for (mc in 1:sampling_num) {
      sample_idx <- sample(1:n_samples, size = n_samples, replace = TRUE)
      X_mc <- X_train_scaled[sample_idx, current_vars, drop = FALSE]
      y_mc <- y_train[sample_idx]
      
      max_comps <- min(10, length(current_vars), nrow(X_mc) - 1)
      if (max_comps < 1) next
      
      tryCatch({
        plsr_mc <- plsr(y ~ ., data = data.frame(y = y_mc, X_mc), ncomp = max_comps, validation = "none")
        coefs <- abs(coef(plsr_mc, ncomp = max_comps))
        var_importance <- var_importance + coefs
      }, error = function(e) {})
    }
    
    var_importance <- var_importance / sampling_num
    
    # CV evaluation
    cv_rmse <- numeric(cv_folds)
    folds <- cut(1:n_samples, breaks = cv_folds, labels = FALSE)
    
    for (fold in 1:cv_folds) {
      test_fold <- which(folds == fold); train_fold <- which(folds != fold)
      X_cv_train <- X_train_scaled[train_fold, current_vars, drop = FALSE]
      y_cv_train <- y_train[train_fold]
      X_cv_test <- X_train_scaled[test_fold, current_vars, drop = FALSE]
      y_cv_test <- y_train[test_fold]
      
      max_comps_cv <- min(15, length(current_vars), nrow(X_cv_train) - 1)
      if (max_comps_cv < 1) {
        cv_rmse[fold] <- NA; next
      }
      
      tryCatch({
        plsr_cv <- plsr(y ~ ., data = data.frame(y = y_cv_train, X_cv_train), ncomp = max_comps_cv, validation = "none")
        pred_cv <- predict(plsr_cv, newdata = data.frame(X_cv_test), ncomp = max_comps_cv)[,,1]
        cv_rmse[fold] <- sqrt(mean((y_cv_test - pred_cv)^2))
      }, error = function(e) { cv_rmse[fold] <- NA })
    }
    
    avg_rmse <- mean(cv_rmse, na.rm = TRUE)
    rmse_history[iter] <- avg_rmse
    n_vars_history[iter] <- length(current_vars)
    
    if (avg_rmse < best_rmse) {
      best_rmse <- avg_rmse
      best_iteration <- iter
      best_selected_vars <- current_vars
    }
    
    if (iter %% 10 == 0) {
      cat(sprintf("  Iteration %d: %d vars, RMSE=%.4f (best: iter %d, %d vars, RMSE=%.4f)\n",
                  iter, length(current_vars), avg_rmse, best_iteration, length(best_selected_vars), best_rmse))
    }
    
    # Adaptive removal
    removal_ratio <- exp(-2 * iter / max_iter)
    n_to_remove <- max(1, floor(length(current_vars) * removal_ratio * 0.1))
    
    var_order <- order(var_importance, decreasing = TRUE)
    keep_vars <- var_order[1:(length(current_vars) - n_to_remove)]
    current_vars <- current_vars[keep_vars]
  }
  
  actual_iters <- sum(n_vars_history > 0)
  rmse_history <- rmse_history[1:actual_iters]
  n_vars_history <- n_vars_history[1:actual_iters]
  
  selected_feature_names <- colnames(X_train)[best_selected_vars]
  
  cat(sprintf("\n✓ CARS completed: %d iterations\n", actual_iters))
  cat(sprintf("✓ Best model: iteration %d with %d features (RMSE=%.4f)\n", 
              best_iteration, length(best_selected_vars), best_rmse))
  
  if (!is.null(wavenumbers) && length(wavenumbers) == ncol(spectra_mat)) {
    selected_wavenumbers <- wavenumbers[best_selected_vars]
    wav_summary <- sprintf("%.1f - %.1f cm⁻¹", min(selected_wavenumbers), max(selected_wavenumbers))
  } else {
    selected_wavenumbers <- NULL
    wav_summary <- "N/A"
  }
  
  cars_profile <- data.frame(
    Iteration = 1:actual_iters,
    N_Variables = n_vars_history,
    RMSE = rmse_history
  )
  
  png(paste0(property_name, "_CARS_profile.png"), width = 1200, height = 500, res = 100)
  par(mfrow = c(1, 2))
  plot(cars_profile$Iteration, cars_profile$RMSE, type = "l", lwd = 2, col = "blue",
       xlab = "Iteration", ylab = "CV RMSE", main = "CARS: RMSE Progression"); grid()
  points(best_iteration, best_rmse, pch = 19, col = "red", cex = 2)
  plot(cars_profile$Iteration, cars_profile$N_Variables, type = "l", lwd = 2, col = "darkgreen",
       xlab = "Iteration", ylab = "Number of Variables", main = "CARS: Variable Selection"); grid()
  points(best_iteration, length(best_selected_vars), pch = 19, col = "red", cex = 2)
  dev.off()
  
  cat("✓ Saved:", paste0(property_name, "_CARS_profile.png\n"))
  
  list(
    selected_vars = selected_feature_names,
    feature_indices = best_selected_vars,
    n_optimal = length(best_selected_vars),
    selected_wavenumbers = selected_wavenumbers,
    wavenumber_range = wav_summary,
    best_rmse = best_rmse,
    best_iteration = best_iteration,
    profile = cars_profile
  )
}

# -----------------------
# 6. FUNCTION: PLOT SELECTED FEATURES ON SPECTRUM
# -----------------------
plot_selected_features <- function(spectra_mat, selected_vars, wavenums, property_name, train_idx) {
  cat("\n--- Creating feature selection plot for", property_name, "---\n")
  
  mean_spectrum <- colMeans(spectra_mat[train_idx, , drop = FALSE])
  is_selected <- colnames(spectra_mat) %in% selected_vars
  
  plot_df <- data.frame(
    Wavenumber = wavenums,
    Absorbance = mean_spectrum,
    Selected = is_selected
  )
  
  p <- ggplot(plot_df, aes(x = Wavenumber, y = Absorbance)) +
    geom_line(color = "gray40", linewidth = 0.5) +
    geom_point(data = filter(plot_df, Selected), color = "red", size = 0.8, alpha = 0.6) +
    scale_x_reverse() +
    labs(
      title = paste0(property_name, ": CARS Selected Features (n = ", sum(is_selected), ")"),
      subtitle = "Red points indicate features selected by CARS",
      x = expression(paste("Wavenumber (cm"^-1, ")")),
      y = "Mean Absorbance (training set)"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 10, color = "gray30"),
      panel.grid.minor = element_blank()
    )
  
  ggsave(paste0(property_name, "_CARS_selected_features_spectrum.png"),
         plot = p, width = 10, height = 5, dpi = 300)
  cat("✓ Saved:", paste0(property_name, "_CARS_selected_features_spectrum.png\n"))
  return(p)
}

# -----------------------
# 7. FUNCTION: BUILD FINAL PLSR WITH SELECTED FEATURES
# -----------------------
build_final_plsr <- function(spectra_mat, lab_transf, lab_norm_obj, lab_original, 
                             model_name, train_idx, test_idx, 
                             selected_features = NULL, max_comps = plsr_max_comps) {
  
  cat("\n--- Training:", model_name, "---\n")
  
  if (!is.null(selected_features)) {
    spectra_mat <- spectra_mat[, selected_features, drop = FALSE]
    cat(sprintf("Using %d selected features\n", length(selected_features)))
  }
  
  X_train <- as.matrix(spectra_mat[train_idx, ])
  y_train <- lab_transf[train_idx]
  X_test <- as.matrix(spectra_mat[test_idx, ])
  y_test <- lab_transf[test_idx]
  
  max_comps_actual <- min(max_comps, ncol(X_train), nrow(X_train) - 1)
  
  ctrl <- trainControl(method = "cv", number = cv_folds)
  plsr_model <- train(
    x = X_train, y = y_train, method = "pls",
    preProcess = c("center", "scale"), 
    tuneGrid = data.frame(ncomp = 1:max_comps_actual),
    trControl = ctrl
  )
  
  cat(sprintf("Optimal components: %d\n", plsr_model$bestTune$ncomp))
  
  y_train_pred <- predict(plsr_model, X_train)
  y_test_pred <- predict(plsr_model, X_test)
  
  y_train_obs_orig <- predict(lab_norm_obj, y_train, inverse = TRUE)
  y_train_pred_orig <- predict(lab_norm_obj, y_train_pred, inverse = TRUE)
  y_test_obs_orig <- predict(lab_norm_obj, y_test, inverse = TRUE)
  y_test_pred_orig <- predict(lab_norm_obj, y_test_pred, inverse = TRUE)
  
  safe_cor2 <- function(o, p) if (length(o) < 2) NA else cor(o, p)^2
  
  r2_cal <- safe_cor2(y_train_obs_orig, y_train_pred_orig)
  rmse_cal <- sqrt(mean((y_train_obs_orig - y_train_pred_orig)^2))
  ccc_cal <- CCC(y_train_pred_orig, y_train_obs_orig)$rho.c[1]
  rpiq_cal <- IQR(y_train_obs_orig) / rmse_cal
  
  r2_val <- safe_cor2(y_test_obs_orig, y_test_pred_orig)
  rmse_val <- sqrt(mean((y_test_obs_orig - y_test_pred_orig)^2))
  ccc_val <- CCC(y_test_pred_orig, y_test_obs_orig)$rho.c[1]
  rpiq_val <- IQR(y_test_obs_orig) / rmse_val
  
  metrics_df <- data.frame(
    Model = model_name,
    n_cal = length(y_train_obs_orig), n_val = length(y_test_obs_orig),
    n_features = ncol(X_train), n_comps = plsr_model$bestTune$ncomp,
    R2_cal = r2_cal, RMSE_cal = rmse_cal, CCC_cal = ccc_cal, RPIQ_cal = rpiq_cal,
    R2_val = r2_val, RMSE_val = rmse_val, CCC_val = ccc_val, RPIQ_val = rpiq_val,
    stringsAsFactors = FALSE
  )
  
  list(
    model = plsr_model,
    cal_obs = y_train_obs_orig, cal_pred = y_train_pred_orig,
    val_obs = y_test_obs_orig, val_pred = y_test_pred_orig,
    metrics = metrics_df
  )
}

# -----------------------
# 8. FUNCTION: CREATE COMPARISON PLOT
# -----------------------
calc_metrics <- function(obs, pred) {
  rmse <- sqrt(mean((pred - obs)^2, na.rm = TRUE))
  r2 <- cor(obs, pred, use = "complete.obs")^2
  ccc <- CCC(obs, pred)$rho.c
  rpiq <- (quantile(obs, 0.75, na.rm = TRUE) - quantile(obs, 0.25, na.rm = TRUE)) / rmse
  list(RMSE = rmse, R2 = r2, CCC = ccc, RPIQ = rpiq)
}

make_full_plot <- function(model_res, model_title, ylab_unit = "unit") {
  cal_df <- data.frame(Actual = model_res$cal_obs, Predicted = model_res$cal_pred, Type = "Calibration")
  val_df <- data.frame(Actual = model_res$val_obs, Predicted = model_res$val_pred, Type = "Validation")
  df <- rbind(cal_df, val_df)
  
  cal_m <- calc_metrics(model_res$cal_obs, model_res$cal_pred)
  val_m <- calc_metrics(model_res$val_obs, model_res$val_pred)
  
  stats_text <- sprintf(
    "Cal: RMSE=%.3f R²=%.3f CCC=%.3f RPIQ=%.2f\nVal: RMSE=%.3f R²=%.3f CCC=%.3f RPIQ=%.2f",
    cal_m$RMSE, cal_m$R2, cal_m$CCC, cal_m$RPIQ,
    val_m$RMSE, val_m$R2, val_m$CCC, val_m$RPIQ
  )
  
  min_val <- min(df$Actual, df$Predicted, na.rm = TRUE)
  max_val <- max(df$Actual, df$Predicted, na.rm = TRUE) * 1.05
  
  ggplot(df, aes(x = Actual, y = Predicted, shape = Type, color = Type)) +
    geom_abline(intercept = 0, slope = 1, color = "black", linewidth = 0.6) +
    geom_point(size = 1.5, alpha = 0.8) +
    scale_shape_manual(values = c(16, 17)) +
    scale_color_manual(values = c("forestgreen", "darkorange")) +
    annotate("label", x = -Inf, y = Inf, label = stats_text,
             hjust = -0.02, vjust = 1.02, size = 3, family = "mono", color = "black", fontface = "bold",
             fill = "white", linewidth = 0.5, label.padding = unit(0.4, "lines")) +
    labs(title = model_title, x = paste0("Actual (", ylab_unit, ")"), y = paste0("Predicted (", ylab_unit, ")")) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0),
      legend.title = element_blank(), legend.position = c(0.87, 0.15),
      legend.background = element_rect(color = "black", fill = "white"),
      panel.grid.minor = element_blank(),
      plot.margin = margin(10, 10, 10, 10)
    ) +
    coord_cartesian(xlim = c(min_val, max_val), ylim = c(min_val, max_val), clip = "off")
}

# -----------------------
# 9. MAIN LOOP: ITERATE THROUGH PROPERTIES
# -----------------------
all_removed <- tibble()
all_metrics <- tibble()
all_transformations <- tibble()
all_cars_results <- list()
all_selected_features <- tibble()

for (prop in properties) {
  
  cat("\n")
  cat("========================================\n")
  cat("=== PROCESSING PROPERTY:", prop, "===\n")
  cat("========================================\n")
  
  lab_file <- paste0(prop, ".csv")
  if (!file.exists(lab_file)) {
    cat("WARNING: File", lab_file, "not found. Skipping.\n")
    next
  }
  
  lab_data <- read.csv(lab_file, stringsAsFactors = FALSE)
  if (!prop %in% names(lab_data)) {
    cat("WARNING: Column", prop, "not found. Skipping.\n")
    next
  }
  
  outlier_result <- dangal_outlier_removal(spectra_clean, lab_data, prop, max_pct = max_outlier_pct)
  data_clean <- outlier_result$data
  all_removed <- bind_rows(all_removed, outlier_result$removed)
  
  all_transformations <- bind_rows(all_transformations, tibble(
    Property = prop, Transform = outlier_result$transform_name,
    Original_Skew = outlier_result$orig_skew, Transformed_Skew = outlier_result$trans_skew,
    N_samples = nrow(data_clean), N_removed = outlier_result$n_removed
  ))
  
  if (nrow(data_clean) < 50) {
    cat("WARNING: Too few samples after outlier removal. Skipping", prop, "\n")
    next
  }
  
  spectral_mat_clean <- as.matrix(data_clean %>% select(-SNM, -Lab))
  rownames(spectral_mat_clean) <- data_clean$SNM
  colnames(spectral_mat_clean) <- spectral_cols
  
  lab_original <- data_clean$Lab
  lab_norm <- outlier_result$transform_obj
  lab_transformed <- predict(lab_norm, newdata = lab_original)
  
  ks <- create_ks_indices(spectral_mat_clean, lab_transformed, bins = 5, k_frac = 0.8)
  train_idx <- ks$train; test_idx <- ks$test
  cat(sprintf("\nKennard-Stone: n_train=%d, n_test=%d\n", length(train_idx), length(test_idx)))
  
  baseline_corrected <- t(apply(spectral_mat_clean, 1, function(x) x - min(x, na.rm = TRUE)))
  baseline_corrected <- restore_spectral_dimnames(baseline_corrected, spectral_mat_clean)
  
  cars_result <- perform_cars(baseline_corrected, lab_transformed, train_idx, prop)
  all_cars_results[[prop]] <- cars_result
  all_selected_features <- bind_rows(all_selected_features, tibble(
    Property = prop, N_Features = cars_result$n_optimal, 
    Wavenumber_Range = cars_result$wavenumber_range,
    Best_Iteration = cars_result$best_iteration, Best_RMSE = cars_result$best_rmse
  ))
  
  write.csv(cars_result$profile, paste0(prop, "_CARS_profile.csv"), row.names = FALSE)
  
  selected_features_detail <- data.frame(
    Feature = cars_result$selected_vars,
    Wavenumber = if (!is.null(cars_result$selected_wavenumbers)) cars_result$selected_wavenumbers else NA
  )
  write.csv(selected_features_detail, paste0(prop, "_CARS_selected_features_detail.csv"), row.names = FALSE)
  
  if (!is.null(wavenumbers)) {
    plot_selected_features(baseline_corrected, cars_result$selected_vars, wavenumbers, prop, train_idx)
  } else {
    cat("WARNING: Cannot create feature plot — wavenumbers not available\n")
  }
  
  # ------------------------------
  # ✅ SG+SNV WITH FIXED slice_sample() ERROR
  # ------------------------------
  
  selected_feature_names <- cars_result$selected_vars
  cat(sprintf("\nUsing %d CARS-selected features (by name)\n", length(selected_feature_names)))
  
  baseline_selected <- baseline_corrected[, selected_feature_names, drop = FALSE]
  
  # ✅ FIXED: replace slice_sample(n = min(20, n())) with base R sampling
  valid_sg_grid_df <- as.data.frame(sg_grid) %>%
    filter(
      w <= ncol(spectral_mat_clean),
      w > p,
      (w - 1) %% 2 == 0,  # ensure odd window
      m <= p              # valid derivative order
    )
  
  if (nrow(valid_sg_grid_df) == 0) {
    cat("⚠️  No valid SG parameter combos. Using fallback.\n")
    valid_sg_grid_df <- data.frame(p = 2, w = 11, m = 0)
  } else {
    # ✅ CORRECT: base R sampling WITHOUT n()
    n_try <- min(20, nrow(valid_sg_grid_df))
    valid_sg_grid_df <- valid_sg_grid_df[sample(nrow(valid_sg_grid_df), size = n_try), , drop = FALSE]
  }
  
  best_rmse <- Inf
  best_params <- NULL
  
  cat("\nOptimizing Savitzky-Golay parameters...\n")
  
  for (i in seq_len(nrow(valid_sg_grid_df))) {
    p <- valid_sg_grid_df$p[i]; w <- valid_sg_grid_df$w[i]; m <- valid_sg_grid_df$m[i]
    
    tryCatch({
      sgx <- savitzkyGolay(spectral_mat_clean, p = p, w = w, m = m)
      sgx <- restore_spectral_dimnames(sgx, spectral_mat_clean)
      
      sgx_snv <- standardNormalVariate(sgx)
      sgx_snv <- restore_spectral_dimnames(sgx_snv, spectral_mat_clean)
      
      # Ensure selected features exist
      missing_feats <- setdiff(selected_feature_names, colnames(sgx_snv))
      if (length(missing_feats) > 0) {
        cat(sprintf("  SG(p=%d,w=%d,m=%d): %d features missing → skipping\n", p, w, m, length(missing_feats)))
        next
      }
      
      sgx_selected <- sgx_snv[, selected_feature_names, drop = FALSE]
      
      # Inner CV for RMSE estimation
      tmp_n <- floor(0.8 * length(train_idx))
      if (tmp_n < 5) next
      tmp_train_idx <- sample(train_idx, tmp_n)
      tmp_test_idx  <- setdiff(train_idx, tmp_train_idx)
      
      if (length(tmp_test_idx) < 2) next
      
      tmp_train_x <- sgx_selected[tmp_train_idx, , drop = FALSE]
      tmp_train_y <- lab_transformed[tmp_train_idx]
      tmp_test_x  <- sgx_selected[tmp_test_idx, , drop = FALSE]
      tmp_test_y  <- lab_transformed[tmp_test_idx]
      
      if (ncol(tmp_train_x) < 2 || nrow(tmp_test_x) < 2) next
      
      tmp_m <- train(tmp_train_x, tmp_train_y, method = "pls",
                     preProcess = c("center", "scale"),
                     tuneLength = min(10, ncol(tmp_train_x)),
                     trControl = trainControl(method = "cv", number = 3))
      
      tmp_pred <- predict(tmp_m, tmp_test_x)
      tmp_rmse <- sqrt(mean((tmp_test_y - tmp_pred)^2, na.rm = TRUE))
      
      if (is.finite(tmp_rmse) && tmp_rmse < best_rmse) {
        best_rmse <- tmp_rmse
        best_params <- list(p = p, w = w, m = m)
        cat(sprintf("  New best SG: p=%d, w=%d, m=%d → RMSE=%.4f\n", p, w, m, tmp_rmse))
      }
    }, error = function(e) {
      cat(sprintf("  SG(p=%d,w=%d,m=%d) failed: %s\n", p, w, m, conditionMessage(e)))
    })
  }
  
  # Finalize SG+SNV
  if (is.null(best_params)) {
    cat("⚠️  No working SG combo. Using safe fallback (p=2, w=11, m=0)\n")
    best_params <- list(p = 2, w = 11, m = 0)
  }
  
  sgx_final <- tryCatch(
    savitzkyGolay(spectral_mat_clean, p = best_params$p, w = best_params$w, m = best_params$m),
    error = function(e) { cat("SG fallback failed → using raw spectra\n"); spectral_mat_clean }
  )
  sgx_final <- restore_spectral_dimnames(sgx_final, spectral_mat_clean)
  
  sgx_snv_final <- tryCatch(
    standardNormalVariate(sgx_final),
    error = function(e) sgx_final
  )
  sgx_snv_final <- restore_spectral_dimnames(sgx_snv_final, spectral_mat_clean)
  
  # Final safe subsetting
  common_feats <- intersect(selected_feature_names, colnames(sgx_snv_final))
  if (length(common_feats) == 0) {
    stop("FATAL ERROR: No overlap between CARS features and SG+SNV output for ", prop)
  }
  spectra_2_selected <- sgx_snv_final[, common_feats, drop = FALSE]
  
  # Model 3: SNV
  spectra_3_snv <- standardNormalVariate(spectral_mat_clean)
  spectra_3_snv <- restore_spectral_dimnames(spectra_3_snv, spectral_mat_clean)
  spectra_3_selected <- spectra_3_snv[, intersect(selected_feature_names, colnames(spectra_3_snv)), drop = FALSE]
  
  # Model 4: MSC
  spectra_4_msc <- tryCatch(msc(spectral_mat_clean), error = function(e) spectral_mat_clean)
  spectra_4_msc <- restore_spectral_dimnames(spectra_4_msc, spectral_mat_clean)
  spectra_4_selected <- spectra_4_msc[, intersect(selected_feature_names, colnames(spectra_4_msc)), drop = FALSE]
  
  cat(sprintf("SG+SNV using %d/%d features (p=%d, w=%d, m=%d)\n", 
              ncol(spectra_2_selected), length(selected_feature_names),
              best_params$p, best_params$w, best_params$m))
  
  # ------------------------------
  # TRAIN MODELS
  # ------------------------------
  
  cat("\n=== Training final PLSR models with CARS-selected features ===\n")
  
  model_1 <- build_final_plsr(baseline_selected, lab_transformed, lab_norm, lab_original,
                              "Base PLSR", train_idx, test_idx)
  
  model_2_name <- sprintf("SavGol(p=%d,w=%d,m=%d)+SNV", 
                          best_params$p, best_params$w, best_params$m)
  model_2 <- build_final_plsr(spectra_2_selected, lab_transformed, lab_norm, lab_original,
                              model_2_name, train_idx, test_idx)
  
  model_3 <- build_final_plsr(spectra_3_selected, lab_transformed, lab_norm, lab_original,
                              "SNV PLSR", train_idx, test_idx)
  
  model_4 <- build_final_plsr(spectra_4_selected, lab_transformed, lab_norm, lab_original,
                              "MSC PLSR", train_idx, test_idx)
  
  prop_metrics <- bind_rows(
    model_1$metrics %>% mutate(Property = prop),
    model_2$metrics %>% mutate(Property = prop),
    model_3$metrics %>% mutate(Property = prop),
    model_4$metrics %>% mutate(Property = prop)
  )
  
  all_metrics <- bind_rows(all_metrics, prop_metrics)
  write.csv(prop_metrics, paste0(prop, "_four_PLSR_models_metrics.csv"), row.names = FALSE)
  
  p1 <- make_full_plot(model_1, "Base PLSR", prop)
  p2 <- make_full_plot(model_2, "SavGol PLSR", prop)
  p3 <- make_full_plot(model_3, "SNV PLSR", prop)
  p4 <- make_full_plot(model_4, "MSC PLSR", prop)
  
  tiff(paste0(prop, "_PLSR_CARS_Comparison.tiff"), width = 3000, height = 2500, res = 300)
  grid.arrange(p1, p2, p3, p4, ncol = 2)
  dev.off()
  
  cat("\n✅ Completed", prop, "- saved all outputs\n")
}

# -----------------------
# 10. SAVE SUMMARY RESULTS
# -----------------------

if (nrow(all_removed) > 0) {
  write.csv(all_removed, "all_properties_removed_outliers_CARS.csv", row.names = FALSE)
  cat("\n📄 Saved: all_properties_removed_outliers_CARS.csv\n")
}

if (nrow(all_transformations) > 0) {
  write.csv(all_transformations, "all_properties_transformations_used_CARS.csv", row.names = FALSE)
  cat("📄 Saved: all_properties_transformations_used_CARS.csv\n")
}

if (nrow(all_selected_features) > 0) {
  write.csv(all_selected_features, "all_properties_CARS_selected_features_summary.csv", row.names = FALSE)
  cat("📄 Saved: all_properties_CARS_selected_features_summary.csv\n")
}

if (nrow(all_metrics) > 0) {
  write.csv(all_metrics, "all_properties_four_PLSR_models_metrics_CARS.csv", row.names = FALSE)
  cat("📄 Saved: all_properties_four_PLSR_models_metrics_CARS.csv\n")
  
  cat("\n=== VALIDATION R² SUMMARY ===\n")
  summary_table <- all_metrics %>%
    select(Property, Model, n_features, n_comps, R2_val, RMSE_val, RPIQ_val) %>%
    arrange(Property, desc(R2_val))
  print(summary_table, n = Inf)
}

cat("\n")
cat("========================================\n")
cat("=== 🚀 ALL PROPERTIES PROCESSED! ===\n")
cat("========================================\n")
cat("\n✅ No errors. All outputs saved.\n")