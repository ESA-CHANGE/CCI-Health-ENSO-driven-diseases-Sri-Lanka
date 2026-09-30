
# ==============================================================================
# SECOND-STAGE META-REGRESSION OF ONI–LEPTOSPIROSIS ASSOCIATIONS
# ==============================================================================
#
# PURPOSE
# -------
# The first-stage DLNM models estimate district-specific, potentially
# non-linear associations between the Oceanic Niño Index (ONI) and
# leptospirosis incidence.
#
# The purpose of this second-stage analysis is to determine whether
# between-district differences in the estimated ONI–leptospirosis association
# are systematically related to district-level environmental and landscape
# characteristics.
#
# The analysis therefore treats the district-level environmental variables as
# META-REGRESSORS of the first-stage ONI association.
#
#
# ---------------------------------------------------------------------------
# FIRST STAGE
# ---------------------------------------------------------------------------
#
# For each district d, the first-stage DLNM provides:
#
#   theta_d = vector of ONI cross-basis coefficients
#
# together with:
#
#   S_d = variance-covariance matrix of theta_d
#
# In the current analysis, the ONI cross-basis has three coefficients.
#
#
# ---------------------------------------------------------------------------
# SECOND STAGE
# ---------------------------------------------------------------------------
#
# A multivariate random-effects meta-analysis is used to pool the
# district-specific coefficient vectors:
#
#   theta_d = beta_0 + u_d
#
# This is the NULL MODEL.
#
# To investigate whether district characteristics explain heterogeneity,
# each candidate district-level variable X is introduced as a meta-regressor:
#
#   theta_d = beta_0 + beta_1 X_d + u_d
#
# Because theta_d is a vector of three first-stage coefficients, beta_1 is
# also a vector of three coefficients.
#
# The joint Wald test therefore tests:
#
#   H0: beta_1,1 = beta_1,2 = beta_1,3 = 0
#
# i.e. whether the complete ONI association profile varies systematically
# with the district-level characteristic.
#
#
# ---------------------------------------------------------------------------
# MODEL COMPARISON
# ---------------------------------------------------------------------------
#
# For each meta-regressor, the univariable model is compared with an
# intercept-only NULL model.
#
# The comparison includes:
#
#   1. AIC
#   2. BIC
#   3. Delta AIC
#   4. Delta BIC
#   5. Likelihood-ratio statistic
#   6. Likelihood-ratio degrees of freedom
#   7. Likelihood-ratio P value
#
# IMPORTANT:
#
# The null model is fitted separately for each covariate using exactly the
# same districts included in that covariate's meta-regression.
#
# This is essential because some district-level covariates contain missing
# values. Comparing a 23-district model with a 25-district null model would
# confound covariate effects with differences in analytical sample.
#
#
# ---------------------------------------------------------------------------
# ESTIMATION METHOD
# ---------------------------------------------------------------------------
#
# REML is retained for the primary second-stage meta-regression, consistent
# with the original analysis.
#
# However, likelihood-based comparison of models with different fixed-effects
# structures should be performed using ML rather than REML.
#
# Therefore:
#
#   REML -> primary parameter estimation, Q, I2 and Wald statistics
#
#   ML   -> null-versus-meta-regression model comparison
#
# After selecting the fixed-effects structure, the final model can be
# estimated with REML.
#
#
# ---------------------------------------------------------------------------
# TARGET DISEASE
# ---------------------------------------------------------------------------
#
# Leptospirosis
#
# =============================================================================


# ==============================================================================
# 5.1 LOAD REQUIRED PACKAGES
# ==============================================================================

library(dplyr)
library(mvmeta)
library(dlnm)


# ==============================================================================
# 5.2 ENVIRONMENT INITIALIZATION AND COVARIATE MATRIX
# ==============================================================================

cat("\n")
cat("============================================================\n")
cat("SECTION 5: SECOND-STAGE ONI–LEPTOSPIROSIS META-REGRESSION\n")
cat("============================================================\n")

cat("\n[1/7] Loading district-level covariate matrix...\n")


# ------------------------------------------------------------------------------
# Covariate matrix
# ------------------------------------------------------------------------------

covariate_path <- "~/Documents/Heidelberg/ESA_climate/analysis/covariate_matrix_Full.csv"

if (!file.exists(covariate_path)) {
  stop(
    paste(
      "ERROR: Covariate matrix file not found at:",
      covariate_path
    )
  )
}

covariate_matrix <- read.csv(
  covariate_path,
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------------------------
# Identify LCCS land-cover variables
# ------------------------------------------------------------------------------

lccs_cols <- grep(
  "^lccs_",
  names(covariate_matrix),
  value = TRUE
)

cat(
  "Found",
  length(lccs_cols),
  "LCCS land-cover variables.\n"
)


# ==============================================================================
# 5.3 DEFINE ADDITIONAL DISTRICT-LEVEL META-REGRESSORS
# ==============================================================================
#
# In addition to land-cover composition, the analysis can evaluate:
#
#   - Human settlement intensity
#   - Rainfall
#   - NDVI
#   - Soil moisture
#
# Long-term summaries and changes are used rather than all annual values
# simultaneously. This keeps the second-stage model parsimonious given the
# relatively small number of districts.
# ==============================================================================


# ------------------------------------------------------------------------------
# Human settlement variables
# ------------------------------------------------------------------------------

hsi_cols <- intersect(
  c(
    "hsi_2016",
    "hsi_2025",
    "hsi_change",
    "hsi_relative_change"
  ),
  names(covariate_matrix)
)


# ------------------------------------------------------------------------------
# Rainfall variables
# ------------------------------------------------------------------------------

rain_cols <- intersect(
  c(
    "rain_monthly_total_mean_2007_2025",
    "rain_change_2007_2025"
  ),
  names(covariate_matrix)
)


# ------------------------------------------------------------------------------
# NDVI variables
# ------------------------------------------------------------------------------

ndvi_cols <- intersect(
  c(
    "ndvi_mean_2007_2025",
    "ndvi_change_2007_2025"
  ),
  names(covariate_matrix)
)


# ------------------------------------------------------------------------------
# Soil-moisture variables
# ------------------------------------------------------------------------------

soil_cols <- intersect(
  c(
    "soil_moisture_mean_2007_2025",
    "soil_moisture_change_2007_2025",
    "soil_moisture_relative_change_2007_2025"
  ),
  names(covariate_matrix)
)


# ------------------------------------------------------------------------------
# Complete set of candidate meta-regressors
# ------------------------------------------------------------------------------

meta_regressors <- unique(
  c(
    lccs_cols,
    hsi_cols,
    rain_cols,
    ndvi_cols,
    soil_cols
  )
)


cat(
  "Total candidate meta-regressors:",
  length(meta_regressors),
  "\n"
)


# ==============================================================================
# 5.4 VERIFY FIRST-STAGE ONI OBJECTS
# ==============================================================================
#
# The second-stage analysis assumes that the following objects have already
# been generated by the first-stage DLNM analysis:
#
#   districts
#   keep_oni
#   coef_oni_all
#   vcov_oni_all
#   b_pred_oni
#   pred_oni_ov
#
# These objects contain the district-specific ONI association estimates and
# their covariance matrices.
# ==============================================================================


required_objects <- c(
  "districts",
  "keep_oni",
  "coef_oni_all",
  "vcov_oni_all",
  "b_pred_oni",
  "pred_oni_ov"
)

missing_objects <- required_objects[
  !sapply(
    required_objects,
    exists,
    inherits = TRUE
  )
]

if (length(missing_objects) > 0) {
  stop(
    paste(
      "ERROR: The following first-stage objects are missing:",
      paste(missing_objects, collapse = ", ")
    )
  )
}


# ------------------------------------------------------------------------------
# Determine the number of first-stage ONI basis coefficients automatically.
# ------------------------------------------------------------------------------

n_basis_coef <- ncol(
  coef_oni_all
)

cat(
  "Number of first-stage ONI basis coefficients:",
  n_basis_coef,
  "\n"
)


# ==============================================================================
# 5.5 EXPLICIT DISTRICT ALIGNMENT
# ==============================================================================
#
# The order of districts in the first-stage coefficient matrix must correspond
# exactly to the order of district-level covariates.
#
# We therefore create an explicit lookup table rather than relying on row
# ordering.
# ==============================================================================


first_stage_lookup <- data.frame(
  district = districts,
  first_stage_row = seq_along(districts),
  keep_oni = keep_oni,
  stringsAsFactors = FALSE
)


loop_metadata_base <- first_stage_lookup %>%
  left_join(
    covariate_matrix %>%
      select(
        district_name,
        district_pcode,
        all_of(meta_regressors)
      ),
    by = c(
      "district" = "district_name"
    )
  )


# ------------------------------------------------------------------------------
# Check the alignment
# ------------------------------------------------------------------------------

cat("\nDistrict alignment check:\n")

print(
  loop_metadata_base %>%
    select(
      district,
      district_pcode,
      first_stage_row,
      keep_oni
    ),
  row.names = FALSE
)


# ==============================================================================
# 5.6 FUNCTION FOR JOINT WALD TEST
# ==============================================================================
#
# The univariable meta-regression has:
#
#   n_basis_coef intercept coefficients
#
# followed by:
#
#   n_basis_coef meta-regression coefficients.
#
# For the current ONI cross-basis:
#
#   3 intercept coefficients
#   3 covariate coefficients
#
# The joint Wald statistic evaluates whether all covariate coefficients are
# simultaneously zero.
#
# The function is written dynamically so that the code remains valid if the
# first-stage ONI cross-basis is changed in the future.
# ==============================================================================


calculate_joint_wald <- function(
    model,
    n_basis_coef
) {
  
  beta_all <- coef(model)
  
  V_all <- vcov(model)
  
  
  # Identify the block corresponding to the meta-regressor.
  
  slope_idx <- (
    n_basis_coef + 1
  ):(2 * n_basis_coef)
  
  
  beta_slopes <- beta_all[
    slope_idx
  ]
  
  
  V_slopes <- V_all[
    slope_idx,
    slope_idx,
    drop = FALSE
  ]
  
  
  # Check numerical validity.
  
  if (
    any(!is.finite(beta_slopes)) ||
    any(!is.finite(V_slopes))
  ) {
    
    return(
      list(
        statistic = NA_real_,
        df = n_basis_coef,
        p_value = NA_real_
      )
    )
  }
  
  
  # Calculate the Wald statistic:
  #
  #       W = beta' V^-1 beta
  
  V_inverse <- tryCatch(
    solve(V_slopes),
    error = function(e) {
      MASS::ginv(V_slopes)
    }
  )
  
  
  wald_stat <- as.numeric(
    t(beta_slopes) %*%
      V_inverse %*%
      beta_slopes
  )
  
  
  wald_p <- pchisq(
    wald_stat,
    df = n_basis_coef,
    lower.tail = FALSE
  )
  
  
  return(
    list(
      statistic = wald_stat,
      df = n_basis_coef,
      p_value = wald_p
    )
  )
}


# ==============================================================================
# 5.7 FUNCTION FOR GENERATING NON-LINEAR ONI PREDICTION FIGURES
# ==============================================================================
#
# This function handles visualization separately from the statistical model
# fitting.
#
# The figure is generated directly in the active graphics device and is NOT
# saved to disk.
#
# For each meta-regressor, the fitted ONI–leptospirosis association is
# reconstructed at:
#
#   - 10th percentile of the covariate
#   - 50th percentile (median)
#   - 90th percentile
#
# This allows the second-stage meta-regression to be visualized as a change
# in the entire non-linear ONI association profile across the environmental
# gradient.
#
# ONI = 0 is used as the centering value, meaning that the relative risks are
# expressed relative to an ONI value of zero.
# ==============================================================================


plot_modulation_curve <- function(
    fit_model,
    model_data,
    covariate_name,
    basis_object,
    overall_prediction,
    make_plot = TRUE
) {
  
  # --------------------------------------------------------------------------
  # Calculate the 10th, 50th and 90th percentiles of the meta-regressor.
  # --------------------------------------------------------------------------
  
  cov_percentiles <- quantile(
    model_data$current_cov,
    probs = c(
      0.10,
      0.50,
      0.90
    ),
    na.rm = TRUE
  )
  
  
  # --------------------------------------------------------------------------
  # Generate predictions at each covariate percentile.
  # --------------------------------------------------------------------------
  
  cov_predictions <- list()
  
  
  for (i in seq_along(cov_percentiles)) {
    
    prediction <- predict(
      fit_model,
      newdata = data.frame(
        current_cov = cov_percentiles[i]
      ),
      vcov = TRUE
    )
    
    
    cov_predictions[[i]] <- crosspred(
      basis = basis_object,
      coef = as.numeric(
        prediction$fit
      ),
      vcov = prediction$vcov,
      model.link = "log",
      cen = 0
    )
  }
  
  
  names(cov_predictions) <- c(
    "Low_10th",
    "Median_50th",
    "High_90th"
  )
  
  
  # --------------------------------------------------------------------------
  # Optional figure.
  #
  # The figure is displayed but not written to disk.
  # --------------------------------------------------------------------------
  
  if (make_plot) {
    
    plot(
      overall_prediction$predvar,
      overall_prediction$allRRfit,
      type = "n",
      ylim = c(0.3, 5.0),
      log = "y",
      las = 1,
      bty = "l",
      xlab = "Oceanic Niño Index (ONI)",
      ylab = "Leptospirosis Relative Risk",
      main = paste(
        "ONI–Leptospirosis association\nmodulated by",
        covariate_name
      )
    )
    
    
    # National pooled ONI association.
    
    lines(
      overall_prediction$predvar,
      overall_prediction$allRRfit,
      lwd = 3,
      col = "#333333"
    )
    
    
    # Low covariate level.
    
    lines(
      cov_predictions$Low_10th$predvar,
      cov_predictions$Low_10th$allRRfit,
      col = "#2B8CBE",
      lwd = 2.5,
      lty = 2
    )
    
    
    # Median covariate level.
    
    lines(
      cov_predictions$Median_50th$predvar,
      cov_predictions$Median_50th$allRRfit,
      col = "#FEB24C",
      lwd = 2.5,
      lty = 2
    )
    
    
    # High covariate level.
    
    lines(
      cov_predictions$High_90th$predvar,
      cov_predictions$High_90th$allRRfit,
      col = "#E31A1C",
      lwd = 2.5,
      lty = 2
    )
    
    
    # Reference relative risk.
    
    abline(
      h = 1,
      lty = 3,
      col = "grey50"
    )
    
    
    legend(
      "topleft",
      bty = "n",
      lwd = c(
        3,
        2.5,
        2.5,
        2.5
      ),
      lty = c(
        1,
        2,
        2,
        2
      ),
      col = c(
        "#333333",
        "#2B8CBE",
        "#FEB24C",
        "#E31A1C"
      ),
      legend = c(
        "National pooled association",
        "10th percentile",
        "50th percentile",
        "90th percentile"
      ),
      cex = 0.75
    )
  }
  
  
  return(
    list(
      percentiles = cov_percentiles,
      predictions = cov_predictions
    )
  )
}


# ==============================================================================
# 5.8 STORAGE OBJECTS
# ==============================================================================
#
# Three types of information are retained:
#
#   1. REML meta-regression models
#   2. ML models used for null comparison
#   3. Prediction objects for visualization
#
# Figures themselves are NOT saved.
# ==============================================================================


meta_model_registry <- list()

null_model_registry <- list()

ml_model_registry <- list()

prediction_curve_list <- list()


# ==============================================================================
# 5.9 RESULT TABLE
# ==============================================================================
#
# The table retains the original diagnostics:
#
#   AIC
#   BIC
#   LogLik
#   Cochran Q
#   Q P value
#   I-squared
#   Wald statistic
#   Wald P value
#
# and adds:
#
#   Null AIC
#   Null BIC
#   Null LogLik
#   Delta AIC
#   Delta BIC
#   LR statistic
#   LR df
#   LR P value
#
# Numeric P values are retained separately from formatted P values.
# ==============================================================================


model_comparison_table <- data.frame(
  Covariate = character(),
  N_districts = integer(),
  
  # Original REML model diagnostics
  AIC = numeric(),
  BIC = numeric(),
  LogLik = numeric(),
  Cochran_Q = numeric(),
  Q_df = integer(),
  Q_p_value = character(),
  I_squared_pct = numeric(),
  Wald_test_stat = numeric(),
  Wald_df = integer(),
  Wald_p_value = character(),
  
  # Numeric Wald P value for downstream screening
  Wald_p_numeric = numeric(),
  
  # Null-model comparison
  Null_AIC = numeric(),
  Null_BIC = numeric(),
  Null_LogLik = numeric(),
  
  Delta_AIC = numeric(),
  Delta_BIC = numeric(),
  
  LR_stat = numeric(),
  LR_df = integer(),
  LR_p_value = character(),
  LR_p_numeric = numeric(),
  
  stringsAsFactors = FALSE
)


# ==============================================================================
# 5.10 MAIN UNIVARIABLE META-REGRESSION LOOP
# ==============================================================================
#
# Every candidate district-level characteristic is evaluated separately.
#
# For each variable:
#
#   STEP A:
#       Identify districts with valid first-stage ONI estimates.
#
#   STEP B:
#       Remove districts with missing values of the current meta-regressor.
#
#   STEP C:
#       Fit the NULL MODEL using those exact districts.
#
#   STEP D:
#       Fit the univariable meta-regression using REML.
#
#   STEP E:
#       Fit corresponding NULL and univariable models using ML for likelihood
#       comparison.
#
#   STEP F:
#       Extract the original REML diagnostics.
#
#   STEP G:
#       Calculate the joint Wald test.
#
#   STEP H:
#       Calculate Delta AIC, Delta BIC and likelihood-ratio statistics.
#
#   STEP I:
#       Generate ONI prediction curves at the 10th, 50th and 90th percentiles.
#
#   STEP J:
#       Display the figure in the active graphics device.
#
# ==============================================================================


cat("\n")
cat("============================================================\n")
cat("UNIVARIABLE META-REGRESSION ANALYSIS\n")
cat("============================================================\n")


for (cov in meta_regressors) {
  
  
  cat("\n")
  cat("------------------------------------------------------------\n")
  cat("Processing meta-regressor:", cov, "\n")
  cat("------------------------------------------------------------\n")
  
  
  # ============================================================================
  # 5.10.1 CHECK COVARIATE VARIANCE
  # ============================================================================
  
  cov_values <- covariate_matrix[[cov]]
  
  if (
    all(is.na(cov_values)) ||
    var(
      cov_values,
      na.rm = TRUE
    ) == 0
  ) {
    
    cat(
      "Skipping:",
      cov,
      "-> zero or undefined variance.\n"
    )
    
    next
  }
  
  
  # ============================================================================
  # 5.10.2 CREATE COMPLETE-CASE ANALYTICAL DATASET
  # ============================================================================
  #
  # Only districts with:
  #
  #   - a valid first-stage ONI model
  #   - a non-missing value of the current meta-regressor
  #
  # are included.
  #
  # The same district set is subsequently used for BOTH the null model and
  # the univariable model.
  # ============================================================================
  
  
  loop_metadata <- loop_metadata_base %>%
    select(
      district,
      district_pcode,
      first_stage_row,
      keep_oni,
      current_cov = all_of(cov)
    ) %>%
    filter(
      keep_oni,
      !is.na(current_cov)
    )
  
  
  n_districts_model <- nrow(
    loop_metadata
  )
  
  
  if (
    n_districts_model < 5
  ) {
    
    cat(
      "Skipping:",
      cov,
      "-> fewer than 5 complete districts.\n"
    )
    
    next
  }
  
  
  first_stage_rows <- loop_metadata$first_stage_row
  
  
  cat(
    "Number of districts:",
    n_districts_model,
    "\n"
  )
  
  
  # ============================================================================
  # 5.10.3 FIT THE NULL MODEL USING REML
  # ============================================================================
  #
  # This model contains no district-level meta-regressor:
  #
  #   theta_d = beta_0 + u_d
  #
  # It represents the baseline random-effects meta-analysis for the exact
  # analytical sample used by the current covariate.
  # ============================================================================
  
  
  null_model_reml <- mvmeta(
    coef_oni_all[
      first_stage_rows,
      ,
      drop = FALSE
    ] ~ 1,
    
    S = vcov_oni_all[
      first_stage_rows
    ],
    
    data = loop_metadata,
    
    method = "reml"
  )
  
  
  # ============================================================================
  # 5.10.4 FIT THE UNIVARIABLE META-REGRESSION USING REML
  # ============================================================================
  #
  # This is the primary model corresponding to the original pipeline.
  #
  #   theta_d = beta_0 + beta_1 X_d + u_d
  #
  # REML is retained for estimation of the random-effects model.
  # ============================================================================
  
  
  fit_model <- mvmeta(
    coef_oni_all[
      first_stage_rows,
      ,
      drop = FALSE
    ] ~ current_cov,
    
    S = vcov_oni_all[
      first_stage_rows
    ],
    
    data = loop_metadata,
    
    method = "reml"
  )
  
  
  meta_model_registry[[cov]] <- fit_model
  
  null_model_registry[[cov]] <- null_model_reml
  
  
  # ============================================================================
  # 5.10.5 FIT NULL AND META-REGRESSION MODELS USING ML
  # ============================================================================
  #
  # These ML models are used ONLY for comparing different fixed-effects
  # structures.
  #
  # The random-effects variance components are therefore estimated under ML
  # for both models.
  #
  # This allows valid likelihood-based comparison between:
  #
  #       NULL:  ~ 1
  #
  # and
  #
  #       META:  ~ current_cov
  # ============================================================================
  
  
  null_model_ml <- mvmeta(
    coef_oni_all[
      first_stage_rows,
      ,
      drop = FALSE
    ] ~ 1,
    
    S = vcov_oni_all[
      first_stage_rows
    ],
    
    data = loop_metadata,
    
    method = "ml"
  )
  
  
  fit_model_ml <- mvmeta(
    coef_oni_all[
      first_stage_rows,
      ,
      drop = FALSE
    ] ~ current_cov,
    
    S = vcov_oni_all[
      first_stage_rows
    ],
    
    data = loop_metadata,
    
    method = "ml"
  )
  
  
  null_model_registry[[paste0(
    cov,
    "_ML"
  )]] <- null_model_ml
  
  
  ml_model_registry[[cov]] <- fit_model_ml
  
  
  # ============================================================================
  # 5.10.6 EXTRACT ORIGINAL REML MODEL DIAGNOSTICS
  # ============================================================================
  #
  # These are the statistics from the original pipeline and should remain the
  # primary descriptive diagnostics for the meta-regression.
  # ============================================================================
  
  
  model_summary <- summary(
    fit_model
  )
  
  
  aic_val <- AIC(
    fit_model
  )
  
  
  bic_val <- BIC(
    fit_model
  )
  
  
  loglik_val <- as.numeric(
    logLik(
      fit_model
    )
  )
  
  
  # ----------------------------------------------------------------------------
  # Multivariate Cochran Q
  # ----------------------------------------------------------------------------
  
  q_stat <- if (
    !is.null(model_summary$qstat$Q)
  ) {
    
    as.numeric(
      model_summary$qstat$Q
    )[1]
    
  } else {
    
    NA_real_
    
  }
  
  
  q_df <- if (
    !is.null(model_summary$qstat$df)
  ) {
    
    as.numeric(
      model_summary$qstat$df
    )[1]
    
  } else {
    
    NA_integer_
    
  }
  
  
  q_p <- if (
    !is.null(model_summary$qstat$pvalue)
  ) {
    
    as.numeric(
      model_summary$qstat$pvalue
    )[1]
    
  } else {
    
    NA_real_
    
  }
  
  
  # ----------------------------------------------------------------------------
  # I-squared
  # ----------------------------------------------------------------------------
  #
  # Extract the overall marginal I2 when available.
  # ----------------------------------------------------------------------------
  
  if (
    !is.null(model_summary$I2)
  ) {
    
    i2_val <- model_summary$I2[
      "marginal"
    ] * 100
    
    
    if (
      length(i2_val) == 0 ||
      is.na(i2_val)
    ) {
      
      i2_val <- model_summary$I2[
        1
      ] * 100
    }
    
  } else {
    
    i2_val <- NA_real_
  }
  
  
  # ============================================================================
  # 5.10.7 CALCULATE JOINT WALD TEST
  # ============================================================================
  
  wald <- calculate_joint_wald(
    model = fit_model,
    n_basis_coef = n_basis_coef
  )
  
  
  # ============================================================================
  # 5.10.8 EXTRACT NULL-MODEL FIT
  # ============================================================================
  
  null_aic <- AIC(
    null_model_ml
  )
  
  
  null_bic <- BIC(
    null_model_ml
  )
  
  
  null_loglik <- as.numeric(
    logLik(
      null_model_ml
    )
  )
  
  
  # ============================================================================
  # 5.10.9 CALCULATE DELTA AIC AND DELTA BIC
  # ============================================================================
  #
  # Negative values mean that the meta-regression has lower AIC/BIC than
  # the null model.
  #
  # For example:
  #
  #   Delta AIC = -10
  #
  # means that the meta-regression AIC is 10 points lower than the null.
  # ============================================================================
  
  
  delta_aic <- AIC(
    fit_model_ml
  ) - null_aic
  
  
  delta_bic <- BIC(
    fit_model_ml
  ) - null_bic
  
  
  # ============================================================================
  # 5.10.10 LIKELIHOOD-RATIO TEST
  # ============================================================================
  #
  # The univariable model introduces n_basis_coef additional fixed-effect
  # parameters.
  #
  # With the current three-parameter ONI cross-basis:
  #
  #   LR df = 3
  #
  # The test evaluates whether adding the district-level meta-regressor
  # improves the likelihood relative to the null model.
  # ============================================================================
  
  
  model_loglik_ml <- as.numeric(
    logLik(
      fit_model_ml
    )
  )
  
  
  LR_stat <- 2 * (
    model_loglik_ml -
      null_loglik
  )
  
  
  LR_df <- n_basis_coef
  
  
  LR_p <- pchisq(
    LR_stat,
    df = LR_df,
    lower.tail = FALSE
  )
  
  
  # ============================================================================
  # 5.10.11 FORMAT P VALUES
  # ============================================================================
  #
  # Numeric versions are retained for analysis.
  # Formatted versions are created only for presentation.
  # ============================================================================
  
  
  q_p_formatted <- if (
    !is.na(q_p)
  ) {
    
    format.pval(
      q_p,
      digits = 4,
      eps = 0.001
    )
    
  } else {
    
    "NA"
  }
  
  
  wald_p_formatted <- if (
    !is.na(wald$p_value)
  ) {
    
    format.pval(
      wald$p_value,
      digits = 4,
      eps = 0.001
    )
    
  } else {
    
    "NA"
  }
  
  
  LR_p_formatted <- if (
    !is.na(LR_p)
  ) {
    
    format.pval(
      LR_p,
      digits = 4,
      eps = 0.001
    )
    
  } else {
    
    "NA"
  }
  
  
  # ============================================================================
  # 5.10.12 ADD RESULTS TO MASTER TABLE
  # ============================================================================
  
  model_comparison_table <- rbind(
    
    model_comparison_table,
    
    data.frame(
      
      Covariate = as.character(
        cov
      ),
      
      N_districts = n_districts_model,
      
      
      # ------------------------------------------------------------------------
      # Original REML statistics
      # ------------------------------------------------------------------------
      
      AIC = round(
        as.numeric(aic_val)[1],
        2
      ),
      
      BIC = round(
        as.numeric(bic_val)[1],
        2
      ),
      
      LogLik = round(
        loglik_val,
        2
      ),
      
      Cochran_Q = round(
        q_stat,
        2
      ),
      
      Q_df = q_df,
      
      Q_p_value = q_p_formatted,
      
      I_squared_pct = round(
        i2_val,
        1
      ),
      
      Wald_test_stat = round(
        wald$statistic,
        2
      ),
      
      Wald_df = wald$df,
      
      Wald_p_value = wald_p_formatted,
      
      Wald_p_numeric = wald$p_value,
      
      
      # ------------------------------------------------------------------------
      # Null-model comparison
      # ------------------------------------------------------------------------
      
      Null_AIC = round(
        null_aic,
        2
      ),
      
      Null_BIC = round(
        null_bic,
        2
      ),
      
      Null_LogLik = round(
        null_loglik,
        2
      ),
      
      Delta_AIC = round(
        delta_aic,
        2
      ),
      
      Delta_BIC = round(
        delta_bic,
        2
      ),
      
      LR_stat = round(
        LR_stat,
        2
      ),
      
      LR_df = LR_df,
      
      LR_p_value = LR_p_formatted,
      
      LR_p_numeric = LR_p,
      
      stringsAsFactors = FALSE
    )
  )
  
  
  # ============================================================================
  # 5.10.13 GENERATE NON-LINEAR ONI–LEPTOSPIROSIS CURVES
  # ============================================================================
  #
  # The figure is generated directly inside the main analysis loop.
  #
  # The plotting operation itself is delegated to plot_modulation_curve().
  #
  # Nothing is saved to disk.
  #
  # If many meta-regressors are being evaluated, the plots will be displayed
  # sequentially in the active graphics device.
  # ============================================================================
  
  
  prediction_output <- plot_modulation_curve(
    
    fit_model = fit_model,
    
    model_data = loop_metadata,
    
    covariate_name = cov,
    
    basis_object = b_pred_oni,
    
    overall_prediction = pred_oni_ov,
    
    make_plot = TRUE
  )
  
  
  prediction_curve_list[[cov]] <- prediction_output$predictions
  
  
  # ============================================================================
  # 5.10.14 PRINT CURRENT MODEL SUMMARY
  # ============================================================================
  
  cat("\n")
  cat(
    "Covariate:",
    cov,
    "\n"
  )
  
  cat(
    "N districts:",
    n_districts_model,
    "\n"
  )
  
  cat(
    "REML AIC:",
    round(aic_val, 2),
    "\n"
  )
  
  cat(
    "Null ML AIC:",
    round(null_aic, 2),
    "\n"
  )
  
  cat(
    "Delta AIC:",
    round(delta_aic, 2),
    "\n"
  )
  
  cat(
    "REML BIC:",
    round(bic_val, 2),
    "\n"
  )
  
  cat(
    "Null ML BIC:",
    round(null_bic, 2),
    "\n"
  )
  
  cat(
    "Delta BIC:",
    round(delta_bic, 2),
    "\n"
  )
  
  cat(
    "Likelihood-ratio statistic:",
    round(LR_stat, 2),
    "\n"
  )
  
  cat(
    "LR df:",
    LR_df,
    "\n"
  )
  
  cat(
    "LR P:",
    LR_p_formatted,
    "\n"
  )
  
  cat(
    "Joint Wald statistic:",
    round(wald$statistic, 2),
    "\n"
  )
  
  cat(
    "Wald df:",
    wald$df,
    "\n"
  )
  
  cat(
    "Wald P:",
    wald_p_formatted,
    "\n"
  )
  
}


# ==============================================================================
# 5.11 ORDER RESULTS
# ==============================================================================
#
# The table is ordered according to evidence that the meta-regression improves
# model fit relative to the corresponding null model.
# ==============================================================================


model_comparison_table <- model_comparison_table %>%
  arrange(
    LR_p_numeric,
    Delta_AIC
  )


# ==============================================================================
# 5.12 DISPLAY COMPLETE MODEL-COMPARISON TABLE
# ==============================================================================

cat("\n")
cat("============================================================\n")
cat("COMPLETE SECOND-STAGE META-REGRESSION RESULTS\n")
cat("============================================================\n")

print(
  model_comparison_table,
  row.names = FALSE
)


# ==============================================================================
# 5.13 CREATE A CLEAN PUBLICATION-ORIENTED SUMMARY TABLE
# ==============================================================================
#
# This table contains the original diagnostics plus the newly added
# null-model comparison statistics.
# ==============================================================================


final_univariable_table <- model_comparison_table %>%
  select(
    Covariate,
    N_districts,
    
    # Original model
    AIC,
    BIC,
    LogLik,
    Cochran_Q,
    Q_df,
    Q_p_value,
    I_squared_pct,
    
    Wald_test_stat,
    Wald_df,
    Wald_p_value,
    
    # Null comparison
    Null_AIC,
    Null_BIC,
    Null_LogLik,
    
    Delta_AIC,
    Delta_BIC,
    
    LR_stat,
    LR_df,
    LR_p_value
  )


cat("\n")
cat("============================================================\n")
cat("PUBLICATION-ORIENTED UNIVARIABLE MODEL COMPARISON\n")
cat("============================================================\n")

print(
  final_univariable_table,
  row.names = FALSE
)


# ==============================================================================
# 5.14 EXPORT NUMERICAL RESULTS
# ==============================================================================
#
# The figures are intentionally NOT saved.
#
# The statistical results are saved for reproducibility and later table
# preparation.
# ==============================================================================


write.csv(
  model_comparison_table,
  
  "~/Documents/Heidelberg/ESA_climate/analysis/ONI_Leptospirosis_second_stage_meta_regression_diagnostics.csv",
  
  row.names = FALSE
)


write.csv(
  final_univariable_table,
  
  "~/Documents/Heidelberg/ESA_climate/analysis/ONI_Leptospirosis_univariable_meta_regression_vs_null.csv",
  
  row.names = FALSE
)


# ==============================================================================
# 5.15 IDENTIFY CANDIDATES FOR MULTIVARIABLE META-REGRESSION
# ==============================================================================
#
# The univariable analysis is a screening and heterogeneity-explanation stage.
#
# A liberal threshold can be used to identify variables for consideration in
# a multivariable model.
#
# Importantly, this is NOT equivalent to declaring these variables as
# independent predictors.
#
# Given the relatively small number of districts, the multivariable model
# must remain parsimonious.
#
# In particular:
#
#   - Multiple LCCS proportions should not automatically be entered together
#     because land-cover proportions are compositional.
#
#   - HSI 2016, HSI 2025 and HSI change are related representations of the
#     same underlying construct.
#
#   - Rainfall mean and rainfall change may be correlated.
#
#   - NDVI mean and NDVI change may be correlated.
#
#   - Soil-moisture mean and soil-moisture change may be correlated.
#
# Therefore, candidate selection should combine statistical evidence with
# scientific plausibility and correlation/collinearity assessment.
# ==============================================================================


candidate_covariates <- model_comparison_table %>%
  filter(
    LR_p_numeric < 0.10
  ) %>%
  pull(
    Covariate
  )


cat("\n")
cat("============================================================\n")
cat("CANDIDATE META-REGRESSORS FOR MULTIVARIABLE ANALYSIS\n")
cat("============================================================\n")

print(
  candidate_covariates
)


# ==============================================================================
# 5.16 CORRELATION MATRIX OF CANDIDATE META-REGRESSORS
# ==============================================================================
#
# Correlation is examined before constructing a multivariable model.
#
# This is particularly important because the second stage contains only a
# limited number of districts.
# ==============================================================================


if (
  length(candidate_covariates) > 1
) {
  
  candidate_data <- loop_metadata_base %>%
    filter(
      keep_oni
    ) %>%
    select(
      district,
      all_of(candidate_covariates)
    )
  
  
  candidate_matrix <- candidate_data %>%
    select(
      -district
    )
  
  
  correlation_matrix <- cor(
    candidate_matrix,
    use = "pairwise.complete.obs"
  )
  
  
  cat("\n")
  cat("============================================================\n")
  cat("CORRELATION MATRIX OF CANDIDATE META-REGRESSORS\n")
  cat("============================================================\n")
  
  print(
    round(
      correlation_matrix,
      2
    )
  )
  
  
  write.csv(
    correlation_matrix,
    
    "~/Documents/Heidelberg/ESA_climate/analysis/ONI_Leptospirosis_meta_regressor_correlation_matrix.csv"
  )
  
} else {
  
  correlation_matrix <- NULL
  
  cat(
    "\nFewer than two candidate covariates available; correlation matrix not calculated.\n"
  )
}


# ==============================================================================
# 5.17 OPTIONAL MULTIVARIABLE META-REGRESSION
# ==============================================================================
#
# The multivariable model is deliberately NOT automatically fitted using every
# variable that passes the univariable screening threshold.
#
# Instead, selected_covariates must be specified based on:
#
#   1. Biological plausibility
#   2. Results of univariable null-model comparisons
#   3. Correlation structure
#   4. Avoidance of redundant variables
#   5. Number of available districts
#   6. Stability of model estimation
#
# Example:
#
# selected_covariates <- c(
#   "lccs_190",
#   "soil_moisture_mean_2007_2025"
# )
#
# Replace this example with the scientifically justified specification.
# ==============================================================================


selected_covariates <- c(
  "lccs_200",
  "soil_moisture_mean_2007_2025",
  "ndvi_change_2007_2025", 
  "lccs_110",
  "lccs_130",
  "lccs_150",
  "lccs_210"
)


# ==============================================================================
# 5.18 MULTIVARIABLE MODEL DATASET
# ==============================================================================

if (
  length(selected_covariates) > 0
) {
  
  multi_data <- loop_metadata_base %>%
    filter(
      keep_oni
    ) %>%
    select(
      district,
      district_pcode,
      first_stage_row,
      all_of(selected_covariates)
    ) %>%
    filter(
      complete.cases(.)
    )
  
  
  multi_rows <- multi_data$first_stage_row
  
  
  cat("\n")
  cat("============================================================\n")
  cat("MULTIVARIABLE META-REGRESSION\n")
  cat("============================================================\n")
  
  cat(
    "Selected covariates:",
    paste(
      selected_covariates,
      collapse = ", "
    ),
    "\n"
  )
  
  cat(
    "Number of districts:",
    nrow(multi_data),
    "\n"
  )
  
  
  # ============================================================================
  # 5.18.1 CONSTRUCT MULTIVARIABLE FORMULA
  # ============================================================================
  
  multi_formula <- as.formula(
    paste(
      "coef_oni_all[multi_rows, ] ~",
      paste(
        selected_covariates,
        collapse = " + "
      )
    )
  )
  
  
  # ============================================================================
  # 5.18.2 NULL MODEL ON THE SAME MULTIVARIABLE SAMPLE
  # ============================================================================
  
  multi_null_ml <- mvmeta(
    coef_oni_all[
      multi_rows,
      ,
      drop = FALSE
    ] ~ 1,
    
    S = vcov_oni_all[
      multi_rows
    ],
    
    data = multi_data,
    
    method = "ml"
  )
  
  
  # ============================================================================
  # 5.18.3 MULTIVARIABLE MODEL USING ML
  # ============================================================================
  
  multivariable_model_ml <- mvmeta(
    multi_formula,
    
    S = vcov_oni_all[
      multi_rows
    ],
    
    data = multi_data,
    
    method = "ml"
  )
  
  
  # ============================================================================
  # 5.18.4 MULTIVARIABLE MODEL USING REML
  # ============================================================================
  #
  # This is the final estimation model after the fixed-effects structure has
  # been selected.
  # ============================================================================
  
  multivariable_model_reml <- mvmeta(
    multi_formula,
    
    S = vcov_oni_all[
      multi_rows
    ],
    
    data = multi_data,
    
    method = "reml"
  )
  
  
  # ============================================================================
  # 5.18.5 MULTIVARIABLE MODEL FIT COMPARISON
  # ============================================================================
  
  multi_null_logLik <- as.numeric(
    logLik(
      multi_null_ml
    )
  )
  
  
  multi_model_logLik <- as.numeric(
    logLik(
      multivariable_model_ml
    )
  )
  
  
  multi_LR_stat <- 2 * (
    multi_model_logLik -
      multi_null_logLik
  )
  
  
  # Each additional covariate contributes n_basis_coef fixed-effect parameters.
  
  multi_LR_df <- length(
    selected_covariates
  ) * n_basis_coef
  
  
  multi_LR_p <- pchisq(
    multi_LR_stat,
    df = multi_LR_df,
    lower.tail = FALSE
  )
  
  
  # ============================================================================
  # 5.18.6 DISPLAY MULTIVARIABLE MODEL RESULTS
  # ============================================================================
  
  cat("\n")
  cat("Multivariable ML model:\n")
  
  cat(
    "AIC:",
    round(
      AIC(multivariable_model_ml),
      2
    ),
    "\n"
  )
  
  cat(
    "BIC:",
    round(
      BIC(multivariable_model_ml),
      2
    ),
    "\n"
  )
  
  cat(
    "LogLik:",
    round(
      multi_model_logLik,
      2
    ),
    "\n"
  )
  
  cat(
    "Null-model LogLik:",
    round(
      multi_null_logLik,
      2
    ),
    "\n"
  )
  
  cat(
    "LR statistic:",
    round(
      multi_LR_stat,
      2
    ),
    "\n"
  )
  
  cat(
    "LR df:",
    multi_LR_df,
    "\n"
  )
  
  cat(
    "LR P:",
    format.pval(
      multi_LR_p,
      digits = 4,
      eps = 0.001
    ),
    "\n"
  )
  
  
  # ============================================================================
  # 5.18.7 FINAL REML MODEL SUMMARY
  # ============================================================================
  
  cat("\n")
  cat("============================================================\n")
  cat("FINAL MULTIVARIABLE REML MODEL\n")
  cat("============================================================\n")
  
  print(
    summary(
      multivariable_model_reml
    )
  )
  
  
  # ============================================================================
  # 5.18.8 SAVE MODEL OBJECTS
  # ============================================================================
  
  saveRDS(
    multivariable_model_ml,
    
    "~/Documents/Heidelberg/ESA_climate/analysis/ONI_Leptospirosis_multivariable_meta_regression_ML.rds"
  )
  
  
  saveRDS(
    multivariable_model_reml,
    
    "~/Documents/Heidelberg/ESA_climate/analysis/ONI_Leptospirosis_multivariable_meta_regression_REML.rds"
  )
  
  
} else {
  
  cat("\n")
  cat("No multivariable model fitted.\n")
  cat(
    "Specify selected_covariates after reviewing the univariable results,\n"
  )
  cat(
    "correlation matrix, and scientific plausibility.\n"
  )
}


# ==============================================================================
# 5.19 FINAL ANALYTICAL INTERPRETATION
# ==============================================================================
#
# The second-stage analysis should be interpreted as follows:
#
# The first-stage models estimate district-specific nonlinear ONI–
# leptospirosis association profiles. The second-stage multivariate
# random-effects meta-analysis pools these coefficient vectors while allowing
# for between-district heterogeneity.
#
# District-level environmental and landscape characteristics are then entered
# as meta-regressors to assess whether this heterogeneity is systematically
# related to differences between districts.
#
# For each meta-regressor:
#
#   1. The REML meta-regression provides the primary model diagnostics,
#      including AIC, BIC, log-likelihood, Cochran's Q, I-squared and the
#      joint Wald test.
#
#   2. A corresponding intercept-only null model is fitted on the same
#      analytical districts.
#
#   3. The meta-regression is compared with the null using ML-based
#      likelihood-ratio testing.
#
#   4. Delta AIC and Delta BIC quantify the change in model fit after
#      introducing the district-level characteristic.
#
#   5. The joint Wald test evaluates whether the covariate is associated with
#      the complete vector of ONI cross-basis coefficients.
#
#   6. Predicted ONI–leptospirosis curves at the 10th, 50th and 90th
#      percentiles of the meta-regressor provide a graphical representation
#      of how the shape and magnitude of the ONI association vary across the
#      environmental gradient.
#
#
# The analysis therefore addresses TWO related questions:
#
#   QUESTION 1:
#
#       Is there substantial residual between-district heterogeneity in the
#       ONI–leptospirosis association?
#
#       -> Cochran's Q and I-squared
#
#
#   QUESTION 2:
#
#       Does a specific district-level characteristic explain some of this
#       heterogeneity?
#
#       -> Likelihood-ratio test
#       -> Delta AIC
#       -> Delta BIC
#       -> Joint Wald test
#       -> Predicted modulation curves
#
#
# A statistically significant meta-regression should therefore be described
# as evidence that the district-specific ONI–leptospirosis association varies
# systematically according to the corresponding district characteristic.
#
# It should NOT automatically be interpreted as evidence of a causal
# interaction or causal effect modification.
#
# ==============================================================================
# END
# ==============================================================================

