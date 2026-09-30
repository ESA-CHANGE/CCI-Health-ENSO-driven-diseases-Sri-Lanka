# ==============================================================================
# ESA-CHANGE | Leptospirosis–Climate Analysis
# ESA-CHANGE | Use Case 6
# ENSO-Driven Climate Variability and Dengue & Leptospirosis in Sri Lanka 
#===============================================================================
# Primary two-stage DLNM model for ENSO and Indian Ocean Dipole
#
# Description:
# This script estimates district-specific associations between monthly
# leptospirosis incidence and two large-scale climate indices:
#   1. Oceanic Niño Index (ONI)
#   2. Dipole Mode Index (DMI)
#
# A first-stage distributed lag non-linear model (DLNM) is fitted separately
# for each district, followed by multivariate random-effects meta-analysis
# to obtain national pooled exposure-response associations.
#
# Input:
#   analysis_df
#     Required variables:
#       district, year, month, time, Leptospirosis, Population, ONI, DMI
#
# Outputs:
#   - District-level model diagnostics
#   - National pooled ONI and DMI exposure-response estimates
#   - Lag-specific pooled ONI and DMI estimates
#   - Overall and lag-specific PDF figures
#
# Methods:
#   First stage: quasipoisson GLM with cross-basis functions for ONI and DMI,
#                seasonal adjustment, long-term temporal adjustment, and
#                population offset.
#   Second stage: multivariate random-effects meta-analysis using REML.
#
# Reference:
#   Gasparrini A, Armstrong B, Kenward MG. Distributed lag non-linear models.
#   Statistics in Medicine. 2010;29:2224–2234.
# ==============================================================================


# ------------------------------------------------------------------------------
# 1. Packages
# ------------------------------------------------------------------------------

library(dplyr)
library(dlnm)
library(mvmeta)
library(splines)

# ------------------------------------------------------------------------------
# 2. Analysis configuration
# ------------------------------------------------------------------------------

lag_limit <- c(0, 6)
lags_to_meta <- c(0, 2, 4, 6)

districts <- unique(analysis_df$district)
n_dist <- length(districts)

# ONI range and exposure-response knots
oni_min <- -1.60
oni_max <- 2.80
oni_knots <- c(-0.5, 0.5)

# DMI range and exposure-response knots
dmi_min <- -0.758
dmi_max <- 0.964
dmi_knots <- c(-0.100, 0.268)


# ------------------------------------------------------------------------------
# 3. QAIC function
# ------------------------------------------------------------------------------

# Quasi-AIC calculated using the Poisson log-likelihood evaluated at the
# fitted mean, with the model-specific dispersion parameter and number of
# estimated coefficients.

fqaic <- function(model) {
  
  loglik <- sum(
    dpois(model$y, model$fitted.values, log = TRUE)
  )
  
  dispersion <- summary(model)$dispersion
  k <- model$rank
  
  -2 * loglik + 2 * k * dispersion
}


# ------------------------------------------------------------------------------
# 4. Prepare analysis data
# ------------------------------------------------------------------------------

# Read the CSV file into R
analysis_df <- read.csv(
  "path/data_Case_6.csv",
  stringsAsFactors = FALSE
)

final_analysis_data <- analysis_df %>%
  arrange(district, year, month) %>%
  filter(
    !is.na(ONI),
    !is.na(DMI)
  )


# ------------------------------------------------------------------------------
# 5. Pre-allocate storage for first-stage estimates
# ------------------------------------------------------------------------------

# Overall cumulative exposure-response coefficients
coef_oni_all <- matrix(
  NA_real_,
  nrow = n_dist,
  ncol = 3,
  dimnames = list(districts, paste0("oni_b", 1:3))
)

vcov_oni_all <- vector("list", n_dist)

coef_dmi_all <- matrix(
  NA_real_,
  nrow = n_dist,
  ncol = 3,
  dimnames = list(districts, paste0("dmi_b", 1:3))
)

vcov_dmi_all <- vector("list", n_dist)


# Lag-specific coefficient and covariance storage
results_oni_by_lag <- setNames(
  lapply(lags_to_meta, function(x) list()),
  paste0("lag", lags_to_meta)
)

results_dmi_by_lag <- setNames(
  lapply(lags_to_meta, function(x) list()),
  paste0("lag", lags_to_meta)
)

first_stage_summaries <- list()


# ------------------------------------------------------------------------------
# 6. First-stage district-specific DLNM models
# ------------------------------------------------------------------------------

message("Fitting first-stage district-specific DLNMs...")

for (i in seq_along(districts)) {
  
  district_id <- districts[i]
  
  sub <- final_analysis_data %>%
    filter(district == district_id)
  
  if (nrow(sub) < 10) {
    next
  }
  
  converged_flag <- TRUE
  
  tryCatch({
    
    # --------------------------------------------------------------------------
    # 6.1 ONI cross-basis
    # --------------------------------------------------------------------------
    
    cb_oni <- crossbasis(
      sub$ONI,
      lag = lag_limit,
      argvar = list(
        fun = "ns",
        knots = oni_knots,
        Boundary.knots = c(oni_min, oni_max)
      ),
      arglag = list(
        fun = "ns",
        knots = c(2, 4)
      )
    )
    
    
    # --------------------------------------------------------------------------
    # 6.2 DMI cross-basis
    # --------------------------------------------------------------------------
    
    cb_dmi <- crossbasis(
      sub$DMI,
      lag = lag_limit,
      argvar = list(
        fun = "ns",
        knots = dmi_knots,
        Boundary.knots = c(dmi_min, dmi_max)
      ),
      arglag = list(
        fun = "ns",
        knots = c(2, 4)
      )
    )
    
    
    # --------------------------------------------------------------------------
    # 6.3 Co-adjusted regression model
    # --------------------------------------------------------------------------
    
    model <- glm(
      Leptospirosis ~ cb_oni + cb_dmi +
        ns(time, df = 15) +
        factor(month) +
        offset(log(Population)),
      family = quasipoisson(),
      data = sub
    )
    
    if (!model$converged) {
      converged_flag <- FALSE
    }
    
    
    # --------------------------------------------------------------------------
    # 6.4 Overall cumulative associations
    # --------------------------------------------------------------------------
    
    reduced_oni <- crossreduce(
      cb_oni,
      model,
      type = "overall",
      cen = 0,
      by = 0.1
    )
    
    reduced_dmi <- crossreduce(
      cb_dmi,
      model,
      type = "overall",
      cen = 0,
      by = 0.05
    )
    
    coef_oni_all[i, ] <- coef(reduced_oni)
    vcov_oni_all[[i]] <- vcov(reduced_oni)
    
    coef_dmi_all[i, ] <- coef(reduced_dmi)
    vcov_dmi_all[[i]] <- vcov(reduced_dmi)
    
    
    # --------------------------------------------------------------------------
    # 6.5 Lag-specific associations
    # --------------------------------------------------------------------------
    
    for (lag_value in lags_to_meta) {
      
      lag_label <- paste0("lag", lag_value)
      
      reduced_oni_lag <- crossreduce(
        cb_oni,
        model,
        type = "lag",
        value = lag_value,
        cen = 0
      )
      
      results_oni_by_lag[[lag_label]][[district_id]] <- list(
        coef = coef(reduced_oni_lag),
        vcov = vcov(reduced_oni_lag)
      )
      
      
      reduced_dmi_lag <- crossreduce(
        cb_dmi,
        model,
        type = "lag",
        value = lag_value,
        cen = 0
      )
      
      results_dmi_by_lag[[lag_label]][[district_id]] <- list(
        coef = coef(reduced_dmi_lag),
        vcov = vcov(reduced_dmi_lag)
      )
    }
    
    
    # --------------------------------------------------------------------------
    # 6.6 District-level model diagnostics
    # --------------------------------------------------------------------------
    
    model_summary <- summary(model)
    
    first_stage_summaries[[district_id]] <- data.frame(
      District = district_id,
      Observations = nrow(sub),
      Null_Deviance = model$null.deviance,
      Residual_Deviance = model$deviance,
      Residual_DF = model$df.residual,
      Parameters_K = model$rank,
      Dispersion_Phi = model_summary$dispersion,
      District_QAIC = fqaic(model),
      Converged = ifelse(converged_flag, "Yes", "No"),
      stringsAsFactors = FALSE
    )
    
  }, error = function(e) {
    
    first_stage_summaries[[district_id]] <<- data.frame(
      District = district_id,
      Observations = nrow(sub),
      Null_Deviance = NA_real_,
      Residual_Deviance = NA_real_,
      Residual_DF = NA_real_,
      Parameters_K = NA_real_,
      Dispersion_Phi = NA_real_,
      District_QAIC = NA_real_,
      Converged = "ERROR",
      stringsAsFactors = FALSE
    )
  })
}


# ------------------------------------------------------------------------------
# 7. First-stage model diagnostics
# ------------------------------------------------------------------------------

stage1_summary <- bind_rows(first_stage_summaries)

message("First-stage modelling completed.")

print(stage1_summary)

message(
  "Districts evaluated: ",
  nrow(stage1_summary)
)

message(
  "Non-converged/error models: ",
  sum(stage1_summary$Converged != "Yes", na.rm = TRUE)
)

message(
  "Sum of district QAIC: ",
  round(sum(stage1_summary$District_QAIC, na.rm = TRUE), 2)
)


# ------------------------------------------------------------------------------
# 8. Second-stage national meta-analysis
# ------------------------------------------------------------------------------

message("Fitting national random-effects meta-analyses...")


# Overall ONI association
keep_oni <- !is.na(coef_oni_all[, 1])

meta_oni_overall <- mvmeta(
  coef_oni_all[keep_oni, , drop = FALSE] ~ 1,
  S = vcov_oni_all[keep_oni],
  method = "reml"
)


# Overall DMI association
keep_dmi <- !is.na(coef_dmi_all[, 1])

meta_dmi_overall <- mvmeta(
  coef_dmi_all[keep_dmi, , drop = FALSE] ~ 1,
  S = vcov_dmi_all[keep_dmi],
  method = "reml"
)

# Examine model summaries 

print(summary(meta_oni_overall))
print(summary(meta_dmi_overall))

# ------------------------------------------------------------------------------
# 9. National pooled exposure-response predictions
# ------------------------------------------------------------------------------

# ONI
x_oni_seq <- seq(oni_min, oni_max, by = 0.1)
b_pred_oni <- onebasis(x_oni_seq, fun = "ns", knots = oni_knots, Boundary.knots = c(oni_min, oni_max))
pred_oni_ov <- crosspred(b_pred_oni, coef = coef(meta_oni_overall), vcov = vcov(meta_oni_overall), model.link = "log", cen = 0, by = 0.1)

# DMI
x_dmi_seq <- seq(dmi_min, dmi_max, by = 0.05)
b_pred_dmi <- onebasis(x_dmi_seq, fun = "ns", knots = dmi_knots, Boundary.knots = c(dmi_min, dmi_max))
pred_dmi_ov <- crosspred(b_pred_dmi, coef = coef(meta_dmi_overall), vcov = vcov(meta_dmi_overall), model.link = "log", cen = 0, by = 0.05)

# ------------------------------------------------------------------------------
# 10. Lag-specific second-stage meta-analysis
# ------------------------------------------------------------------------------

meta_oni_lag <- list()
meta_dmi_lag <- list()


mv_lag_oni_results <- list()
mv_lag_dmi_results <- list()

for(l in lags_to_meta) {
  lag_lab <- paste0("lag", l)
  
  # --- Compile ONI Lag Parameters ---
  list_oni <- results_oni_by_lag[[lag_lab]]
  list_oni <- list_oni[!sapply(list_oni, is.null)]
  if(length(list_oni) >= 3) {
    c_mat <- do.call(rbind, lapply(list_oni, function(x) x$coef))
    k <- apply(c_mat, 1, function(x) all(abs(x) < 10))
    if(sum(k) >= 3) {
      v_list <- lapply(list_oni[k], function(x) { m <- x$vcov; dimnames(m) <- list(paste0("b",1:3), paste0("b",1:3)); m })
      mv_lag_oni_results[[lag_lab]] <- mvmeta(c_mat[k, , drop=FALSE] ~ 1, S = v_list, method = "reml")
    }
  }
  
  # --- Compile DMI Lag Parameters ---
  list_dmi <- results_dmi_by_lag[[lag_lab]]
  list_dmi <- list_dmi[!sapply(list_dmi, is.null)]
  if(length(list_dmi) >= 3) {
    c_mat <- do.call(rbind, lapply(list_dmi, function(x) x$coef))
    k <- apply(c_mat, 1, function(x) all(abs(x) < 10))
    if(sum(k) >= 3) {
      v_list <- lapply(list_dmi[k], function(x) { m <- x$vcov; dimnames(m) <- list(paste0("b",1:3), paste0("b",1:3)); m })
      mv_lag_dmi_results[[lag_lab]] <- mvmeta(c_mat[k, , drop=FALSE] ~ 1, S = v_list, method = "reml")
    }
  }
}


# ------------------------------------------------------------------------------
# 11.1 Visualization pooled overall cumulative association
# ------------------------------------------------------------------------------
# --- Plot 4.1: Adjusted ONI Effects ---
plot(pred_oni_ov$predvar, pred_oni_ov$allRRfit, type = "n", ylim = c(.2, 6), las = 1, bty = "l",
     xlab = "Oceanic Niño Index (ONI)", ylab = "Relative Risk", main = "ONI Cumulative Risk")
polygon(c(pred_oni_ov$predvar, rev(pred_oni_ov$predvar)), c(pred_oni_ov$allRRlow, rev(pred_oni_ov$allRRhigh)),
        col = rgb(0.1, 0.1, 0.8, 0.2), border = NA)
abline(h = 1, lty = 2, col = "grey40")
lines(pred_oni_ov$predvar, pred_oni_ov$allRRfit, col = "darkblue", lwd = 3)
# Reference lines
abline(h = 1, lty = 2, col = "grey40") # Null risk baseline
# Truncated Vertical Anomaly Threshold Lines (stopping exactly at y = 3.5)
lines(c(-0.5, -0.5), c(0.2, 5), lty = 3, col = "darkgrey", lwd = 1.5)
lines(c(0.5, 0.5), c(0.2, 5), lty = 3, col = "darkgrey", lwd = 1.5)
# Left-pointing arrow at x = -0.5 for La Niña (code = 2 indicates arrow head at the end point)
arrows(x0 = -0.5, y0 = 4, x1 = -1.2, y1 = 4, length = 0.1, lwd = 1.5, col = "darkblue", code = 2)
text(x = -0.85, y = 4.3, labels = "La Niña", col = "darkblue", cex = 0.8, font = 2)

# Right-pointing arrow at x = 0.5 for El Niño
arrows(x0 = 0.5, y0 = 4, x1 = 1.2, y1 = 4, length = 0.1, lwd = 1.5, col = "darkred", code = 2)
text(x = 0.85, y = 4.3, labels = "El Niño", col = "darkred", cex = 0.8, font = 2)
# Legend placement
legend("topleft", legend = c("National Pooled Average", "95% Confidence Interval"),
       col = c("darkblue", NA), lwd = c(3, NA), 
       fill = c(NA, rgb(0.1, 0.1, 0.8, 0.15)), border = c(NA, NA),
       bty = "n", cex = 0.9, y.intersp = 1.2)



# --- Plot 4.2: Adjusted DMI Effects ---
plot(pred_dmi_ov$predvar, pred_dmi_ov$allRRfit, type = "n", ylim = c(.2, 6), las = 1, bty = "l",
     xlab = "Dipole Mode Index (DMI)", ylab = "Relative Risk", main = "DMI Cumulative Risk")
polygon(c(pred_dmi_ov$predvar, rev(pred_dmi_ov$predvar)), c(pred_dmi_ov$allRRlow, rev(pred_dmi_ov$allRRhigh)),
        col = rgb(0.1, 0.8, 0.1, 0.2), border = NA) # Green shading for DMI
abline(h = 1, lty = 2, col = "grey40")
# Truncated Vertical Anomaly Threshold Lines at neutral zero boundaries
# Standard IOD thresholds typically look at -0.4 and +0.4 bounds
lines(c(-0.4, -0.4), c(0.2, 5), lty = 3, col = "darkgrey", lwd = 1.5)
lines(c(0.4, 0.4), c(0.2, 5), lty = 3, col = "darkgrey", lwd = 1.5)
# Directional Operational Vectors for Indian Ocean Dynamics at y = 2.5
# Left-pointing arrow for Negative IOD phase
arrows(x0 = -0.4, y0 = 4, x1 = -0.75, y1 = 4, length = 0.08, lwd = 1.5, col = "navyblue", code = 2)
text(x = -0.58, y = 4.3, labels = "Neg IOD", col = "navyblue", cex = 0.75, font = 2)

# Right-pointing arrow for Positive IOD phase
arrows(x0 = 0.4, y0 = 4, x1 = 0.75, y1 = 4, length = 0.08, lwd = 1.5, col = "darkgoldenrod", code = 2)
text(x = 0.58, y = 4.3, labels = "Pos IOD", col = "darkgoldenrod", cex = 0.75, font = 2)

lines(pred_dmi_ov$predvar, pred_dmi_ov$allRRfit, col = "darkgreen", lwd = 3)
# Legend placement
# Legend placement
legend("topleft", legend = c("National Pooled Average", "95% Confidence Interval"),
       col = c("darkgreen", NA), lwd = c(3, NA), 
       fill = c(NA, rgb(0.1, 0.8, 0.1, 0.15)), border = c(NA, NA),
       bty = "n", cex = 0.9, y.intersp = 1.2)
       
# ------------------------------------------------------------------------------
# 11.1 Visualization lag-sepcific summaries
# ------------------------------------------------------------------------------
# Generates a Single-Column Matrix of Lag-Specific Risk Maps for ONI
# Adjusts mfrow to c(4, 1) for a stacked, single-column configuration
# ==============================================================================


par(mfrow = c(2, 2), mar = c(4, 4.5, 2.5, 2), oma = c(0, 0, 2, 0))

for(l in lags_to_meta) {
  lag_lab <- paste0("lag", l)
  
  # Isolate and extract ONI Tracks
  m_oni <- mv_lag_oni_results[[lag_lab]]
  
  if(!is.null(m_oni)) {
    # Calculate cross-predictions for the specific lag layer
    p_oni <- crosspred(b_pred_oni, coef = coef(m_oni), vcov = vcov(m_oni), 
                       model.link = "log", cen = 0)
    
    # Initialize blank plot frame with y-axis log-bounds
    plot(p_oni$predvar, p_oni$allRRfit, type = "n", ylim = c(0.2, 4), las = 1, 
         bty = "l", xlab = "Oceanic Niño Index (ONI)", ylab = "Relative Risk (RR)", 
         main = paste(" Lag", l, "Months"))
    
    # Render 95% Confidence Interval Band (Subtle Red/Firebrick shading)
    polygon(c(p_oni$predvar, rev(p_oni$predvar)), c(p_oni$allRRlow, rev(p_oni$allRRhigh)), 
            col = rgb(0.8, 0.2, 0.2, 0.1), border = NA)
    
    # Draw horizontal null risk baseline (RR = 1.0)
    abline(h = 1, lty = 2, col = "grey40")
    
    # Truncated Vertical Anomaly Threshold Lines (stopping exactly at y = 6)
    lines(c(-0.5, -0.5), c(0.2, 3), lty = 3, col = "darkgrey", lwd = 1.5)
    lines(c(0.5, 0.5), c(0.2, 3), lty = 3, col = "darkgrey", lwd = 1.5)
    
    # Directional Operational Vectors at y = 4
    # Left-pointing arrow at x = -0.5 for La Niña
    arrows(x0 = -0.5, y0 = 2.5, x1 = -1.2, y1 = 2.5, length = 0.08, lwd = 1.5, col = "darkblue", code = 2)
    text(x = -0.85, y = 2.7, labels = "La Niña", col = "darkblue", cex = 0.75, font = 2)
    
    # Right-pointing arrow at x = 0.5 for El Niño
    arrows(x0 = 0.5, y0 = 2.5, x1 = 1.2, y1 = 2.5, length = 0.08, lwd = 1.5, col = "darkred", code = 2)
    text(x = 0.85, y = 2.7, labels = "El Niño", col = "darkred", cex = 0.75, font = 2)
    
    # Plot the pooled lag-specific effect estimate line
    lines(p_oni$predvar, p_oni$allRRfit, col = "firebrick", lwd = 2.5)
    
    # Legend placement in top left corner of each plot
    legend("topleft", legend = c(paste("Pooled Average (Lag", l, ")"), "95% Confidence Interval"),
           col = c("firebrick", NA), lwd = c(2.5, NA), 
           fill = c(NA, rgb(0.8, 0.2, 0.2, 0.1)), border = c(NA, NA),
           bty = "n", cex = 0.85, y.intersp = 1.1)
    
  } else {
    # Gracefully skip and denote lags with missing data fields
    plot.new()
    title(main = paste("ONI Lag", l, "— (Insufficient Data Available)"))
  }
}

# Add overall master title above the column stack
mtext("Lag-Specific Association Between ONI and Leptospirosis", 
      outer = TRUE, cex = 1.2, font = 2)


# ==============================================================================

# ==============================================================================
# RENDER MATRIX: LAG-SPECIFIC ASSOCIATIONS FOR DMI
# ==============================================================================

# Configure a 2x2 multi-frame layout to match the ONI matrix structural layout
par(mfrow = c(2, 2), mar = c(4, 4.5, 2.5, 2), oma = c(0, 0, 2, 0))

for(l in lags_to_meta) {
  lag_lab <- paste0("lag", l)
  
  # Isolate and extract DMI Tracks
  m_dmi <- mv_lag_dmi_results[[lag_lab]]
  
  if(!is.null(m_dmi)) {
    # Calculate cross-predictions for the specific DMI lag layer
    # Note: b_pred_dmi must be predefined in your workspace matching DMI ranges
    p_dmi <- crosspred(b_pred_dmi, coef = coef(m_dmi), vcov = vcov(m_dmi), 
                       model.link = "log", cen = 0)
    
    # Initialize blank plot frame with y-axis bounds set to match your main DMI plot scale
    plot(p_dmi$predvar, p_dmi$allRRfit, type = "n", ylim = c(0.6, 2), las = 1, 
         bty = "l", xlab = "Dipole Mode Index (DMI)", ylab = "Relative Risk (RR)", 
         main = paste(" Lag", l, "Months"))
    
    # Render 95% Confidence Interval Band (Subtle Forest Green shading)
    polygon(c(p_dmi$predvar, rev(p_dmi$predvar)), c(p_dmi$allRRlow, rev(p_dmi$allRRhigh)), 
            col = rgb(0.13, 0.55, 0.13, 0.1), border = NA)
    
    # Draw horizontal null risk baseline reference line (RR = 1.0)
    abline(h = 1, lty = 2, col = "grey40")
    
    # Truncated Vertical Anomaly Threshold Lines at neutral zero boundaries
    # Standard IOD thresholds typically look at -0.4 and +0.4 bounds
    lines(c(-0.4, -0.4), c(0.2, 1.5), lty = 3, col = "darkgrey", lwd = 1.5)
    lines(c(0.4, 0.4), c(0.2, 1.5), lty = 3, col = "darkgrey", lwd = 1.5)
    
    # Directional Operational Vectors for Indian Ocean Dynamics at y = 2.5
    # Left-pointing arrow for Negative IOD phase
    arrows(x0 = -0.4, y0 = 1.6, x1 = -0.75, y1 = 1.6, length = 0.08, lwd = 1.5, col = "navyblue", code = 2)
    text(x = -0.58, y = 1.65, labels = "Neg IOD", col = "navyblue", cex = 0.75, font = 2)
    
    # Right-pointing arrow for Positive IOD phase
    arrows(x0 = 0.4, y0 = 1.6, x1 = 0.75, y1 = 1.6, length = 0.08, lwd = 1.5, col = "darkgoldenrod", code = 2)
    text(x = 0.58, y = 1.65, labels = "Pos IOD", col = "darkgoldenrod", cex = 0.75, font = 2)
    
    # Plot the pooled lag-specific effect estimate line in rich forest green
    lines(p_dmi$predvar, p_dmi$allRRfit, col = "forestgreen", lwd = 2.5)
    
    # Legend placement matching your ONI configuration
    legend("topleft", legend = c(paste("Pooled Average (Lag", l, ")"), "95% Confidence Interval"),
           col = c("forestgreen", NA), lwd = c(2.5, NA), 
           fill = c(NA, rgb(0.13, 0.55, 0.13, 0.1)), border = c(NA, NA),
           bty = "n", cex = 0.85, y.intersp = 1.1)
    
  } else {
    # Gracefully skip and denote lags with missing data fields
    plot.new()
    title(main = paste("DMI Lag", l, "— (Insufficient Data Available)"))
  }
}

# Add overall master title above the column stack
mtext("Lag-Specific Association Between DMI and Leptospirosis", 
      outer = TRUE, cex = 1.2, font = 2)


# ------------------------------------------------------------------------------
# 12. Extract key national ONI and DMI estimates
# ------------------------------------------------------------------------------

# Maximum cumulative ONI-associated relative risk
max_oni_index <- which.max(pred_oni_overall$allRRfit)

max_oni_rr <- pred_oni_overall$allRRfit[max_oni_index]
max_oni_lci <- pred_oni_overall$allRRlow[max_oni_index]
max_oni_uci <- pred_oni_overall$allRRhigh[max_oni_index]
oni_at_max <- pred_oni_overall$predvar[max_oni_index]


# ONI value at which the lower 95% CI first exceeds 1
significant_oni <- which(
  pred_oni_overall$allRRlow > 1 &
    pred_oni_overall$predvar > 0
)

if (length(significant_oni) > 0) {
  
  oni_sig_index <- significant_oni[1]
  
  oni_at_significance <- pred_oni_overall$predvar[oni_sig_index]
  oni_sig_rr <- pred_oni_overall$allRRfit[oni_sig_index]
  oni_sig_lci <- pred_oni_overall$allRRlow[oni_sig_index]
  
} else {
  
  oni_at_significance <- NA_real_
  oni_sig_rr <- NA_real_
  oni_sig_lci <- NA_real_
}


# Minimum cumulative ONI-associated relative risk during negative ONI values
negative_oni <- which(
  pred_oni_overall$predvar < 0
)

if (length(negative_oni) > 0) {
  
  minimum_oni_index <- negative_oni[
    which.min(pred_oni_overall$allRRfit[negative_oni])
  ]
  
  minimum_oni_rr <- pred_oni_overall$allRRfit[minimum_oni_index]
  oni_at_minimum <- pred_oni_overall$predvar[minimum_oni_index]
  
} else {
  
  minimum_oni_rr <- NA_real_
  oni_at_minimum <- NA_real_
}


# Maximum cumulative DMI-associated relative risk
max_dmi_index <- which.max(pred_dmi_overall$allRRfit)

max_dmi_rr <- pred_dmi_overall$allRRfit[max_dmi_index]
max_dmi_lci <- pred_dmi_overall$allRRlow[max_dmi_index]
max_dmi_uci <- pred_dmi_overall$allRRhigh[max_dmi_index]
dmi_at_max <- pred_dmi_overall$predvar[max_dmi_index]


# ------------------------------------------------------------------------------
# 13. Save key model outputs
# ------------------------------------------------------------------------------

national_estimates <- data.frame(
  Index = c("ONI", "DMI"),
  Maximum_RR = c(max_oni_rr, max_dmi_rr),
  Lower_95CI = c(max_oni_lci, max_dmi_lci),
  Upper_95CI = c(max_oni_uci, max_dmi_uci),
  Exposure_at_Maximum = c(oni_at_max, dmi_at_max)
)


print(stage1_summary,n="inf")
print(national_estimates,n="inf")

# ------------------------------------------------------------------------------
# END
# ------------------------------------------------------------------------------


