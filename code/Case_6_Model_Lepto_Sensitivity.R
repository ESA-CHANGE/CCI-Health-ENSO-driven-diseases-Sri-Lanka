# ==============================================================================
# Sensitivity analysis: specification of the distributed lag models
# ==============================================================================
#
# Purpose:
# This analysis evaluates the robustness of the primary two-stage DLNM results
# to alternative choices for:
#   1. degrees of freedom for the long-term temporal trend;
#   2. exposure-response knot placement for ONI;
#   3. exposure-response knot placement for DMI;
#   4. maximum lag window; and
#   5. placement of knots in the lag-response function.
#
# The default specification corresponds to the model used in the primary
# analysis. Default knot locations were selected a priori based on inspection
# of the observed time series, exposure distributions, frequency histograms,
# and correlation structure of the climate indices.
#
# All sensitivity models retain the same analytical framework as the primary
# analysis: district-specific first-stage quasipoisson regression models,
# adjustment for month and population, and second-stage multivariate
# random-effects meta-analysis using REML.
#
# ==============================================================================


# ------------------------------------------------------------------------------
# 0. Packages
# ------------------------------------------------------------------------------

library(dplyr)
library(dlnm)
library(mvmeta)
library(splines)


# ------------------------------------------------------------------------------
# 1. Sensitivity-analysis specifications
# ------------------------------------------------------------------------------

# Maximum lag windows evaluated.
#
# The primary analysis uses a 0-6 month lag window. Shorter lag windows are
# evaluated to determine whether the estimated associations depend on including
# longer delayed effects.

lag_scenarios <- list(
  Lag2 = c(0, 2),
  Lag3 = c(0, 3),
  Lag4 = c(0, 4),
  Lag5 = c(0, 5),
  Lag6 = c(0, 6)
)


# Degrees of freedom for the long-term temporal trend.
#
# The primary specification uses 15 degrees of freedom, corresponding to the
# temporal adjustment used in the main analysis. Lower and higher values are
# evaluated to assess sensitivity to the degree of temporal control.

df_grid <- c(5, 10, 15, 19)


# ONI exposure-response knot specifications.
#
# Default:
# Knots at -0.5 and 0.5, selected from inspection of the observed ONI
# time series and exposure distribution and used in the primary analysis.
#
# Wider:
# A broader central exposure range is represented by knots at -0.75 and 0.75.
#
# Quantile:
# Knots are placed at the 33rd and 67th percentiles of the observed ONI
# distribution.

oni_grid <- list(
  Default = c(-0.5, 0.5),
  Wider = c(-0.75, 0.75),
  Quantile = quantile(
    analysis_df$ONI,
    probs = c(0.33, 0.67),
    na.rm = TRUE
  )
)


# DMI exposure-response knot specifications.
#
# Default:
# Knots at -0.10 and 0.268, selected from inspection of the observed DMI
# time series and exposure distribution and used in the primary analysis.
#
# Symmetric:
# A symmetric alternative specification using knots at -0.40 and 0.40.
#
# Quantile:
# Knots are placed at the 33rd and 67th percentiles of the observed DMI
# distribution.

dmi_grid <- list(
  Default = c(-0.10, 0.268),
  Symmetric = c(-0.40, 0.40),
  Quantile = quantile(
    analysis_df$DMI,
    probs = c(0.33, 0.67),
    na.rm = TRUE
  )
)


# Alternative specifications for the placement of knots in the lag-response
# function.
#
# Standard:
# Logarithmically spaced knots, as used in the primary analysis.
#
# Early:
# Knots concentrated earlier in the lag period.
#
# Late:
# Knots shifted towards the later part of the lag period.
#
# Knot locations are subsequently checked to ensure that they fall strictly
# within the evaluated lag range.

lag_knot_grid <- list(
  Standard = function(maxlag) {
    logknots(maxlag, nk = 2)
  },
  Early = function(maxlag) {
    c(
      max(1, floor(maxlag / 3)),
      max(2, floor(2 * maxlag / 3))
    )
  },
  Late = function(maxlag) {
    c(
      max(1, ceiling(maxlag / 2)),
      maxlag - 1
    )
  }
)


# Boundary knots for the exposure-response functions.
#
# These boundaries are held constant across all sensitivity analyses so that
# the sensitivity assessment isolates the effect of internal knot placement.

oni_boundary <- c(-1.6, 2.8)
dmi_boundary <- c(-0.758, 0.964)


# ------------------------------------------------------------------------------
# 2. Data preparation
# ------------------------------------------------------------------------------

districts <- unique(analysis_df$district)
n_dist <- length(districts)

final_analysis_data <- analysis_df %>%
  arrange(district, year, month) %>%
  filter(
    !is.na(ONI),
    !is.na(DMI)
  )


# ------------------------------------------------------------------------------
# 3. QAIC function
# ------------------------------------------------------------------------------

# QAIC is calculated from the Poisson log-likelihood evaluated at the fitted
# values, with the quasipoisson dispersion parameter used to account for
# overdispersion.

fqaic <- function(model) {
  loglik <- sum(
    dpois(
      model$y,
      model$fitted.values,
      log = TRUE
    )
  )
  
  dispersion <- summary(model)$dispersion
  k <- model$rank
  
  -2 * loglik + 2 * k * dispersion
}


# ------------------------------------------------------------------------------
# 4. Storage for sensitivity-analysis results
# ------------------------------------------------------------------------------

results <- list()
counter <- 1


# ------------------------------------------------------------------------------
# 5. Complete sensitivity grid
# ------------------------------------------------------------------------------

for (df_now in df_grid) {
  
  for (oni_name in names(oni_grid)) {
    
    for (dmi_name in names(dmi_grid)) {
      
      for (lag_name in names(lag_scenarios)) {
        
        current_lag <- lag_scenarios[[lag_name]]
        maxlag <- current_lag[2]
        
        for (lag_knot_name in names(lag_knot_grid)) {
          
          current_lag_knots <- lag_knot_grid[[lag_knot_name]](maxlag)
          
          current_lag_knots <- unique(current_lag_knots)
          
          current_lag_knots <- current_lag_knots[
            current_lag_knots > 0 &
              current_lag_knots < maxlag
          ]
          
          if (length(current_lag_knots) == 0) {
            current_lag_knots <- NULL
          }
          
          cat(
            "\n====================================================\n",
            "Temporal df:", df_now,
            "| ONI knots:", oni_name,
            "| DMI knots:", dmi_name,
            "| Lag:", lag_name,
            "| Lag knots:", lag_knot_name,
            "\n"
          )
          
          
          # --------------------------------------------------------------------
          # Storage for district-specific first-stage estimates
          # --------------------------------------------------------------------
          
          coef_oni <- matrix(
            NA_real_,
            nrow = n_dist,
            ncol = 3
          )
          
          coef_dmi <- matrix(
            NA_real_,
            nrow = n_dist,
            ncol = 3
          )
          
          vcov_oni <- vector(
            "list",
            n_dist
          )
          
          vcov_dmi <- vector(
            "list",
            n_dist
          )
          
          qaic_vec <- rep(
            NA_real_,
            n_dist
          )
          
          dev_vec <- rep(
            NA_real_,
            n_dist
          )
          
          convergence <- rep(
            TRUE,
            n_dist
          )
          
          
          # ====================================================================
          # 5.1 First-stage district-specific models
          # ====================================================================
          
          for (i in seq_along(districts)) {
            
            sub <- final_analysis_data %>%
              filter(
                district == districts[i]
              )
            
            if (nrow(sub) < 10) {
              next
            }
            
            tryCatch({
              
              # ONI cross-basis.
              cb_oni <- crossbasis(
                sub$ONI,
                lag = current_lag,
                argvar = list(
                  fun = "ns",
                  knots = oni_grid[[oni_name]],
                  Boundary.knots = oni_boundary
                ),
                arglag = list(
                  fun = "ns",
                  knots = current_lag_knots
                )
              )
              
              
              # DMI cross-basis.
              cb_dmi <- crossbasis(
                sub$DMI,
                lag = current_lag,
                argvar = list(
                  fun = "ns",
                  knots = dmi_grid[[dmi_name]],
                  Boundary.knots = dmi_boundary
                ),
                arglag = list(
                  fun = "ns",
                  knots = current_lag_knots
                )
              )
              
              
              # District-specific regression model.
              #
              # The model simultaneously estimates the associations of ONI
              # and DMI while controlling for the long-term temporal trend,
              # calendar month, and population size.
              
              model <- glm(
                Leptospirosis ~
                  cb_oni +
                  cb_dmi +
                  ns(time, df = df_now) +
                  factor(month) +
                  offset(log(Population)),
                family = quasipoisson(),
                data = sub
              )
              
              
              if (!model$converged) {
                convergence[i] <- FALSE
              }
              
              
              qaic_vec[i] <- fqaic(model)
              dev_vec[i] <- model$deviance
              
              
              # Reduce the cross-basis functions to obtain the overall
              # exposure-response association across the complete lag period,
              # centred at an index value of zero.
              
              cr_oni <- crossreduce(
                cb_oni,
                model,
                type = "overall",
                cen = 0
              )
              
              cr_dmi <- crossreduce(
                cb_dmi,
                model,
                type = "overall",
                cen = 0
              )
              
              
              # Store district-specific coefficients and covariance matrices
              # for second-stage meta-analysis.
              
              coef_oni[i, ] <- coef(cr_oni)
              vcov_oni[[i]] <- vcov(cr_oni)
              
              coef_dmi[i, ] <- coef(cr_dmi)
              vcov_dmi[[i]] <- vcov(cr_dmi)
              
            }, error = function(e) {
              
              convergence[i] <<- FALSE
              
            })
          }
          
          
          # ====================================================================
          # 5.2 Second-stage multivariate meta-analysis
          # ====================================================================
          
          keep_oni <- !is.na(coef_oni[, 1])
          keep_dmi <- !is.na(coef_dmi[, 1])
          
          if (
            sum(keep_oni) < 2 ||
            sum(keep_dmi) < 2
          ) {
            
            cat(
              "Skipped: insufficient districts for meta-analysis.\n"
            )
            
            next
          }
          
          
          # Pool the district-specific ONI exposure-response coefficients
          # using a multivariate random-effects meta-analysis.
          
          meta_oni <- mvmeta(
            coef_oni[keep_oni, ] ~ 1,
            S = vcov_oni[keep_oni],
            method = "reml"
          )
          
          
          # Pool the district-specific DMI exposure-response coefficients
          # using the same second-stage model.
          
          meta_dmi <- mvmeta(
            coef_dmi[keep_dmi, ] ~ 1,
            S = vcov_dmi[keep_dmi],
            method = "reml"
          )
          
          
          # ====================================================================
          # 5.3 Store model-comparison statistics
          # ====================================================================
          
          results[[counter]] <- data.frame(
            Time_df = df_now,
            ONI_knots = oni_name,
            DMI_knots = dmi_name,
            Lag_window = lag_name,
            Lag_knots = lag_knot_name,
            Stage1_QAIC = sum(
              qaic_vec,
              na.rm = TRUE
            ),
            Mean_residual_deviance = mean(
              dev_vec,
              na.rm = TRUE
            ),
            Non_converged = sum(
              !convergence
            ),
            ONI_LogLik = as.numeric(
              logLik(meta_oni)
            ),
            ONI_AIC = AIC(
              meta_oni
            ),
            ONI_BIC = BIC(
              meta_oni
            ),
            DMI_LogLik = as.numeric(
              logLik(meta_dmi)
            ),
            DMI_AIC = AIC(
              meta_dmi
            ),
            DMI_BIC = BIC(
              meta_dmi
            )
          )
          
          counter <- counter + 1
        }
      }
    }
  }
}


# ------------------------------------------------------------------------------
# 6. Assemble and rank sensitivity-analysis results
# ------------------------------------------------------------------------------

sensitivity_matrix <- bind_rows(
  results
)

sensitivity_matrix <- sensitivity_matrix %>%
  arrange(
    Stage1_QAIC,
    ONI_AIC
  )


print(sensitivity_matrix)


# ------------------------------------------------------------------------------
# 7. Identify the five best-fitting fully converged specifications
# ------------------------------------------------------------------------------

# Models are first restricted to specifications in which all contributing
# district-specific models converged. Among these models, specifications are
# ranked using the summed first-stage QAIC, mean residual deviance, and the
# second-stage ONI and DMI AIC values.

best5_models <- sensitivity_matrix %>%
  filter(
    Non_converged == 0
  ) %>%
  arrange(
    Stage1_QAIC,
    Mean_residual_deviance,
    ONI_AIC,
    DMI_AIC
  ) %>%
  slice_head(
    n = 5
  )

print(best5_models)