###############################################################################
# Duration of the effect of IRS interventions on malaria incidence (Kinshasa)
#
# Method (see manuscript, Methods):
#   1. Fit a generalized additive model (GAM; mgcv) of weekly malaria incidence
#      (cases per 100,000 population) as a
#      function of year, annual and semi-annual harmonic seasonal terms, and a
#      separate smooth of the time index for each IRS intervention.
#   2. Check model fit visually (diagnostic and residual plots).
#   3. Compute the first derivative of each intervention smooth (gratia).
#   4. Define the "plateau" as the range where the 95% CI of the derivative
#      includes zero, i.e. the intervention no longer has an active effect.
#   5. Report the start, midpoint and end of each plateau in months.
#
# Input : data/Palu_hebdo_analyse_p_FINAL_Kin.csv (weekly surveillance data)
# Output: diagnostic plots, derivative plots, smooth plots with plateau,
#         and a summary table of plateau timing (months).
###############################################################################


# ---- 0. Setup ---------------------------------------------------------------

rm(list = ls(all.names = TRUE))

library(MMWRweek)  # MMWR week -> calendar date
library(mgcv)      # GAM fitting
library(gratia)    # derivatives, smooth estimates, diagnostics
library(ggplot2)   # plotting

# Weeks elapsed since 1960 (Stata-style weekly date, 52 weeks per year)
stataweekdate <- function(year, week) {
  (year - 1960) * 52 + week - 1
}

# Convert weeks to months (1 month ~ 4.345 weeks)
weeks_to_months <- function(x) round(x / 4.345, 1)


# ---- 1. Import and prepare data ---------------------------------------------

datapalu <- read.csv("data/Palu_hebdo_analyse_p_FINAL_Kin.csv")

# Impute missing weekly case counts in Kinshasa with the mean of the two
# adjacent weeks (3872 and 3691 cases)
datapalu$Kin.DSNIS[is.na(datapalu$Kin.DSNIS)] <- (3872 + 3691) / 2

# Incidence per 100,000 population
datapalu$KinInc100k <- (datapalu$Kin.DSNIS / datapalu$DenomKin) * 100000

# Weekly time index and calendar date
datapalu$Date     <- with(datapalu, stataweekdate(year, week))
datapalu          <- datapalu[order(datapalu$Date), ]
datapalu$CalDate  <- MMWRweek2Date(MMWRyear = datapalu$year,
                                   MMWRweek = datapalu$week,
                                   MMWRday  = 1)

# Variables used in the analysis
datapalukin <- datapalu[, c("Date", "CalDate", "year", "week",
                            "Kin.DSNIS", "DenomKin", "KinInc100k",
                            "IRS2019sep", "IRS2020sep", "IRS2021sep",
                            "TimeI1", "TimeI2", "TimeI3")]

# Seasonal harmonic terms: annual (52-week) and semi-annual (26-week) cycles
datapalukin$sin52 <- sin(2 * pi * datapalukin$Date / 52)
datapalukin$cos52 <- cos(2 * pi * datapalukin$Date / 52)
datapalukin$sin26 <- sin(2 * pi * datapalukin$Date / 26)
datapalukin$cos26 <- cos(2 * pi * datapalukin$Date / 26)

# Exclude the weeks during which IRS was being implemented
# (coded NA in the IRS indicator variables)
datapalukin_model <- subset(datapalukin,
                            !is.na(IRS2019sep) &
                            !is.na(IRS2020sep) &
                            !is.na(IRS2021sep))


# ---- 2. Fit the GAM ---------------------------------------------------------
# TimeI1, TimeI2, TimeI3 = weeks since each IRS intervention (time indices)

# Response: weekly malaria incidence per 100,000 population
gam_mod_KIN <- gam(
  KinInc100k ~ s(year, k = 1) +
    sin52 + cos52 +            # annual seasonality
    sin26 + cos26 +            # semi-annual seasonality
    s(TimeI1) + s(TimeI2) + s(TimeI3),
  data = datapalukin_model
)

# Alternative (count model with population offset; effects on log-incidence scale):
# gam_mod_KIN <- gam(
#   Kin.DSNIS ~ s(year, k = 1) + sin52 + cos52 + sin26 + cos26 +
#     s(TimeI1) + s(TimeI2) + s(TimeI3) + offset(log(DenomKin / 100000)),
#   family = nb(), data = datapalukin_model
# )

summary(gam_mod_KIN)


# ---- 3. Model diagnostics (visual inspection) -------------------------------

# Residual diagnostics: QQ plot, residuals vs linear predictor,
# histogram of residuals, observed vs fitted
appraise(gam_mod_KIN)

# Basis dimension check and convergence information
gam.check(gam_mod_KIN)

# Partial effect of each smooth with partial residuals
draw(gam_mod_KIN, residuals = TRUE)


# ---- 4. Derivatives and plateau of each intervention effect -----------------

# Function: for one intervention smooth, compute the first derivative,
# identify the plateau (95% CI of the derivative includes zero) and plot it.
analyse_plateau <- function(model, time_var, n = 200) {

  smooth_term <- paste0("s(", time_var, ")")

  # First derivative of the smooth with 95% confidence interval
  deriv_df <- derivatives(model, select = smooth_term, n = n,
                          type = "central", level = 0.95)

  # Plateau: points where the derivative is not significantly different from 0
  plateau_pts <- deriv_df[deriv_df$.lower_ci < 0 & deriv_df$.upper_ci > 0, ]

  if (nrow(plateau_pts) == 0) {
    warning("No plateau found for ", smooth_term,
            ": the derivative is significantly different from 0 everywhere.")
    return(NULL)
  }

  # Warn if the non-significant region is not one contiguous block
  # (e.g. the derivative also crosses zero at a peak or trough of the curve)
  in_plateau <- deriv_df$.lower_ci < 0 & deriv_df$.upper_ci > 0
  if (sum(diff(in_plateau) == 1) + in_plateau[1] > 1) {
    warning("The plateau for ", smooth_term, " is made of several separate ",
            "regions; start/end span all of them. Check the derivative plot.")
  }

  plateau_start <- min(plateau_pts[[time_var]])
  plateau_end   <- max(plateau_pts[[time_var]])
  plateau_mid   <- mean(c(plateau_start, plateau_end))

  # Plot of the derivative
  p_deriv <- draw(deriv_df) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "red") +
    labs(title = paste("First derivative of", smooth_term))

  # Plot of the smooth with the plateau region
  smooth_est <- smooth_estimates(model, select = smooth_term)
  y_top      <- max(smooth_est$.estimate + 1.96 * smooth_est$.se)

  p_smooth <- ggplot(smooth_est, aes(x = .data[[time_var]], y = .estimate)) +
    annotate("rect", xmin = plateau_start, xmax = plateau_end,
             ymin = -Inf, ymax = Inf, fill = "lightblue", alpha = 0.3) +
    geom_ribbon(aes(ymin = .estimate - 1.96 * .se,
                    ymax = .estimate + 1.96 * .se), alpha = 0.2) +
    geom_line() +
    geom_vline(xintercept = c(plateau_start, plateau_end),
               colour = "black", linetype = "dotted") +
    geom_vline(xintercept = plateau_mid, colour = "purple", linetype = "dashed") +
    annotate("text", x = plateau_start, y = y_top, size = 3.5,
             vjust = -0.5, hjust = -0.1,
             label = paste0("Start ≈ ", weeks_to_months(plateau_start), " mo")) +
    annotate("text", x = plateau_end, y = y_top, size = 3.5,
             vjust = -0.5, hjust = 1.1,
             label = paste0("End ≈ ", weeks_to_months(plateau_end), " mo")) +
    annotate("text", x = plateau_mid, y = y_top, size = 3.5,
             vjust = 1.5, hjust = 0.5, colour = "purple",
             label = paste0("Midpoint ≈ ", weeks_to_months(plateau_mid), " mo")) +
    labs(title = paste("Smooth of", time_var, "with plateau region"),
         x = paste(time_var, "(weeks since IRS)"),
         y = "Estimated effect on incidence (per 100,000)") +
    theme_minimal()

  print(p_deriv)
  print(p_smooth)

  data.frame(intervention  = time_var,
             start_weeks   = plateau_start,
             mid_weeks     = plateau_mid,
             end_weeks     = plateau_end,
             start_months  = weeks_to_months(plateau_start),
             mid_months    = weeks_to_months(plateau_mid),
             end_months    = weeks_to_months(plateau_end))
}

# Run for each IRS intervention
plateau_results <- do.call(rbind, lapply(c("TimeI1", "TimeI2", "TimeI3"),
                                         function(v) analyse_plateau(gam_mod_KIN, v)))

plateau_results
