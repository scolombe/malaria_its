# =============================================================================
# Malaria interrupted time-series analysis — Burundi
# Negative binomial regression assessing the effect of indoor residual
# spraying (IRS) on malaria incidence in Kinyinya
# Sensitivity analysis: simple ITS model, with LLIN as co-variate
# =============================================================================

rm(list = ls(all = TRUE))

# 1. Libraries -----------------------------------------------------------------
# Only what is actually used.

library(MASS)      # glm.nb()
library(sandwich)  # NeweyWest()
library(broom)     # tidy()
library(pander)    # printing tables
library(MMWRweek)  # MMWRweek2Date()
library(dplyr)
library(tidyr)
library(ggplot2)

# 2. Settings ------------------------------------------------------------------

per       <- 1e5          # incidence expressed per 100,000
nw_lag    <- 2            # lag for Newey-West standard errors
z         <- qnorm(0.975) # for 95% confidence intervals

# 3. Helper functions ----------------------------------------------------------

#' Model-based table (exponentiated coefficients, profile-likelihood 95% CIs)
glmtidy <- function(model, caption = "") {
  pander(
    tidy(model, exponentiate = TRUE, conf.int = TRUE),
    caption     = caption,
    split.cell  = 10,
    split.table = Inf
  )
}

#' Add sine/cosine seasonal terms for the given periods (in weeks)
add_harmonics <- function(df, time, periods) {
  for (p in periods) {
    df[[paste0("sin", p)]] <- sin(2 * pi * time / p)
    df[[paste0("cos", p)]] <- cos(2 * pi * time / p)
  }
  df
}

# 4. Import and prepare data ---------------------------------------------------
data_file <- "data/Palu_hebdo_analyse_p_FINAL_Kin_paper.csv"

raw <- read.csv(data_file)

kin <- raw %>%
  mutate(week_start = MMWRweek2Date(MMWRyear = year, MMWRweek = week, MMWRday = 1)) %>%
  arrange(week_start)

# Impute missing Kinyinya case counts as the mean of the neighbouring weeks.
missing_kin <- kin %>% filter(is.na(Kin.DSNIS)) %>% select(year, week)
if (nrow(missing_kin) > 0) {
  message("Imputing Kin.DSNIS for: ",
          paste0(missing_kin$year, "-w", missing_kin$week, collapse = ", "))
  kin$Kin.DSNIS <- round(zoo::na.approx(kin$Kin.DSNIS, na.rm = FALSE))
}

kin <- kin %>%
  mutate(
    # Time index in weeks since the first week, based on real dates.
    time   = as.numeric(week_start - min(week_start)) / 7,
    KinInc = Kin.DSNIS / DenomKin,
    GihInc = Gih.DSNIS / DenomGih,
    logc   = log(GihInc),     # control-site incidence (offset)
    logpop = log(DenomKin)    # Kinyinya population (offset)
  ) %>%
  add_harmonics(time = .$time, periods = c(26, 52))

# Exclude the weeks during which IRS was being implemented (coded NA)
model_data <- kin

# 5. Model: uncontrolled -------------------------------------
# Counts as the outcome; population as offsets.
# This models log(Kinyinya incidence).

fit <- glm.nb(
  Kin.DSNIS ~ outbreak*(time + sin26 + cos26 + sin52 + cos52) +
    IRS2019sep + TimeI1 +
    IRS2020sep + TimeI2 +
    IRS2021sep + TimeI3 +
    CDM14 + CDM17 + CDM19 + CDM22 + TimeCDM14 + TimeCDM17  +  TimeCDM19 + TimeCDM22 +
    offset(logpop),
  data = model_data
)
 

summary(fit)
AIC(fit)
if (!is.null(fit$th.warn)) warning("theta estimation: ", fit$th.warn)

# Model-based table (profile-likelihood CIs)
glmtidy(fit, "Uncontrolled model (model-based CIs)")

# 6. Table with Newey-West (autocorrelation-robust) CIs ------------------------
# Note: excluded IRS weeks create gaps, so the NW lag is approximate there.

vc_nw <- NeweyWest(fit, lag = nw_lag, prewhite = FALSE)
se_nw <- sqrt(diag(vc_nw))
b     <- coef(fit)

table_nw <- tibble(
  Coefficient = names(b),
  IRR         = exp(b),
  Lower_CI    = exp(b - z * se_nw),
  Upper_CI    = exp(b + z * se_nw),
  p_value     = 2 * pnorm(-abs(b / se_nw))
) %>%
  mutate(across(c(IRR, Lower_CI, Upper_CI), ~ round(.x, 3)),
         p_value = signif(p_value, 3))

pander(table_nw, caption = "Uncontrolled model (Newey-West CIs)")


