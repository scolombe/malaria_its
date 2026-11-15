# =============================================================================
# Malaria interrupted time-series analysis — Burundi
# Negative binomial regression assessing the effect of indoor residual
# spraying (IRS) on malaria incidence in Ryansoro
#Sensitivity analysis: no control time series.
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

# Key dates. 
irs_start    <- as.Date(c("2018-02-26", "2018-08-20", "2021-06-14"))
irs_end      <- as.Date(c("2018-03-12", "2018-09-03", "2021-07-28"))
effect_end   <- as.Date(c("2018-09-12", "2019-03-03", "2022-03-28"))
llin_start   <- as.Date(c("2014-06-01", "2017-06-01", "2019-12-01", "2022-09-01"))
llin_end     <- as.Date(c("2014-06-30", "2017-06-30", "2019-12-31", "2022-09-30"))

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
data_file <- "data/Palu_hebdo_analyse_p_FINAL_Ryansoro.csv"
raw <- read.csv(data_file)

rya <- raw %>%
  mutate(week_start = MMWRweek2Date(MMWRyear = year, MMWRweek = week, MMWRday = 1)) %>%
  arrange(week_start)

# Impute missing Ryansoro case counts as the mean of the neighbouring weeks.
missing_rya <- rya %>% filter(is.na(Rya.DSNIS)) %>% select(year, week)
if (nrow(missing_rya) > 0) {
  message("Imputing Rya.DSNIS for: ",
          paste0(missing_rya$year, "-w", missing_rya$week, collapse = ", "))
  rya$Rya.DSNIS <- round(zoo::na.approx(rya$Rya.DSNIS, na.rm = FALSE))
}

model_data <- rya %>%
  mutate(
    # Time index in weeks since the first week, based on real dates.
    time   = as.numeric(week_start - min(week_start)) / 7,
    RyaInc = Rya.DSNIS / DenomRya,
    KibInc = Kib.DSNIS / DenomKib,
    logc   = log(KibInc),     # control-site incidence (offset)
    logpop = log(DenomRya)    # Ryansoro population (offset)
  ) %>%
  add_harmonics(time = .$time, periods = c(26, 52))

# 5. Model: uncontrolled, with LLIN as variable -------------------------------------
# Counts as the outcome; population as offsets.
# This models log(Ryansoro incidence).

fit <- glm.nb(
  Rya.DSNIS ~ time + CDM14 + CDM17 + CDM19 + CDM22 + TimeCDM14 + TimeCDM17 + TimeCDM19 + TimeCDM22+
    IRS2018feb + TimeI1 +
    IRS2018aug + TimeI2 +
    IRS2021jul + TimeI3 +
    sin26 + cos26 + sin52 + cos52 +
    offset(logpop),
  data = model_data
)

summary(fit)
AIC(fit)
if (!is.null(fit$th.warn)) warning("theta estimation: ", fit$th.warn)

# Model-based table (profile-likelihood CIs)
glmtidy(fit, "Simple ITS model (model-based CIs)")

# 6. Table with Newey-West (autocorrelation-robust) CIs ------------------------

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

pander(table_nw, caption = "Simple ITS model (Newey-West CIs)")
