# =============================================================================
# Malaria interrupted time-series analysis — Ryansoro
# Negative binomial regression assessing the effect of indoor residual
# spraying (IRS) on malaria incidence in Ryansoro, using Kibuye as a control
# site (entered as an offset).
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
data_file <- "xxx"
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

# 5. Model: with control site, without CDM -------------------------------------
# Counts as the outcome; population and control incidence as offsets.
# This models log(Ryansoro incidence / Kibuye incidence).

fit <- glm.nb(
  Rya.DSNIS ~ time +
    IRS2018feb + TimeI1 +
    IRS2018aug + TimeI2 +
    IRS2021jul + TimeI3 +
    sin26 + cos26 + sin52 + cos52 +
    offset(logc + logpop),
  data = model_data
)

summary(fit)
AIC(fit)
if (!is.null(fit$th.warn)) warning("theta estimation: ", fit$th.warn)

# Model-based table 
glmtidy(fit, "Controled model (model-based CIs)")

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

pander(table_nw, caption = "Controled model (Newey-West CIs)")

# 7. Fitted incidence and trend ------------------------------------------------

# Fitted incidence with 95% CI: predict on the link scale (offsets included),
# build the CI there, back-transform, and divide by the population.
pred <- predict(fit, newdata = model_data, type = "link", se.fit = TRUE)

model_data <- model_data %>%
  mutate(
    fit_inc = exp(pred$fit                   - logpop) * per,
    lwr_inc = exp(pred$fit - z * pred$se.fit - logpop) * per,
    upr_inc = exp(pred$fit + z * pred$se.fit - logpop) * per
  )

# Trend: same model with seasonality removed (sin/cos = 0) and the control
# incidence fixed at its mean (on the log scale). Computed from the fitted
# coefficients, so it always matches the current model.
a  <- mean(log(model_data$KibInc), na.rm = TRUE)
cf <- coef(fit)

model_data$trend <- exp(
  cf["(Intercept)"] +
    cf["time"]       * model_data$time +
    cf["IRS2018feb"] * model_data$IRS2018feb +
    cf["IRS2018aug"] * model_data$IRS2018aug +
    cf["IRS2021jul"] * model_data$IRS2021jul +
    cf["TimeI1"]     * model_data$TimeI1 +
    cf["TimeI2"]     * model_data$TimeI2 +
    cf["TimeI3"]     * model_data$TimeI3 +
    a
)
model_data$trend_inc <- model_data$trend * per   # per 100,000 for the plot

# 8. Plot ----------------------------------------------------------------------

# All weeks
plot_data <- rya %>%
  transmute(week_start, obs_inc = RyaInc * per) %>%
  left_join(model_data %>% select(week_start, fit_inc, lwr_inc, upr_inc, trend_inc),
            by = "week_start")

# The IRS effect is assumed to start when implementation ends
periods <- tibble(
  type  = c(rep("IRS implementation", length(irs_start)),
            rep("LLIN distribution", length(llin_start)),
            rep("Estimated duration of IRS effect", length(irs_end))),
  start = c(irs_start, llin_start, irs_end),
  end   = c(irs_end,   llin_end,   effect_end)
)

line_cols <- c("Observed incidence"  = "black",
               "Estimated incidence" = "#4EADED",
               "Estimated trend"     = "#ED684E")

period_fills <- c("IRS implementation"               = alpha("black", 0.8),
                  "LLIN distribution"                = alpha("purple", 0.3),
                  "Estimated duration of IRS effect" = alpha("lightblue", 0.4))

trendplot <- ggplot(plot_data, aes(x = week_start)) +
  geom_rect(data = periods, inherit.aes = FALSE,
            aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf, fill = type)) +
  geom_ribbon(aes(ymin = lwr_inc, ymax = upr_inc), fill = "#4EADED", alpha = 0.25) +
  geom_point(aes(y = obs_inc, colour = "Observed incidence"), size = 0.8) +
  geom_line(aes(y = fit_inc, colour = "Estimated incidence"), linewidth = 0.8) +
  geom_line(aes(y = trend_inc, colour = "Estimated trend"), linetype = 5, linewidth = 1.2) +
  scale_colour_manual(
    values = line_cols, breaks = names(line_cols),
    guide = guide_legend(order = 1,
                         override.aes = list(shape = c(16, NA, NA),
                                             linetype = c(0, 1, 5)))
  ) +
  scale_fill_manual(values = period_fills, breaks = names(period_fills),
                    guide = guide_legend(order = 2)) +
  scale_x_date(date_breaks = "20 weeks", date_labels = "%W-%Y") +
  labs(x = "Reporting week",
       y = "Malaria incidence per 100,000",
       caption = "Shaded blue band: 95% CI of the estimated incidence") +
  theme_bw() +
  theme(
    legend.title       = element_blank(),
    legend.position    = "top",
    legend.direction   = "horizontal",
    legend.box         = "vertical",
    axis.text.x        = element_text(angle = 45, vjust = 1, hjust = 1),
    panel.grid.major.x = element_blank(),
    axis.ticks         = element_line(colour = "grey"),
    text               = element_text(size = 14)
  )

trendplot
# ggsave("output/trendplot_Ryansoro.png", trendplot, width = 12, height = 7, dpi = 300)

# 9. Model diagnostics ---------------------------------------------------------

# Residuals and fitted values, aligned with the rows actually used in the fit
# (glm.nb silently drops rows with missing values).
diag_data <- model_data[as.integer(rownames(model.frame(fit))), ]
diag_data$res_dev  <- residuals(fit, type = "deviance")
diag_data$res_pear <- residuals(fit, type = "pearson")
diag_data$mu       <- fitted(fit)   # expected weekly case counts

res <- diag_data$res_dev

# 9a. Periodogram--------
TSA::periodogram(res)
title(main = "Periodogram of deviance residuals")

# 9b. Distribution of residuals ------------------------------------------------
# Not an assumption of the NB model, but large departures flag outlier weeks.
hist(res, probability = TRUE, breaks = 22,
     main = "Histogram of deviance residuals with normal curve",
     xlab = "Deviance residuals")
curve(dnorm(x, mean = mean(res), sd = sd(res)), add = TRUE, col = "green", lwd = 2)

shapiro.test(res)

# 9c. Mean-variance relationship -----------------------------------------------
# Weeks are grouped by their fitted mean; within each group the observed
# variance, mean((y - mu)^2), is compared with the variance the model assumes:
# Poisson = mu, negative binomial = mu + mu^2 / theta.
n_groups <- 20

mv <- diag_data %>%
  mutate(group = ntile(mu, n_groups)) %>%
  group_by(group) %>%
  summarise(mean_fit = mean(mu),
            var_obs  = mean((Rya.DSNIS - mu)^2),
            .groups  = "drop")

plot(mv$mean_fit, mv$var_obs, log = "xy", pch = 16,
     xlab = "Mean fitted weekly cases (group)",
     ylab = "Observed variance (group)",
     main = "Mean-variance plot")
curve(x + x^2 / fit$theta, add = TRUE, col = "red", lwd = 2)
curve(x * 1, add = TRUE, col = "grey40", lty = 2)
legend("topleft", bty = "n",
       legend = c("Negative binomial: mu + mu^2/theta", "Poisson: mu"),
       col = c("red", "grey40"), lty = c(1, 2), lwd = c(2, 1))

# Pearson dispersion: should be close to 1 if the NB variance is adequate
sum(diag_data$res_pear^2) / df.residual(fit)

# 9d. Residuals over time: leftover trend --------------------------------------
res_trend <- lm(res_dev ~ time, data = diag_data)
summary(res_trend)

plot(diag_data$week_start, diag_data$res_dev, pch = 16, cex = 0.6,
     xlab = "Week", ylab = "Deviance residuals",
     main = "Residuals over time with linear trend")
abline(h = 0, lty = 2, col = "grey")
lines(diag_data$week_start, fitted(res_trend), col = "red", lwd = 2)

# 9e. Autocorrelation ----------------------------------------------------------
# Residuals placed on the full weekly calendar so the lags correspond to real weeks.
res_full <- rya %>%
  select(week_start) %>%
  left_join(diag_data %>% select(week_start, res_dev), by = "week_start") %>%
  pull(res_dev)

pacf(res_full, lag.max = 50, na.action = na.pass,
     main = "Partial autocorrelation of deviance residuals")
legend("topright", legend = "+/- 95% CI", col = "blue", lty = 2, lwd = 1, bty = "n")

# Ljung-Box test 
# Significant autocorrelation supports the use of Newey-West standard errors.
sapply(c(1, 4, 10, 26), function(l)
  Box.test(res, lag = l, type = "Ljung-Box")$p.value)

# 10. Parallel trends: Ryansoro vs Kibuye before the intervention --------------
# The model uses Kibuye's incidence as an offset, which assumes that, without
# IRS, log incidence in both sites would move in parallel. We check this by
# testing for a trend in the difference of log incidence in the pre-period.
# Incidence (not raw counts) is used, consistent with the model offsets:
# with counts, diverging population growth would show up as a false trend.

# Pre-intervention definition: weeks before the first IRS campaign
pre_def <- "TimeI1"

par_data <- model_data %>%
  mutate(
    log_rya = log(RyaInc * per),   # log incidence per 100,000, Ryansoro
    log_kib = log(KibInc * per),   # log incidence per 100,000, Kibuye
    d       = log_rya - log_kib    # log incidence ratio Ryansoro / Kibuye
  ) %>%
  filter(is.finite(d))             # drops weeks with 0 cases (log(0) = -Inf)

pre_data <- par_data %>% filter(.data[[pre_def]] == 0)

# 10a. Trend in the log difference, pre-period only ----------------------------
testdiff <- lm(d ~ time, data = pre_data)
summary(testdiff)

# Newey-West test for the slope (weekly change in the log incidence ratio)
se_trend <- sqrt(diag(NeweyWest(testdiff, prewhite = FALSE, lag = nw_lag)))["time"]
z_trend  <- coef(testdiff)["time"] / se_trend
p_trend  <- 2 * pnorm(-abs(z_trend))
c(slope = unname(coef(testdiff)["time"]),
  se_nw = unname(se_trend),
  p     = unname(p_trend))

# 10b. Log difference pre-intervention -----------------------------------------
plot(pre_data$week_start, pre_data$d, pch = 16, cex = 0.6,
     main = "Difference between log of incidence in Ryansoro and Kibuye, pre-intervention",
     ylab = "Difference in log(incidence)", xlab = "Week")
lines(pre_data$week_start, fitted(testdiff), col = "red", lwd = 2)
legend("topleft", bty = "n", legend = "Linear trend",
       col = "red", lty = 1, lwd = 2)

# 10c. Log incidence in both sites, pre-period ---------------------------------
yl <- range(c(pre_data$log_rya, pre_data$log_kib))
plot(pre_data$week_start, pre_data$log_kib, type = "l", lwd = 2, col = "darkgreen",
     ylim = yl,
     main = "Log incidence in Ryansoro (blue) and Kibuye (green), pre-intervention",
     ylab = "Log incidence per 100,000", xlab = "Week")
lines(pre_data$week_start, pre_data$log_rya, lwd = 2, col = "blue")
legend("topright", bty = "n", legend = c("Ryansoro", "Kibuye"),
       col = c("blue", "darkgreen"), lwd = 2)

# 11. Cases averted by the 2021 IRS campaign -----------------------------------
# Cases averted = counterfactual prediction (no 2021 campaign) minus the model's
# fitted cases, summed over the period during which IRS is assumed effective.
# Bounds: counterfactual +/- 2 SE (response scale), summed over the weeks.
# Everything here is in weekly CASE COUNTS (not incidence per 100,000).

# Window: end of the 2021 implementation to the end of its assumed effect
window_start <- irs_end[3]      # 2021-07-28
window_end   <- effect_end[3]   # 2022-03-28

# Counterfactual: the 2021 campaign never happened
cf_data <- rya %>% mutate(IRS2021jul = 0, TimeI3 = 0)

# Both predictions are made on the full weekly dataset (rya), so the rows line
# up week by week. 
fit_cases <- predict(fit, newdata = rya,     type = "response")
cf_pred   <- predict(fit, newdata = cf_data, type = "response", se.fit = TRUE)

averted <- rya %>%
  transmute(
    week_start,
    fit_cases = as.numeric(fit_cases),
    in_window = week_start >= window_start & week_start <= window_end,
    cf_cases  = if_else(in_window, as.numeric(cf_pred$fit), NA_real_),
    cf_lwr    = cf_cases - 2 * as.numeric(cf_pred$se.fit),
    cf_upr    = cf_cases + 2 * as.numeric(cf_pred$se.fit),
    evit      = cf_cases - fit_cases,
    evit_lwr  = cf_lwr   - fit_cases,
    evit_upr  = cf_upr   - fit_cases
  )

### Nombre de cas évités
cases_averted <- c(
  averted = sum(averted$evit,     na.rm = TRUE),
  lower   = sum(averted$evit_lwr, na.rm = TRUE),
  upper   = sum(averted$evit_upr, na.rm = TRUE)
)
round(cases_averted)

# Plot: fitted cases (black) vs counterfactual without the 2021 campaign (pink)
zoom_x  <- as.Date(c("2014-03-01", "2024-12-31"))   # narrow this for a real close-up
zoom_y  <- c(0, 5000)
arrow_y <- 0.95 * zoom_y[2]

closeup_econom <- ggplot(averted, aes(x = week_start)) +
  annotate("rect", xmin = irs_start[3], xmax = irs_end[3], ymin = -Inf, ymax = Inf,
           colour = "gray", fill = "gray", alpha = 0.3) +
  annotate("rect", xmin = window_start, xmax = window_end, ymin = -Inf, ymax = Inf,
           colour = "lightblue", fill = "lightblue", alpha = 0.2) +
  geom_ribbon(aes(ymin = cf_lwr, ymax = cf_upr), fill = "pink", alpha = 0.5) +
    geom_line(aes(y = fit_cases), colour = "black", linewidth = 0.5) +
  geom_line(aes(y = cf_cases),  colour = "pink",  linewidth = 1) +
  geom_line(aes(y = cf_upr),    colour = "pink") +
  geom_line(aes(y = cf_lwr),    colour = "pink") +
  annotate("segment", x = window_start, xend = window_end, y = arrow_y, yend = arrow_y,
           arrow = arrow(length = unit(0.3, "cm"))) +
  scale_x_date(date_breaks = "10 weeks", date_labels = "%W-%Y") +
  labs(x = "", y = "Number of malaria cases") +
  coord_cartesian(xlim = zoom_x, ylim = zoom_y) +
  theme_bw(base_size = 24) +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    axis.ticks  = element_line(colour = "grey"),
    plot.title  = element_text(size = 16, face = "bold")
  )+
coord_cartesian(xlim=c(as.Date("2021-01-01"), as.Date("2022-12-31")), ylim=c(0,7500))

closeup_econom

# ggsave("output/cases_averted_2021_zoom.png", closeup_econom, width = 10, height = 4.5, dpi = 300)

