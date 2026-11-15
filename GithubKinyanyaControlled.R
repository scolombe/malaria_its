# =============================================================================
# Malaria interrupted time-series analysis — Burundi
# Negative binomial regression assessing the effect of indoor residual
# spraying (IRS) on malaria incidence in Kinyinya, using gihofi as a control
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

# 5. Model: with control site, without CDM -------------------------------------
# Counts as the outcome; population and control incidence as offsets.
# This models log(Kinyinya incidence / Gihofi incidence).

fit <- glm.nb(
  Kin.DSNIS ~ time +
    IRS2019sep + TimeI1 +
    IRS2020sep + TimeI2 +
    IRS2021sep + TimeI3 +
    sin26 + cos26 + sin52 + cos52 +
    offset(logc + logpop),
  data = model_data
)

summary(fit)
AIC(fit)
if (!is.null(fit$th.warn)) warning("theta estimation: ", fit$th.warn)

# Model-based table (profile-likelihood CIs)
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
    fit_inc = exp(pred$fit                  - logpop) * per,
    lwr_inc = exp(pred$fit - z * pred$se.fit - logpop) * per,
    upr_inc = exp(pred$fit + z * pred$se.fit - logpop) * per
  )

# Trend: same model with seasonality removed (sin/cos = 0) and the control
# incidence fixed at its mean (on the log scale). Computed from the fitted
# coefficients, so it always matches the current model.
a  <- mean(log(model_data$GihInc), na.rm = TRUE)
cf <- coef(fit)

model_data$trend <- exp(
  cf["(Intercept)"] +
    cf["time"]       * model_data$time +
    cf["IRS2019sep"] * model_data$IRS2019sep +
    cf["IRS2020sep"] * model_data$IRS2020sep +
    cf["IRS2021sep"] * model_data$IRS2021sep +
    cf["TimeI1"]     * model_data$TimeI1 +
    cf["TimeI2"]     * model_data$TimeI2 +
    cf["TimeI3"]     * model_data$TimeI3 +
    a
)
model_data$trend_inc <- model_data$trend * per   # per 100,000 for the plot

# 8. Plot ----------------------------------------------------------------------

# All weeks
plot_data <- kin %>%
  transmute(week_start, obs_inc = KinInc * per) %>%
  left_join(model_data %>% select(week_start, fit_inc, lwr_inc, upr_inc, trend_inc),
            by = "week_start")

PID2020min <- MMWRweek2Date(MMWRyear = 2020, MMWRweek = 40, MMWRday = 1)
PID2020max <- MMWRweek2Date(MMWRyear = 2020, MMWRweek = 42, MMWRday = 1)

periods <- tibble(
  type = c(rep("IRS implementation", 3),
           rep("LLIN distribution", 4),
           rep("Estimated duration of IRS effect", 3)),
  start = as.Date(c("2019-09-02", as.character(PID2020min), "2021-08-16",
                    "2014-06-01", "2017-06-01", "2019-12-01", "2022-09-01",
                    "2019-10-09", "2020-10-11", "2021-10-08")),
  end   = as.Date(c("2019-10-09", as.character(PID2020max), "2021-10-08",
                    "2014-06-30", "2017-06-30", "2019-12-31", "2022-09-30",
                    "2020-05-09", "2021-05-11", "2022-05-08"))
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
# ggsave("output/trendplot_kinyinya.png", trendplot, width = 12, height = 7, dpi = 300)


# 9. Model diagnostics ---------------------------------------------------------

# Residuals and fitted values, aligned with the rows actually used in the fit
# (glm.nb silently drops rows with missing values).
diag_data <- model_data[as.integer(rownames(model.frame(fit))), ]
diag_data$res_dev  <- residuals(fit, type = "deviance")
diag_data$res_pear <- residuals(fit, type = "pearson")
diag_data$mu       <- fitted(fit)   # expected weekly case counts

res <- diag_data$res_dev

# 9a. Periodogram ------
TSA::periodogram(res)
title(main = "Periodogram of deviance residuals")

# 9b. Distribution of residuals ----------------------------------------------
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
            var_obs  = mean((Kin.DSNIS - mu)^2),
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

# 9d. Residuals over time: leftover trend -------------------------------------
res_trend <- lm(res_dev ~ time, data = diag_data)
summary(res_trend)

plot(diag_data$week_start, diag_data$res_dev, pch = 16, cex = 0.6,
     xlab = "Week", ylab = "Deviance residuals",
     main = "Residuals over time with linear trend")
abline(h = 0, lty = 2, col = "grey")
lines(diag_data$week_start, fitted(res_trend), col = "red", lwd = 2)

# 9e. Autocorrelation ---------------------------------------------------------
# Residuals placed on the full weekly calendar so the lags correspond to real weeks.
res_full <- kin %>%
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


# 10. Parallel trends: Kinyinya vs Gihofi before the intervention --------------
# The model uses Gihofi's incidence as an offset, which assumes that, without
# IRS, log incidence in both sites would move in parallel. We check this by
# testing for a trend in the difference of log incidence in the pre-period.
# Incidence (not raw counts) is used, consistent with the model offsets:
# with counts, diverging population growth would show up as a false trend.

# Pre-intervention definition.

pre_def <- "TimeI1"

par_data <- model_data %>%
  mutate(
    logk = log(KinInc * per),   # log incidence per 100,000, Kinyinya
    logg = log(GihInc * per),   # log incidence per 100,000, Gihofi
    d    = logk - logg          # log incidence ratio Kinyinya / Gihofi
  ) %>%
  filter(is.finite(d))          # drops weeks with 0 cases (log(0) = -Inf)

pre_data <- par_data[par_data[[pre_def]] == 0, ]

# 10a. Trend in the log difference, pre-period only ---------------------------
testdiff <- lm(d ~ time, data = pre_data)
summary(testdiff)

# Newey-West test for the slope (weekly change in the log incidence ratio)
se_trend <- sqrt(diag(NeweyWest(testdiff, prewhite = FALSE, lag = nw_lag)))["time"]
z_trend  <- coef(testdiff)["time"] / se_trend
p_trend  <- 2 * pnorm(-abs(z_trend))
c(slope = unname(coef(testdiff)["time"]),
  se_nw = unname(se_trend),
  p     = unname(p_trend))

# 10b. Log difference pre-intervention ------------------------------------
plot(pre_data$week_start, pre_data$d, pch = 16, cex = 0.6,
     main = "Difference between log incidence in Kinyinya and Gihofi, pre-intervention",
     ylab = "Difference in log incidence", xlab = "Week")
lines(pre_data$week_start, fitted(testdiff), col = "red", lwd = 2)
legend("topleft", bty = "n", legend = "Linear trend",
       col = "red", lty = 1, lwd = 2)

# 10c. Log incidence in both sites, pre-period ---------------------------------
yl <- range(c(pre_data$logk, pre_data$logg))
plot(pre_data$week_start, pre_data$logg, type = "l", lwd = 2, col = "darkgreen",
     ylim = yl,
     main = "Log incidence in Kinyinya (blue) and Gihofi (green), pre-intervention",
     ylab = "Log incidence per 100,000", xlab = "Week")
lines(pre_data$week_start, pre_data$logk, lwd = 2, col = "blue")
legend("topright", bty = "n", legend = c("Kinyinya", "Gihofi"),
       col = c("blue", "darkgreen"), lwd = 2)

# 11. Checking impact of Gihofi's IRS ------------------------------------------

# Pre-intervention definition
pre_def <- "IRSGov"

par_data <- model_data %>%
  mutate(
    logk = log(KinInc * per),   # log incidence per 100,000, Kinyinya
    logg = log(GihInc * per),   # log incidence per 100,000, Gihofi
    d    = logk - logg          # log incidence ratio Kinyinya / Gihofi
  ) %>%
  filter(is.finite(d))          # drops weeks with 0 cases (log(0) = -Inf)

pre_data <- par_data %>% filter(.data[[pre_def]] == 0)

# 11a. Trend in the log difference, pre-period only ---------------------------
testdiff <- lm(d ~ time, data = pre_data)
summary(testdiff)

# Newey-West test for the slope (weekly change in the log incidence ratio)
se_trend <- sqrt(diag(NeweyWest(testdiff, prewhite = FALSE, lag = nw_lag)))["time"]
z_trend  <- coef(testdiff)["time"] / se_trend
p_trend  <- 2 * pnorm(-abs(z_trend))
c(slope = unname(coef(testdiff)["time"]),
  se_nw = unname(se_trend),
  p     = unname(p_trend))

# 11b. Log difference pre-intervention ------------------------------------
plot(pre_data$week_start, pre_data$d, pch = 16, cex = 0.6,
     main = "Difference between log incidence in Kinyinya and Gihofi, pre-intervention",
     ylab = "Difference in log incidence", xlab = "Week")
lines(pre_data$week_start, fitted(testdiff), col = "red", lwd = 2)
legend("topleft", bty = "n", legend = "Linear trend",
       col = "red", lty = 1, lwd = 2)

# 11c. Log incidence in both sites, pre-period ---------------------------------
yl <- range(c(pre_data$logk, pre_data$logg))
plot(pre_data$week_start, pre_data$logg, type = "l", lwd = 2, col = "darkgreen",
     ylim = yl,
     main = "Log incidence in Kinyinya (blue) and Gihofi (green), pre-intervention",
     ylab = "Log incidence per 100,000", xlab = "Week")
lines(pre_data$week_start, pre_data$logk, lwd = 2, col = "blue")
legend("topright", bty = "n", legend = c("Kinyinya", "Gihofi"),
       col = c("blue", "darkgreen"), lwd = 2)


# 12. Cases averted by the 2021 IRS campaign -----------------------------------
# Cases averted = counterfactual prediction (no 2021 campaign) minus the model's
# fitted cases, summed over the period during which IRS is assumed effective.
# Bounds: counterfactual +/- 2 SE (response scale), summed over the weeks.

irs21     <- as.Date(c("2021-08-16", "2021-10-08"))  # IRS implementation 2021
win_start <- as.Date("2021-10-08")                   # start of effect period
win_end   <- as.Date("2022-06-08")                   # end of effect period

# Factual: fitted cases, interventions as they happened
model_data$fitted_cases <- predict(fit, newdata = model_data, type = "response")

# Counterfactual: no 2021 campaign
datanew2021 <- model_data
datanew2021$IRS2021sep <- 0
datanew2021$TimeI3     <- 0

make_prediction <- function(model, newdata) {
  pred <- predict(model, newdata, type = "response", se.fit = TRUE)
  data.frame(
    week_start = newdata$week_start,
    fit        = pred$fit,
    lwr        = pred$fit - 2 * pred$se.fit,
    upr        = pred$fit + 2 * pred$se.fit
  )
}

pred2021 <- make_prediction(fit, datanew2021)

# 12a. Number of cases averted -------------------------------------------------
# Factual and counterfactual are matched by week, and only the effect period
# is summed.
CE <- pred2021 %>%
  left_join(model_data %>% select(week_start, fitted_cases), by = "week_start") %>%
  filter(week_start >= win_start, week_start <= win_end) %>%
  mutate(
    averted     = fit - fitted_cases,
    averted_lwr = lwr - fitted_cases,
    averted_upr = upr - fitted_cases
  )

cases_averted     <- sum(CE$averted,     na.rm = TRUE)
cases_averted_lwr <- sum(CE$averted_lwr, na.rm = TRUE)
cases_averted_upr <- sum(CE$averted_upr, na.rm = TRUE)

c(cases_averted = cases_averted, lower = cases_averted_lwr, upper = cases_averted_upr)

# 12b. Plots: main plot and close-up ------------------------------------------

# Counterfactual kept only in the effect period (set to NA outside it)
pred2021_win <- pred2021
outside_win  <- pred2021_win$week_start < win_start | pred2021_win$week_start > win_end
pred2021_win[outside_win, c("fit", "lwr", "upr")] <- NA

# All weeks (lines break at excluded weeks)
cases_plot_data <- kin %>%
  select(week_start) %>%
  left_join(model_data %>% select(week_start, fitted_cases), by = "week_start") %>%
  left_join(pred2021_win, by = "week_start")

# For the close-up: continuous black line. 
cases_plot_data <- cases_plot_data %>%
  left_join(kin %>% select(week_start, Kin.DSNIS), by = "week_start") %>%
  mutate(cases_line = ifelse(is.na(fitted_cases), Kin.DSNIS, fitted_cases))

# Height of the arrow marking the effect period, just above the data in it
in_win  <- with(cases_plot_data, week_start >= win_start & week_start <= win_end)
y_arrow <- 1.02 * max(c(cases_plot_data$fitted_cases[in_win],
                        cases_plot_data$upr[in_win]), na.rm = TRUE)

# Layers and styling shared by both plots
averted_layers <- list(
  annotate("rect", xmin = irs21[1], xmax = irs21[2], ymin = -Inf, ymax = Inf,
           fill = "gray", colour = "gray", alpha = 0.3),
  annotate("rect", xmin = win_start, xmax = win_end, ymin = -Inf, ymax = Inf,
           fill = "lightblue", colour = "lightblue", alpha = 0.2),
  geom_ribbon(aes(ymin = lwr, ymax = upr), fill = "pink", alpha = 0.5),
  geom_line(aes(y = fitted_cases, colour = "Estimated cases"), linewidth = 0.5),
  geom_line(aes(y = fit, colour = "Expected cases without IRS"), linewidth = 1),
  geom_line(aes(y = upr, colour = "Expected cases without IRS")),
  geom_line(aes(y = lwr, colour = "Expected cases without IRS")),
  annotate("segment", x = win_start, xend = win_end, y = y_arrow, yend = y_arrow,
           arrow = arrow(length = unit(0.3, "cm"))),
  scale_colour_manual(values = c("Estimated cases"            = "black",
                                 "Expected cases without IRS" = "pink")),
  labs(x = "Reporting week", y = "Number of malaria cases"),
  theme_bw(),
  theme(legend.title    = element_blank(),
        legend.position = "top",
        axis.text.x     = element_text(angle = 45, vjust = 1, hjust = 1),
        axis.ticks      = element_line(colour = "grey"),
        text            = element_text(size = 12),
        plot.title      = element_text(size = 16, face = "bold"))
)

# Main plot: whole study period
cas_evites_2021 <- ggplot(cases_plot_data, aes(x = week_start)) +
  averted_layers +
  scale_x_date(date_breaks = "2 months", date_labels = "%b-%Y") +
  ggtitle(paste0("Estimated number of malaria cases averted in Kinyinya, ",
                 format(min(kin$week_start), "%Y"), "-",
                 format(max(kin$week_start), "%Y")))

cas_evites_2021

# Close-up: around the 2021 campaign, cases only (no title, no legend)
y_max_close <- max(c(cases_plot_data$fitted_cases[in_win],
                     cases_plot_data$upr[in_win]), na.rm = TRUE)

closeup_econom <- ggplot(cases_plot_data, aes(x = week_start)) +
  annotate("rect", xmin = irs21[1], xmax = irs21[2], ymin = -Inf, ymax = Inf,
           fill = "gray", colour = "gray", alpha = 0.3) +
  annotate("rect", xmin = win_start, xmax = win_end, ymin = -Inf, ymax = Inf,
           fill = "lightblue", colour = "lightblue", alpha = 0.2) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), fill = "pink", alpha = 0.5) +
  geom_line(aes(y = cases_line), colour = "black", linewidth = 0.5) +
  geom_line(aes(y = fit), colour = "pink", linewidth = 1) +
  geom_line(aes(y = upr), colour = "pink") +
  geom_line(aes(y = lwr), colour = "pink") +
  annotate("segment", x = win_start, xend = win_end,
           y = 1.15 * y_max_close, yend = 1.15 * y_max_close,
           arrow = arrow(length = unit(0.3, "cm"))) +
  scale_x_date(date_breaks = "10 weeks", date_labels = "%W-%Y") +
  coord_cartesian(xlim = as.Date(c("2021-03-01", "2022-12-31")),
                  ylim = c(0, 1.35 * y_max_close)) +
  labs(x = NULL, y = "Number of malaria cases") +
  theme_bw(base_size = 14) +
  theme(legend.position = "none",
        axis.text.x     = element_text(angle = 45, vjust = 1, hjust = 1),
        axis.ticks      = element_line(colour = "grey"))

closeup_econom
# Wide format:
# ggsave("output/cases_averted_2021_zoom.png", closeup_econom, width = 10, height = 4.5, dpi = 300)
