library(survival)
library(tidyverse)

if (system("whoami", intern = TRUE) == "unicph\\hpl802") {
  path <- "~/Github/GSD-GPC/R/"
} else if (system("whoami", intern = TRUE) == "wwu\\danzerm") {
  path <- "~/../Desktop/projects/GSD-GPC/R/"
} else if (system("whoami", intern = TRUE) == "hpl802") {
  path <- "."
}

sapply(list.files(path, full.names = TRUE, pattern = "*.R"), source)

# Set sample size for two groups
n_E <- 100
n_C <- 100

# Set rates for exponential distributions in both groups
rate_E <- 1
rate_C <- 1

# Set design parameters for rNTB
eps <- 0.5 # time horizon
tau <- 0.0 # equivalence boundary

# Trial schedule
a <- 1.5 # accrual duratio
f <- 0.5 # follow-up duration
t1 <- 1 # first analysis
t2 <- a + f # final analysis

# Number of simulations
runs <- 10000

# Save Kaplan-Meier integrals for both analyses
# Kaplan-Meier integrals compute the probabilistic index
# "EC" means \int S_E dS_C
# "CE" means \int S_E dS_C
km_integral_values_EC_interim <- rep(NA, runs)
km_integral_values_EC_final <- rep(NA, runs)
km_integral_values_CE_interim <- rep(NA, runs)
km_integral_values_CE_final <- rep(NA, runs)

# Save variance estimates for the Kaplan Meier integrals
km_integral_var_EC_interim <- rep(NA, runs)
km_integral_var_EC_final <- rep(NA, runs)
km_integral_var_CE_interim <- rep(NA, runs)
km_integral_var_CE_final <- rep(NA, runs)

# Save variance estimates for rNTBs
rntb_var_interim <- rep(NA, runs)
rntb_var_final <- rep(NA, runs)
rntb_cov <- rep(NA, runs)

for (run in seq(runs)) {
  cat("\r", run)

  # Simulate survival times
  t_E <- rexp(n_E, rate = rate_E)
  t_C <- rexp(n_C, rate = rate_C)

  # Simulate censoring times
  rec_E <- runif(n_E, min = 0, max = a)
  rec_C <- runif(n_C, min = 0, max = a)

  # First analysis data
  interim_df <- data.frame(
    time = c(pmin(t_C, pmax(t1 - rec_C, 0)), pmin(t_E, pmax(t1 - rec_E, 0))),
    status = c(t_C <= pmax(t1 - rec_C, 0), t_E <= pmax(t1 - rec_E, 0)),
    group = c(rep(0, n_C), rep(1, n_E))
  )
  interim_df <- interim_df[interim_df$time > 0, ]

  # Final analysis data
  final_df <- data.frame(
    time = c(pmin(t_C, pmax(t2 - rec_C, 0)), pmin(t_E, pmax(t2 - rec_E, 0))),
    status = c(t_C <= pmax(t2 - rec_C, 0), t_E <= pmax(t2 - rec_E, 0)),
    group = c(rep(0, n_C), rep(1, n_E))
  )

  # Compute Kaplan-Meier estimates via survfit
  surv_int_E <- survfit(
    formula = Surv(time, status) ~ 1,
    data = subset(interim_df, group == 1)
  )
  surv_int_C <- survfit(
    formula = Surv(time, status) ~ 1,
    data = subset(interim_df, group == 0)
  )
  surv_fin_E <- survfit(
    formula = Surv(time, status) ~ 1,
    data = subset(final_df, group == 1)
  )
  surv_fin_C <- survfit(
    formula = Surv(time, status) ~ 1,
    data = subset(final_df, group == 0)
  )

  # Compute Kaplan-Meier integrals from estimates
  km_integral_values_EC_interim[run] <- -km_integral(
    surv_int_E,
    surv_int_C,
    lower = 0,
    upper = eps
  )
  km_integral_values_EC_final[run] <- -km_integral(
    surv_fin_E,
    surv_fin_C,
    lower = 0,
    upper = eps
  )

  km_integral_values_CE_interim[run] <- -km_integral(
    surv_int_C,
    surv_int_E,
    lower = 0,
    upper = eps
  )
  km_integral_values_CE_final[run] <- -km_integral(
    surv_fin_C,
    surv_fin_E,
    lower = 0,
    upper = eps
  )

  # Extract estimates at end of time window
  survsum_int_C <- summary(surv_int_C, times = eps)
  km_int_C_eps <- survsum_int_C$surv
  nasd_int_C_eps <- survsum_int_C$std.chaz
  survsum_int_E <- summary(surv_int_E, times = eps)
  km_int_E_eps <- survsum_int_E$surv
  nasd_int_E_eps <- survsum_int_E$std.chaz
  survsum_fin_C <- summary(surv_fin_C, times = eps)
  km_fin_C_eps <- survsum_fin_C$surv
  nasd_fin_C_eps <- survsum_fin_C$std.chaz
  survsum_fin_E <- summary(surv_fin_E, times = eps)
  km_fin_E_eps <- survsum_fin_E$surv
  nasd_fin_E_eps <- survsum_fin_E$std.chaz

  # Compute variance estimates for Kaplan-Meier integrals
  km_integral_var_CE_interim[run] <-
    km_int_C_eps^2 *
    (km_int_E_eps * nasd_int_E_eps)^2 -
    2 *
      km_int_C_eps *
      km_int_E_eps *
      var_single_integral(
        surv_int_E,
        surv_int_C,
        lower = 0,
        upper = eps,
        vartype = "Greenwood"
      ) +
    var_double_integral(
      surv_int_E,
      surv_int_C,
      lower = 0,
      upper = eps,
      vartype = "Greenwood"
    ) +
    var_double_integral(
      surv_int_C,
      surv_int_E,
      lower = 0,
      upper = eps,
      vartype = "Greenwood"
    )

  km_integral_var_EC_interim[run] <-
    km_int_E_eps^2 *
    (km_int_C_eps * nasd_int_C_eps)^2 -
    2 *
      km_int_E_eps *
      km_int_C_eps *
      var_single_integral(
        surv_int_C,
        surv_int_E,
        lower = 0,
        upper = eps,
        vartype = "Greenwood"
      ) +
    var_double_integral(
      surv_int_C,
      surv_int_E,
      lower = 0,
      upper = eps,
      vartype = "Greenwood"
    ) +
    var_double_integral(
      surv_int_E,
      surv_int_C,
      lower = 0,
      upper = eps,
      vartype = "Greenwood"
    )

  km_integral_var_CE_final[run] <-
    km_fin_C_eps^2 *
    (km_fin_E_eps * nasd_fin_E_eps)^2 -
    2 *
      km_fin_C_eps *
      km_fin_E_eps *
      var_single_integral(
        surv_fin_E,
        surv_fin_C,
        lower = 0,
        upper = eps,
        vartype = "Greenwood"
      ) +
    var_double_integral(
      surv_fin_E,
      surv_fin_C,
      lower = 0,
      upper = eps,
      vartype = "Greenwood"
    ) +
    var_double_integral(
      surv_fin_C,
      surv_fin_E,
      lower = 0,
      upper = eps,
      vartype = "Greenwood"
    )

  km_integral_var_EC_final[run] <-
    km_fin_E_eps^2 *
    (km_fin_C_eps * nasd_fin_C_eps)^2 -
    2 *
      km_fin_E_eps *
      km_fin_C_eps *
      var_single_integral(
        surv_fin_C,
        surv_fin_E,
        lower = 0,
        upper = eps,
        vartype = "Greenwood"
      ) +
    var_double_integral(
      surv_fin_C,
      surv_fin_E,
      lower = 0,
      upper = eps,
      vartype = "Greenwood"
    ) +
    var_double_integral(
      surv_fin_E,
      surv_fin_C,
      lower = 0,
      upper = eps,
      vartype = "Greenwood"
    )

  # Compute variance estimates rot rNTB
  rntb_var_interim[run] <-
    km_int_C_eps^2 *
    (km_int_E_eps * nasd_int_E_eps)^2 -
    4 *
      km_int_C_eps *
      km_int_E_eps *
      var_single_integral(
        surv_int_E,
        surv_int_C,
        lower = 0,
        upper = eps,
        vartype = "Greenwood"
      ) +
    4 *
      var_double_integral(
        surv_int_C,
        surv_int_E,
        lower = 0,
        upper = eps,
        vartype = "Greenwood"
      ) +
    km_int_E_eps^2 * (km_int_C_eps * nasd_int_C_eps)^2 -
    4 *
      km_int_E_eps *
      km_int_C_eps *
      var_single_integral(
        surv_int_C,
        surv_int_E,
        lower = 0,
        upper = eps,
        vartype = "Greenwood"
      ) +
    4 *
      var_double_integral(
        surv_int_E,
        surv_int_C,
        lower = 0,
        upper = eps,
        vartype = "Greenwood"
      )

  rntb_cov[run] <- rntb_var_final[run] <-
    km_fin_C_eps^2 *
    (km_fin_E_eps * nasd_fin_E_eps)^2 -
    4 *
      km_fin_C_eps *
      km_fin_E_eps *
      var_single_integral(
        surv_fin_E,
        surv_fin_C,
        lower = 0,
        upper = eps,
        vartype = "Greenwood"
      ) +
    4 *
      var_double_integral(
        surv_fin_C,
        surv_fin_E,
        lower = 0,
        upper = eps,
        vartype = "Greenwood"
      ) +
    km_fin_E_eps^2 * (km_fin_C_eps * nasd_fin_C_eps)^2 -
    4 *
      km_fin_E_eps *
      km_fin_C_eps *
      var_single_integral(
        surv_fin_C,
        surv_fin_E,
        lower = 0,
        upper = eps,
        vartype = "Greenwood"
      ) +
    4 *
      var_double_integral(
        surv_fin_E,
        surv_fin_C,
        lower = 0,
        upper = eps,
        vartype = "Greenwood"
      )
}

# Compare empirical variances of Kaplan-Meier integrals with estimated variances
var(sqrt(n_E) * km_integral_values_CE_interim)
mean(n_E * km_integral_var_CE_interim)
var(sqrt(n_E) * km_integral_values_EC_interim)
mean(n_E * km_integral_var_EC_interim)
var(sqrt(n_E) * km_integral_values_CE_final)
mean(n_E * km_integral_var_CE_final)
var(sqrt(n_E) * km_integral_values_EC_final)
mean(n_E * km_integral_var_EC_final)

# Compute rNTB estimates from probabilistic indices
rntb_interim <- km_integral_values_EC_interim - km_integral_values_CE_interim
rntb_final <- km_integral_values_EC_final - km_integral_values_CE_final

# Compare empirical variances of rNTB estimates with estimated variances
var(sqrt(n_E) * rntb_interim)
mean(n_E * rntb_var_interim)
var(sqrt(n_E) * rntb_final)
mean(n_E * rntb_var_final)

# Compute empirical levels of corresponding tests
mean(rntb_interim / sqrt(rntb_var_interim) < qnorm(0.025))
mean(rntb_final / sqrt(rntb_var_final) < qnorm(0.025))
