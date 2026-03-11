# Computes true probabilistic index for two exponential distributions
true_km_integral_exp <- function(
  rate_integrand,
  rate_integrator,
  lower,
  upper
) {
  helper_function <- function(u) {
    pexp(u, rate = rate_integrand, lower.tail = FALSE) *
      (-dexp(u, rate = rate_integrator))
  }
  return(integrate(helper_function, lower = lower, upper = upper)$value)
}

# Computes Kaplan-Meier integrals
km_integral <- function(integrand, integrator, lower, upper) {
  # integrand and integrator are both survfit objects
  # lower and upper are reals
  # We compute the integral S_integrand d(S_integrator)
  #   from lower to upper where S are KM-estimates

  if (typeof(lower) == "list") {
    lower <- unlist(lower)
  }
  if (typeof(upper) == "list") {
    upper <- unlist(upper)
  }

  # Collect jump times of integrator in the interval
  all_jump_times <- c(0, integrator$time)
  jump_indices <- which(all_jump_times > lower & all_jump_times <= upper)

  # If no jumps occur in this time window, integral is 0
  if (length(jump_indices) == 0) {
    return(0)
  }

  # Find jump sizes at these times
  all_values <- c(1, integrator$surv)
  jump_sizes <- all_values[jump_indices] - all_values[jump_indices - 1]

  # Find values of integrand at jump times of integrator
  integrand_val <- summary(integrand, times = all_jump_times[jump_indices])$surv

  return(jump_sizes %*% integrand_val)
}

# NOTE:
# std.err of survfit object is the Greenwood estimator of s.d. of the cum. haz.
# std.chaz of survfit object is the standard estimator of s.d. of the cum. haz.
# std.err of summary.survfit object is the estimator of s.d. of the survival function
#  = product of std.err of survfit object and surv of survfit object
# std.chaz of summary.survfit object is the standard estimator of s.d. of the cum. haz.

# TODO:
# So far, we always use the standard estimator. Variance is often underestimated.
# Should also try Greenwood instead

# Evaluates simple integrals for variance estimation
var_single_integral <- function(
  integrand,
  integrator,
  lower,
  upper,
  vartype = "std"
) {
  if (typeof(lower) == "list") {
    lower <- unlist(lower)
  }
  if (typeof(upper) == "list") {
    upper <- unlist(upper)
  }

  # Collect jump times of integrator in the interval
  all_jump_times <- c(0, integrator$time)
  jump_indices <- which(all_jump_times > lower & all_jump_times <= upper)
  jump_times <- all_jump_times[jump_indices]

  # If no jumps occur in this time window, integral is 0
  if (length(jump_indices) == 0) {
    return(0)
  }

  # Find jump sizes at these times
  all_values <- c(1, integrator$surv)
  jump_sizes <- all_values[jump_indices] - all_values[jump_indices - 1]

  survsum <- summary(integrand, times = jump_times)
  if (vartype == "std") {
    # Find values of integrand at jump times of integrator
    integrand_surv_val <- survsum$surv
    integrand_haz_sd_val <- survsum$std.chaz^2
  } else if (vartype == "Greenwood") {
    # Find values of integrand at jump times of integrator
    integrand_surv_val <- survsum$surv
    # Recover Greenwood estimator
    integrand_haz_sd_val <- (survsum$std.err / survsum$surv)^2
  }
  integrand_val <- integrand_haz_sd_val * integrand_surv_val

  return(jump_sizes %*% integrand_val)
}

# Evaluates double integrals for variance estimation
var_double_integral <- function(
  integrand,
  integrator,
  lower,
  upper,
  vartype = "std"
) {
  if (typeof(lower) == "list") {
    lower <- unlist(lower)
  }
  if (typeof(upper) == "list") {
    upper <- unlist(upper)
  }

  # Collect jump times of integrator in the interval
  all_jump_times <- c(0, integrator$time)
  jump_indices <- which(all_jump_times > lower & all_jump_times <= upper)
  jump_times <- all_jump_times[jump_indices]
  num_jump_times <- length(jump_times)
  jump_time_matrix <- outer(jump_times, jump_times, FUN = pmin)

  # If no jumps occur in this time window, integral is 0
  if (length(jump_indices) == 0) {
    return(0)
  }

  # Find jump sizes at these times
  all_values <- c(1, integrator$surv)
  jump_sizes <- all_values[jump_indices] - all_values[jump_indices - 1]

  if (vartype == "std") {
    # Find values of integrand at jump times of integrator
    integrand_surv_val <- summary(integrand, times = jump_times)$surv
    integrand_haz_sd_val <- matrix(
      summary(integrand, times = jump_time_matrix)$std.chaz^2,
      nrow = num_jump_times,
      ncol = num_jump_times
    )
  } else if (vartype == "Greenwood") {
    # Find values of integrand at jump times of integrator
    survsum <- summary(integrand, times = jump_times)
    integrand_surv_val <- survsum$surv

    # Recover Greenwood estimator
    greenwood_var_vec <- (survsum$std.err / survsum$surv)^2
    integrand_haz_sd_val <- outer(
      1:num_jump_times,
      1:num_jump_times,
      FUN = function(i, j) greenwood_var_vec[pmin(i, j)]
    )
  }

  integrand_surv_sd_val <- integrand_haz_sd_val *
    outer(integrand_surv_val, integrand_surv_val, "*")
  return(jump_sizes %*% integrand_surv_sd_val %*% jump_sizes)
}
