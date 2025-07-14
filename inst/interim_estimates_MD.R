library(survival)
library(BuyseTest)
library(DelayedGSD)

source("R/plot.R")
source("R/runTrial.R")
source("R/simTrial.R")
source("R/subset.R")
source("R/summary.R")

n.E <- 500
n.C <- 500

# Patients followed-up until 0.5 years
h <- 0.5

# Recruitment over period of 1 year
a <- 1

# Interim anlysis after 0.75 years:
#  - 1/4 of patients have complete follow-up (recruited before 0.25)
#  - 1/2 of patients already recruited but incomplete follow-up (recruited before 0.75)
t_int <- 1
a_int.comp <- t_int - h

runs <- 1000
my_seeds <- runif(runs)

estimates <- matrix(NA, ncol = 4, nrow = runs)
colnames(estimates) <- c("all", "complete_interim", "all_interim_raw", "all_interim_change")

# Raw pairwise comparisons
# Can be used for:
#  - Comparison of complete observations
#  - 'First summand' in new procedure (including pipeline patients)
raw_comparison_surv <- function(df){

  n_C <- sum(df$group == "control")
  n_E <- sum(df$group == "experimental")

  df <- df[order(-df$statusSurv, df$timeSurv),]
  df$c.countSurv <- cumsum((df$group == "control") * (df$statusSurv))
  e.sup_count <- sum(df$c.countSurv[df$group == "experimental"])
  return(e.sup_count/(n_C * n_E))

}

# Compute empirical probability of possible order changes for pipeline patients
order_changes <- function(df){

  n_C <- sum(df$group == "control")
  n_E <- sum(df$group == "experimental")

  recdates.C <- df$timeInclusion[df$group == "control"]
  # Extract recruitment dates of pipeline patients
  pipe_recdates.C <- recdates.C[recdates.C > a_int.comp & recdates.C <= t_int]

  surv.E <- survfit(formula = Surv(timeSurv, statusSurv) ~ 1,
                    data = df[df$group == "experimental",])
  surv.C <- survfit(formula = Surv(timeSurv, statusSurv) ~ 1,
                    data = df[df$group == "control",])

  km_final.E <- summary(surv.E, times = h)$surv
  km_final.C <- summary(surv.C, times = h)$surv

  km_recdates.C <- summary(surv.C, times = t_int - pipe_recdates.C)$surv
  return(sum((1 - (km_final.C/km_recdates.C)) * km_final.E)/n_C)

}

for(run in 1:runs){

  # Simulate trial
  dfH0.delay <- simTrial(n.E = n.E, n.C = n.C, scale.censoring.C = Inf, admin.censoring = h)

  # complete observations at interim analysis (no additional administrative censoring necessary)
  dfH0.delay_int.comp <- dfH0.delay[dfH0.delay$timeInclusion <= a_int.comp,]
  # recruited before interim analysis (apply additional administrative censoring)
  dfH0.delay_int.rec <- dfH0.delay[dfH0.delay$timeInclusion <= t_int,]
  dfH0.delay_int.rec$statusSurv <- dfH0.delay_int.rec$statusSurv *
    ifelse(dfH0.delay_int.rec$timeInclusion + dfH0.delay_int.rec$timeSurv <= t_int, 1, 0)
  dfH0.delay_int.rec$timeSurv <- pmin(dfH0.delay_int.rec$timeSurv, pmax(0, t_int - dfH0.delay_int.rec$timeInclusion))
  dfH0.delay_int.rec$statusTox <- ifelse(dfH0.delay_int.rec$statusTox == 0, 0,
                                         ifelse(dfH0.delay_int.rec$statusTox == 1 & dfH0.delay_int.rec$timeInclusion + dfH0.delay_int.rec$timeTox <= t_int, 1,
                                                ifelse(dfH0.delay_int.rec$statusSurv == 1, 2, 0)))


  estimates[run, "all"] <- raw_comparison_surv(dfH0.delay)
  estimates[run, "complete_interim"] <- raw_comparison_surv(dfH0.delay_int.comp)
  estimates[run, "all_interim_raw"] <- raw_comparison_surv(dfH0.delay_int.rec)
  estimates[run, "all_interim_change"] <- order_changes(dfH0.delay_int.rec)

}

1/var(estimates[, "all"])
1/var(estimates[, "complete_interim"])
1/var(estimates[, "all_interim_raw"] + estimates[, "all_interim_change"])
