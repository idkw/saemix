#####################################################################
# Bioequivalence 2x2 crossover trial — Simulation + fit with IOV
#####################################################################

# Load the installed package
library(saemix)

# =================================================================
# 1. Simulate data
# =================================================================
set.seed(12345)

# Design
n_per_seq <- 24
n_subjects <- n_per_seq * 2  # 48 subjects
times <- c(0.5, 1, 1.5, 2, 3, 4, 6, 8, 12, 24)
dose <- 100  # mg

# True population parameters
ka_pop  <- 1.5    # absorption rate (h^-1)
V_pop   <- 30     # volume (L)
CL_pop  <- 5      # clearance (L/h)
Frel_pop <- 1.10  # relative bioavailability test vs reference (+10%)

# Variability
omega_ka  <- 0.30  # IIV on log(ka)
omega_V   <- 0.20  # IIV on log(V)
omega_CL  <- 0.15  # IIV on log(CL)
psi_ka    <- 0.20  # IOV on log(ka)
sigma_prop <- 0.10 # proportional residual error (10%)

dat_list <- vector("list", n_subjects * 2 * length(times))
idx <- 0L
for (i in 1:n_subjects) {
  seq <- ifelse(i <= n_per_seq, "RT", "TR")

  # IIV
  ka_i  <- ka_pop * exp(rnorm(1, 0, omega_ka))
  V_i   <- V_pop  * exp(rnorm(1, 0, omega_V))
  CL_i  <- CL_pop * exp(rnorm(1, 0, omega_CL))
  k_i   <- CL_i / V_i

  for (period in 1:2) {
    trt <- if (seq == "RT") as.integer(period == 2) else as.integer(period == 1)

    # IOV on ka
    ka_ik <- ka_i * exp(rnorm(1, 0, psi_ka))
    F_ik  <- ifelse(trt == 1, Frel_pop, 1.0)

    for (t in times) {
      idx <- idx + 1L
      cpred <- F_ik * dose * ka_ik / (V_i * (ka_ik - k_i)) *
               (exp(-k_i * t) - exp(-ka_ik * t))
      cobs <- max(cpred * (1 + rnorm(1, 0, sigma_prop)), 0.01)
      dat_list[[idx]] <- data.frame(Id = i, Time = t, Dose = dose,
                                    Concentration = round(cobs, 4),
                                    Treatment = trt, Period = period)
    }
  }
}
dat <- do.call(rbind, dat_list)

cat("Dataset: ", nrow(dat), "observations,", n_subjects, "subjects,",
    "2 periods, 10 times\n")
cat("  Subjects 1-24: sequence RT,  subjects 25-48: sequence TR\n\n")

# =================================================================
# 2. Create saemix objects
# =================================================================

saemix.data <- saemixData(
  name.data      = dat,
  header         = TRUE,
  name.group     = "Id",
  name.predictors = c("Dose", "Time", "Treatment"),
  name.response  = "Concentration",
  name.occ       = "Period",
  name.X         = "Time",
  units          = list(x = "hr", y = "mg/L")
)

# Structural model: 1-compartment with relative bioavailability
model_1cpt_be <- function(psi, id, xidep) {
  dose <- xidep[, 1]
  time <- xidep[, 2]
  trt  <- xidep[, 3]

  ka   <- psi[id, 1]
  V    <- psi[id, 2]
  CL   <- psi[id, 3]
  Frel <- psi[id, 4]

  k <- CL / V
  Fbio <- 1 + trt * (Frel - 1)
  ypred <- Fbio * dose * ka / (V * (ka - k)) * (exp(-k * time) - exp(-ka * time))
  return(ypred)
}

saemix.model <- saemixModel(
  model      = model_1cpt_be,
  psi0       = matrix(c(1.5, 30, 5, 1.1), nrow = 1,
                      dimnames = list(NULL, c("ka", "V", "CL", "Frel"))),
  transform.par    = c(1, 1, 1, 1),           # all log-normal
  covariance.model = diag(c(1, 1, 1, 0)),     # IIV on ka, V, CL; not on Frel
  covariance.model.iov = diag(c(1, 0, 0, 0)), # IOV on ka only
  error.model      = "proportional"
)

# =================================================================
# 3. Run SAEM with IOV
# =================================================================
cat("Running SAEM algorithm with IOV...\n")
saemix.fit <- saemix(saemix.model, saemix.data,
                     list(seed = 12345, nbiter.saemix = c(300, 100),
                          save = FALSE, save.graphs = FALSE, print = FALSE,
                          displayProgress = FALSE))

# =================================================================
# 4. Results
# =================================================================
cat("\n========================================\n")
cat("  RESULTS\n")
cat("========================================\n\n")
print(saemix.fit@results)

cat("\n--- Comparison with true values ---\n")
cat(sprintf("  ka   : true=%.2f  est=%.2f\n", ka_pop, saemix.fit@results@fixed.psi[1]))
cat(sprintf("  V    : true=%.1f  est=%.1f\n", V_pop, saemix.fit@results@fixed.psi[2]))
cat(sprintf("  CL   : true=%.1f  est=%.1f\n", CL_pop, saemix.fit@results@fixed.psi[3]))
cat(sprintf("  Frel : true=%.2f  est=%.2f\n", Frel_pop, saemix.fit@results@fixed.psi[4]))
cat(sprintf("  omega2.ka : true=%.3f  est=%.3f\n", omega_ka^2, diag(saemix.fit@results@omega)[1]))
cat(sprintf("  omega2.V  : true=%.3f  est=%.3f\n", omega_V^2, diag(saemix.fit@results@omega)[2]))
cat(sprintf("  omega2.CL : true=%.3f  est=%.3f\n", omega_CL^2, diag(saemix.fit@results@omega)[3]))
cat(sprintf("  psi2.ka (IOV) : true=%.3f  est=%.3f\n", psi_ka^2, diag(saemix.fit@results@psi.iov)[1]))
cat(sprintf("  sigma_prop : true=%.3f  est=%.3f\n", sigma_prop, saemix.fit@results@respar[2]))

cat("\nOccasion fixed effects (beta.occ):\n")
print(saemix.fit@results@beta.occ)