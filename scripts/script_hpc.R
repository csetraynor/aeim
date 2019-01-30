roxygen2::roxygenise(clean = TRUE)

Sys.time()
betas01_t = c(trt = -0.22)
betas02_t = c(trt = -0.4)
betas12_t = c(trt = -0.33)

fn <- function(t, x, betas, ...){
  (-1 + 0.02 * t - 0.03 * t ^ 2 + 0.005 * t ^ 3) + (x * betas)
}

cens = c(5.25, 7.75)

set.seed(9911)
covs <- data.frame(id = 1:2000, trt = stats::rbinom(2000, 1L, 0.5))

sim_rp <- rsimid(
  betas01 = betas01_t,
  betas02 = betas02_t,
  betas12 = betas12_t,
  loghazard01 = fn,
  loghazard02 = fn,
  loghazard12 = fn,
  x = covs,
  cens = cens
)
Sys.time()

sim_rp$time_diff = sim_rp$os_time - sim_rp$df_time

saveRDS(sim_rp, "rp1var.RDS")

##################################################################

Sys.time()
betas01_t = c(trt = -0.22, age = -0.11)
betas02_t = c(trt = -0.4, age = 0.22)
betas12_t = c(trt = -0.33, age = 0.18)

fn <- function(t, x, betas, ...){
  (-1 + 0.02 * t - 0.03 * t ^ 2 + 0.005 * t ^ 3) + (x * betas)
}

cens = c(5.25, 7.75)

set.seed(9911)
covs <- data.frame(id = 1:2000, trt = stats::rbinom(2000, 1L, 0.5),
                   age = rnorm(2000, mean = 50, sd = 10))

sim_rp <- rsimid(
  betas01 = betas01_t,
  betas02 = betas02_t,
  betas12 = betas12_t,
  loghazard01 = fn,
  loghazard02 = fn,
  loghazard12 = fn,
  x = covs,
  cens = cens
)
Sys.time()

library(survival)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = 2)

sim_rp$time_diff = sim_rp$os_time - sim_rp$df_time

saveRDS(sim_rp, "rp2var.RDS")

roxygen2::roxygenise(clean = TRUE)

Sys.time()
betas01_t = c(trt = -0.22)
betas02_t = c(trt = -0.4)
betas12_t = c(trt = -0.33)

fn <- function(t, x, betas, ...){
  (-1 + 0.02 * t - 0.03 * t ^ 2 + 0.005 * t ^ 3) + (x * betas)
}

cens = c(5.25, 7.75)

set.seed(9911)
covs <- data.frame(id = 1:20000, trt = stats::rbinom(20000, 1L, 0.5))

sim_rp <- rsimid(
  betas01 = betas01_t,
  betas02 = betas02_t,
  betas12 = betas12_t,
  loghazard01 = fn,
  loghazard02 = fn,
  loghazard12 = fn,
  x = covs,
  cens = cens
)
Sys.time()

sim_rp$time_diff = sim_rp$os_time - sim_rp$df_time

saveRDS(sim_rp, "rp10var.RDS")

##################################################################

Sys.time()
betas01_t = c(trt = -0.22, age = -0.11)
betas02_t = c(trt = -0.4, age = 0.22)
betas12_t = c(trt = -0.33, age = 0.18)

fn <- function(t, x, betas, ...){
  (-1 + 0.02 * t - 0.03 * t ^ 2 + 0.005 * t ^ 3) + (x * betas)
}

cens = c(5.25, 7.75)

set.seed(9911)
covs <- data.frame(id = 1:20000, trt = stats::rbinom(20000, 1L, 0.5),
                   age = rnorm(20000, mean = 50, sd = 10))

sim_rp <- rsimid(
  betas01 = betas01_t,
  betas02 = betas02_t,
  betas12 = betas12_t,
  loghazard01 = fn,
  loghazard02 = fn,
  loghazard12 = fn,
  x = covs,
  cens = cens
)
Sys.time()

sim_rp$time_diff = sim_rp$os_time - sim_rp$df_time

saveRDS(sim_rp, "rp20var.RDS")
