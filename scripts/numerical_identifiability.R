betas01_t = c(trt = -0.2)
betas02_t = c(trt =-0.4)
betas12_t = c(trt =-0.5)
lambdas01_t = 0.1
lambdas02_t = 0.2
lambdas12_t = 0.3
gammas01_t = 1.5
gammas02_t = 1
gammas12_t = 2
cens = c(4.5, 5.5)

set.seed(9911)
covs <- data.frame(id = 1:2000, trt = stats::rbinom(2000, 1L, 0.5))

sim_wei <- rsimid(
  dist01 = "weibull",
  dist02 = "weibull",
  dist12 = "weibull",
  betas01 = betas01_t,
  betas02 = betas02_t,
  betas12 = betas12_t,
  lambdas01 = lambdas01_t,
  lambdas02 = lambdas02_t,
  lambdas12 = lambdas12_t,
  gammas01 = gammas01_t,
  gammas02 = gammas02_t,
  gammas12 = gammas12_t,
  x = covs,
  cens = cens
)

library(survival)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = 2)
stanfit <- idm_stan(formula01 = Surv(time=df_time,event=df_event)~trt,
                     formula02 = Surv(time=os_time,event=os_event)~trt,
                     formula12 = Surv(time=time_diff,event=os_event)~trt,
                     data = sim_wei,
                     basehaz01 = "weibull",
                     basehaz02 = "weibull",
                     basehaz12 = "weibull",
                     prior01           = rstanarm::normal(),
                     prior_intercept01 = rstanarm::normal(),
                     prior_aux01       = rstanarm::normal(),
                     prior02           = rstanarm::normal(),
                     prior_intercept02 = rstanarm::normal(),
                     prior_aux02       = rstanarm::normal(),
                     prior12           = rstanarm::normal(),
                     prior_intercept12 = rstanarm::normal(),
                     prior_aux12       = rstanarm::normal(),
                    iter = 2000,
                    chains = 4,
                    control = list(adapt_delta = 0.99)
)




stanpars <- c(if (standata$has_intercept01) "alpha01",
              if (standata$K01)             "beta01",
              if (standata$nvars01)         "aux01",
              if (standata$has_intercept02) "alpha02",
              if (standata$K02)             "beta02",
              if (standata$nvars02)         "aux02",
              if (standata$has_intercept12) "alpha12",
              if (standata$K12)             "beta12",
              if (standata$nvars12)         "aux12")

stanfile <-  "src/stan_files/MS.stan"
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = 2)
fit <- stan(file = stanfile,
            data = standata,
            pars = stanpars,
            iter = 2000,
            chains = 4,
            control = list(adapt_delta = 0.99))
s_elapsed <- sum(get_elapsed_time(fit))


print(fit)
rstan::traceplot(fit, 'lp__')

rstan::traceplot(fit, c('beta01','beta02', 'beta12'), ncol = 1)

library(cowplot)

## ----sim-extract-alpha---------------------------------------------------
pp_lambda01 <- exp( rstan::extract(fit,'alpha01')$alpha01 )
pp_gamma01 <- rstan::extract(fit,'aux01')$aux01
pp_beta01 <- rstan::extract(fit,'beta01')$beta01

## ----plot-alpha-vs-test--------------------------------------------------
ggplot(data.frame(lambda = pp_lambda01, gamma = pp_gamma01)) +
  geom_density(aes(x = lambda)) +
  geom_vline(aes(xintercept = lambdas01_t), colour = 'red') +
  ggtitle('Posterior distribution of alpha\nshowing true value in red')

## ----plot-mu-vs-test-----------------------------------------------------
ggplot(data.frame(lambda = pp_lambda01, gamma = pp_gamma01)) +
  geom_density(aes(x = gamma)) +
  geom_vline(aes(xintercept = gammas01_t), colour = 'red') +
  ggtitle('Posterior distribution of mu\nshowing true value in red')

## ----plot-mu-vs-alpha----------------------------------------------------
p1 <- ggplot(data.frame(lambda = pp_lambda01, gamma = pp_gamma01)) +
  geom_density2d(aes(x = lambda, y = gamma), linetype = "dashed") +
  geom_point(aes(x = lambdas01_t, y = gammas01_t), colour = 'red', size = 4) +
  ggtitle('Posterior distributions of lambda and gamma\nshowing true parameter values')

## ----plot-beta-vs-test-----------------------------------------------------
p2 <- ggplot(data.frame(beta = pp_beta01)) +
  geom_density(aes(x = beta)) +
  geom_vline(aes(xintercept = betas01_t[1]), colour = 'red') +
  ggtitle('Posterior distribution of beta\nshowing true value')

cowplot::plot_grid(p1, p2, labels = c("Theta", "beta"))

## ------------------------------------------------------------------------
mean(pp_lambda01 >= lambdas01_t)

## ------------------------------------------------------------------------
mean(pp_gamma01 >= gammas01_t)

## ------------------------------------------------------------------------
mean(pp_lambda01 >= lambdas01_t & pp_gamma01 >= gammas01_t)



