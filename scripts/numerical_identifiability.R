roxygen2::roxygenise(clean = TRUE)

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

sim_wei$time_diff = sim_wei$os_time - sim_wei$df_time

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
                    iter = 1000,
                    chains = 2,
                    control = list(adapt_delta = 0.8)
)


stanfit2 <- idm_stan(formula01 = Surv(time=df_time,event=df_event)~trt,
                    formula02 = Surv(time=os_time,event=os_event)~trt,
                    formula12 = Surv(time=time_diff,event=os_event)~trt,
                    data = sim_wei,
                    basehaz01 = "ms",
                    basehaz02 = "ms",
                    basehaz12 = "ms",
                    prior01           = rstanarm::normal(),
                    prior_intercept01 = rstanarm::normal(),
                    prior_aux01       = rstanarm::normal(),
                    prior02           = rstanarm::normal(),
                    prior_intercept02 = rstanarm::normal(),
                    prior_aux02       = rstanarm::normal(),
                    prior12           = rstanarm::normal(),
                    prior_intercept12 = rstanarm::normal(),
                    prior_aux12       = rstanarm::normal(),
                    iter = 1000,
                    chains = 2,
                    control = list(adapt_delta = 0.8)
)

stanfit3 <- idm_stan(formula01 = Surv(time=df_time,event=df_event)~trt,
                     formula02 = Surv(time=os_time,event=os_event)~trt,
                     formula12 = Surv(time=time_diff,event=os_event)~trt,
                     data = sim_wei,
                     basehaz01 = "exp",
                     basehaz02 = "exp",
                     basehaz12 = "exp",
                     prior01           = rstanarm::normal(),
                     prior_intercept01 = rstanarm::normal(),
                     prior_aux01       = rstanarm::normal(),
                     prior02           = rstanarm::normal(),
                     prior_intercept02 = rstanarm::normal(),
                     prior_aux02       = rstanarm::normal(),
                     prior12           = rstanarm::normal(),
                     prior_intercept12 = rstanarm::normal(),
                     prior_aux12       = rstanarm::normal(),
                     iter = 1000,
                     chains = 2,
                     control = list(adapt_delta = 0.8)
)


loo1 <- loo(stanfit, cores = 2, k_threshold = 0.7)
loo2 <- loo(stanfit2, cores = 2, k_threshold = 0.7)
loo3 <- loo(stanfit3, cores = 2, k_threshold = 0.7)
compare_models(loo1, loo2, loo3)



ps <- posterior_fit(stanfit, standardise = TRUE,  times = list(0,0,0), extrapolate = TRUE, control = list(edist = 5),
                    type = "cumhaz")

ps2 <- posterior_fit(stanfit,  times = 0, ids = c(7,13,15),
    extrapolate = TRUE, condition = FALSE, control = list(edist = 5),
    type = "cumhaz")

ps
plot(ps,
     labels = c("0 -> 1", "0 -> 2", "1 -> 2"))


ps2 <- posterior_fit(stanfit2, type = "haz")
ps2
plot(ps2,
     ids = lapply(seq_along(ps2), function(x) 1:6),
     xlab = list("Time(years)", "Time(years)", "Time since event 1 (years)"),
     labels = c("0 -> 1", "0 -> 2", "1 -> 2"))
# ps2 <- posterior_fit(stanfit2,  type = "surv")
sps3 <- posterior_fit(stanfit3,
                      type = "cumhaz" )

plot(sps3,
     ids = lapply(seq_along(ps2), function(x) 1:6),
     xlab = list("Time(years)", "Time(years)", "Time since event 1 (years)"),
     labels = c("0 -> 1", "0 -> 2", "1 -> 2"))
sps3


library(mstate)

# transition matrix for illness-death model
tmat <- trans.illdeath()
# data in wide format, for transition 1 this is dataset E1 of
# Therneau & Grambsch (T&G)
tg <- data.frame(illt=c(1,1,6,6,8,9),ills=c(1,0,1,1,0,1),
                 dt=c(5,1,9,7,8,12),ds=c(1,1,1,1,1,1),
                 x1=c(1,1,1,0,0,0),x2=c(6:1))
# data in long format using msprep
tglong <- msprep(time=c(NA,"illt","dt"),status=c(NA,"ills","ds"),
                 data=tg,keep=c("x1","x2"),trans=tmat)
# expanded covariates
tglong <- expand.covs(tglong,c("x1","x2"))
# Cox model with different covariate
cx <- coxph(Surv(Tstart,Tstop,status)~x1.1+x2.2+strata(trans),
            data=tglong,method="breslow")
# new data, to check whether results are the same for transition 1 as T&G
newdata <- data.frame(trans=1:3,x1.1=c(0,0,0),x2.2=c(0,1,0),strata=1:3)
fit <- msfit(cx,newdata,trans=tmat)
tv <- unique(fit$Haz$time)
# mssample
set.seed(1234)

library(dplyr)
haz = as.data.frame(ps2) %>%
  filter(id == 1) %>%
  mutate(Haz = median,
         trans = transition) %>%
  select(time, Haz, trans)

library(mstate)
# transition matrix for illness-death model
tmat <- trans.illdeath()


mssample(Haz=haz,
         trans=tmat,
         clock = "reset",
         M=10,
         tvec = tvec,
         output="data",do.trace=25)


set.seed(1234)
mssample(Haz=fit$Haz,trans=tmat,tvec=tv,M=100)

mstate::transMat(
  list(c(2,3),c(),c(2)),
  names = c("Disease-free", "Death", "Relapse")
  )


tt <- ps[[1]]
plot(tt)


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



