roxygen2::roxygenise(clean = TRUE)
data("sclc_demo")
sclc_demo$TIMEDIFF <- sclc_demo$OS_TIME - sclc_demo$PFS_TIME
sclc_demo$STATUS[sclc_demo$STATUS == 1] <- 0
sclc_demo$STATUS[sclc_demo$STATUS == 3] <- 1
sclc_demo$STATUS[sclc_demo$STATUS == 2] <- 1
# sclc_demo <- sclc_demo[sclc_demo$OS_TIME > 0, ]
# sclc_demo <- sclc_demo[sclc_demo$PFS_TIME > 0, ]

sclc_demo$TIMEDIFF[sclc_demo$PFS_STATUS == 1]
library(rstan)
options(mc.cores=2)

sclc_demo$TRTARM <- "G"
sclc_demo$TRTARM[sclc_demo$TRT_ARM == 1] <- "CE"

stanfit_2 <- idm_stan(
  formula01 = Surv(time=PFS_TIME,event=PFS_STATUS)~TRTARM,
  formula02 = Surv(time=OS_TIME,event=STATUS)~1,
  formula12 = Surv(time=TIMEDIFF,event=STATUS)~TRTARM,
  data = sclc_demo,
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
  prior_aux12       = rstanarm::normal()
)

plot(stanfit)

loo1 <- loo(stanfit, cores = 2, k_threshold = 0.7)
loo2 <- loo(stanfit_1, cores = 2, k_threshold = 0.7)
loo3 <- loo(stanfit_2, cores = 2, k_threshold = 0.7)
as.data.frame( compare_models(loo1, loo2, loo3) )

stanfit_2
msfit <- posterior_fit(stanfit_2,
                       times = 0,
                       extrapolate = TRUE,
                       condition = FALSE,
                       control = list(edist = 65),
                       type = "cumhaz")

plotfit <- posterior_fit(stanfit_2,
                       times = 0,
                       extrapolate = TRUE,
                       standardise = TRUE,
                       type = "cumhaz")

plot(plotfit)



posterior_fit(stanfit2,
              newdata = newdata,
              times = 0,
              extrapolate = TRUE, condition = FALSE, control = list(edist = 5), type = "cumhaz")

#### To conduct the analysis of survival brier score

library(dplyr)
haz <-  as.data.frame(msfit) %>%
  filter(id == 1) %>%
  dplyr::mutate(Haz = median,
                trans = transition) %>%
  select(time, Haz, trans)

library(mstate)
# transition matrix for illness-death model
tmat <- trans.illdeath()

tv <- sort(unique(haz$time))
out <- mssample(Haz=haz,
                trans=tmat,
                clock = "reset",
                M=1000,
                output = "state",
                tvec = tv,
                do.trace=500)




print(fit)
rstan::traceplot(fit, 'lp__')

rstan::traceplot(fit, c('beta01','beta02', 'beta12'), ncol = 1)
