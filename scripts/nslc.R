roxygen2::roxygenise(clean = TRUE)

scl <- sas7bdat::read.sas7bdat("/media/mtr/A5C2-009E/SCL/c9732_demographic.sas7bdat")

scl <- within(scl, {
  PD_STATUS <- ifelse(is.nan(PD_TIME), 0, 1)
  OS_STATUS <- ifelse(STATUS == 2, 1, 0)
  PD_TIME = PFS_TIME
  TIMEDIFF = OS_TIME - PD_TIME
  TRTARM = as.factor(TRT_ARM)
})

scl$TIMEDIFF[scl$PD_STATUS == 1]

library(caret)
train.index <- createDataPartition(scl$PD_STATUS & scl$OS_STATUS, p = .7, list = FALSE)
train <- scl[ train.index,]
test  <- scl[-train.index,]

nrow(test)
nrow(train)


library(rstan)
options(mc.cores=2)


stanfit_2 <- idm_stan(
  formula01 = Surv(time=PFS_TIME,event=PD_STATUS)~TRTARM,
  formula02 = Surv(time=OS_TIME,event=OS_STATUS)~TRTARM,
  formula12 = Surv(time=TIMEDIFF,event=OS_STATUS)~TRTARM,
  data = train,
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

newdata

newdata <- handle_newdata(formula01 = Surv(time=PFS_TIME,event=PFS_STATUS)~TRTARM,
                     formula02 = Surv(time=OS_TIME,event=OS_STATUS)~1,
                    formula12 = Surv(time=TIMEDIFF,event=OS_STATUS)~TRTARM,
                            data = test,
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

stanfit_2


plotfit <- posterior_fit(stanfit_2,
                       times = 0,
                       extrapolate = TRUE,
                       standardise = TRUE,
                       type = "cumhaz")

plot(plotfit)

posterior_fit(stanfit_2,
              newdata = newdata,
              times = 0,
              extrapolate = TRUE, condition = FALSE, control = list(edist = 5), type = "cumhaz")

#### To conduct the analysis of survival brier score

msfit <- posterior_fit(stanfit_2,
                       newdata = newdata,
                       times = 0,
                       extrapolate = TRUE,
                       condition = FALSE,
                       control = list(edist = 65),
                       type = "cumhaz")

msevent <- function(object){

  tmat <- mstate::trans.illdeath()

  out <- list()

  for(i in unique(object$id)){
    haz <-  as.data.frame(msfit) %>%
      dplyr::filter(id == i) %>%
      dplyr::mutate(Haz = median,
                    trans = transition) %>%
      dplyr::select(time, Haz, trans)

    tv <- sort(unique(haz$time))
    out[[i]] <- mssample(Haz=haz,
                    trans=tmat,
                    clock = "reset",
                    M=1000,
                    output = "state",
                    tvec = tv,
                    do.trace=500)
  }

  out
}

fitevent <- msevent(msfit)



print(fit)
rstan::traceplot(fit, 'lp__')

rstan::traceplot(fit, c('beta01','beta02', 'beta12'), ncol = 1)
