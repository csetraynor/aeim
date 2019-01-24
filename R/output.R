plot.prevalence.msm <- function(x, mintime=NULL, maxtime=NULL, timezero=NULL, initstates=NULL,
                                interp=c("start","midpoint"), censtime=Inf, subset=NULL,
                                covariates="population", misccovariates="mean",
                                piecewise.times=NULL, piecewise.covariates=NULL, xlab="Times",ylab="Prevalence (%)",
                                lwd.obs=1, lwd.exp=1, lty.obs=1, lty.exp=2,
                                col.obs="blue", col.exp="red", legend.pos=NULL,...){
  if (!inherits(x, "msm")) stop("expected x to be a msm model")
  time <- model.extract(x$data$mf, "time")
  if (is.null(mintime)) mintime <- min(time)
  if (is.null(maxtime)) maxtime <- max(time)
  t <- seq(mintime, maxtime, length=100)
  obs <- observed.msm(x, t, interp, censtime, subset)
  expec <- expected.msm(x, t, timezero=timezero, initstates=initstates, covariates=covariates, misccovariates=misccovariates,
                        piecewise.times=piecewise.times, piecewise.covariates=piecewise.covariates, risk=obs$risk, subset=subset, ci="none")[[2]]
  states <- seq(length=x$qmodel$nstates)
  S <- length(states)
  ncols <- ceiling(sqrt(S))
  nrows <- if (floor(sqrt(S))^2 < S && S <= floor(sqrt(S))*ceiling(sqrt(S))) floor(sqrt(S)) else ceiling(sqrt(S))
  par(mfrow=c(nrows, ncols))
  for (i in states) {
    plot(t, obs$obsperc[,i], type="l", ylim=c(0, 100), xlab=xlab, ylab=ylab, lwd=lwd.obs, lty=lty.obs, col=col.obs,
         main=rownames(x$qmodel$qmatrix)[i],...)
    lines(t, expec[,i], lwd=lwd.exp, lty=lty.exp, col=col.exp)
  }
  if (!is.numeric(legend.pos) || length(legend.pos) != 2)
    legend.pos <- c(0.4*maxtime, 40)
  legend(x=legend.pos[1], y=legend.pos[2], legend=c("Observed","Expected"), lty=c(lty.obs,lty.exp), lwd=c(lwd.obs,lwd.exp), col=c(col.obs,col.exp))
  invisible()
}



plot_idm <- function(fit, formula01, formula02, formula12, data){

  formula <- list(formula01, formula02, formula12)

  formula <- lapply(formula, function(f) parse_formula(f, data = data))

  data01 <- handle_state01(formula[[1]], formula[[2]], data)
  data02 <- handle_state02(formula[[1]], formula[[2]], data)

  formula01 <- formula[[1]]
  formula02 <- formula[[2]]
  formula12 <- formula[[3]]


  km01 <- survival::survfit(as.formula(paste("Surv(time = time01, event = delta01)~", formula01$rhs , collapse = "" ) ), data = data01)

  p01 <- survminer::ggsurvplot(km01, risk.table = TRUE) + ggtitle(expression('S'['01']*'(t)') )

  km02 <- survival::survfit(as.formula(paste("Surv(time = time02, event = delta02)~", formula02$rhs , collapse = "" ) ), data = data02)

  p02 <- survminer::ggsurvplot(km02, risk.table = TRUE) + ggtitle( expression('S'['02']*'(t)'))
  data12 = data[data[[formula[[1]]$dvar ]] == 1,]
  km12 <- survival::survfit(as.formula(paste( c(formula12$lhs, "~", formula12$rhs ), collapse = "" ) ) , data =  data12)

  p12 <- survminer::ggsurvplot(km12, risk.table = TRUE) + ggtitle( expression('S'['12']*'(t)') ) +xlab("Time since non-terminal event")

  est_fit <- survminer::arrange_ggsurvplots(list(p01, p02, p12), ncol = 3, nrow = 1 )


  pars = list(
    pp_lambda01 = exp( rstan::extract(fit,'alpha01')$alpha01 ),
    pp_gamma01 = rstan::extract(fit,'aux01')$aux01,
    pp_beta01 = rstan::extract(fit,'beta01')$beta01,
    pp_lambda02 = exp( rstan::extract(fit,'alpha02')$alpha02 ),
    pp_gamma02 = rstan::extract(fit,'aux02')$aux02,
    pp_beta02 = rstan::extract(fit,'beta02')$beta02,
    pp_lambda12 = exp( rstan::extract(fit,'alpha12')$alpha12 ),
    pp_gamma12 = rstan::extract(fit,'aux12')$aux12,
    pp_beta12 = rstan::extract(fit,'beta12')$beta12
  )


  pp_newdata <-
    purrr::pmap(pars, function( pp_lambda01_i,pp_gamma01_i,pp_beta01_i,
                                pp_lambda02_i,pp_gamma02_i,pp_beta02_i,
                                pp_lambda12_i,pp_gamma12_i,pp_beta12_i){

      pp_beta01_i <- c(pp_beta01_i)
      names(pp_beta01_i) <- colnames(covs)[-1]
      pp_beta02_i <- c(pp_beta02_i)
      names(pp_beta02_i) <- colnames(covs)[-1]
      pp_beta12_i <- c(pp_beta12_i)
      names(pp_beta12_i) <- colnames(covs)[-1]
      rsimid(
        dist01 = "weibull",
        dist02 = "weibull",
        dist12 = "weibull",
        betas01 = pp_beta01_i,
        betas02 = pp_beta02_i,
        betas12 = pp_beta12_i,
        lambdas01 = pp_lambda01_i,
        lambdas02 = pp_lambda02_i,
        lambdas12 = pp_lambda12_i,
        gammas01 = pp_gamma01_i,
        gammas02 = pp_gamma02_i,
        gammas12 = pp_gamma12_i,
        x = covs
      )
    }

    )

}



handle_state01 <- function(formula01, formula02, data){
  N <- nrow(data) # number of individuals
  ids <- seq(N)
  dataR <- data[ ,c(formula01$allvars)]
  dataD <- data[ ,c(formula02$allvars)]

  time01 <- pmin(dataR[, formula01$tvar_end], dataD[, formula02$tvar_end])

  yesR <- dataR[, formula01$tvar_end] <= dataD[, formula02$tvar_end]

  delta01 <- rep(NA, N)
  delta01[yesR] <- dataR[yesR, formula01$dvar]
  delta01[!yesR] <- 0

  data01 <- data
  data01$time01 <- time01
  data01$delta01 <- delta01

  data01
}

handle_state02 <- function(formula01, formula02, data){
  N <- nrow(data) # number of individuals
  ids <- seq(N)
  dataR <- data[ ,c(formula01$allvars)]
  dataD <- data[ ,c(formula02$allvars)]

  time02 <- pmin(dataR[, formula01$tvar_end], dataD[, formula02$tvar_end])

  yesD <- dataD[, formula02$tvar_end] <= dataR[, formula01$tvar_end]

  delta02 <- rep(NA, N)
  delta02[yesD] <- dataD[yesD, formula02$dvar]
  delta02[!yesD] <- 0

  data02 <- data
  data02$time02 <- time02
  data02$delta02 <- delta02

  data02
}



