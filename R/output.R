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

# plot method for stanidm ----------------------------------------------
#' Plot stanidm
#'
#' @rdname plot.stanidm
#' @export
#' @templateVar cigeomArg ci_geom_args
#' @param prob A scalar between 0 and 1 specifying the width to
plot.stanidm <- function(x, plotfun = "basehaz", pars = NULL,
                          regex_pars = NULL, ..., prob = 0.95,
                          limits = c("ci", "none"),
                          ci_geom_args = NULL,
                         labels = "auto") {

  validate_stanidm_object(x)

  limits <- match.arg(limits)

  if (plotfun %in% c("basehaz", "tde")) {

    stanpars <- extract_pars(x)

    secondpars <- lapply(seq_along(x$basehaz), function(i){
      has_intercept <- check_for_intercept(x$basehaz[[i]])
      t_min <- min(x$entrytime[[i]])
      t_max <- max(x$eventtime[[i]])
      times <- seq(t_min, t_max, by = (t_max - t_min) / 200)
      nlist(has_intercept, t_min, t_max, times)
    })


    if (plotfun == "basehaz") {
      if (!is.null(pars))
        warning2("'pars' is ignored when plotting the baseline hazard.")
      if (!is.null(regex_pars))
        warning2("'regex_pars' is ignored when plotting the baseline hazard.")

      plotdat <- lapply(seq_along(x$basehaz), function(i){

        args <- nlist(times     = secondpars[[i]]$times,
                      basehaz   = get_basehaz(x)[[i]],
                      aux       = stanpars[[i]]$aux,
                      intercept = stanpars[[i]]$alpha)
        basehaz <- do.call(evaluate_basehaz, args)
        basehaz <- median_and_bounds(basehaz, prob, na.rm = TRUE)
        plotdat <- data.frame(times = secondpars[[i]]$times, basehaz)

      })
      ylab <- "Baseline hazard rate"
      xlab <- "Time"
    } else if (plotfun == "tde") {
      stop2("tde not implemented")
    }
    geom_defs <- list(color = "black")  # default plot args
    geom_args <- set_geom_args(geom_defs, ...)
    geom_ylab <- ggplot2::ylab(ylab)
    geom_xlab <- ggplot2::xlab(xlab)
    geom_maps <- list(aes_string(x = "times", y = "med"), method = "loess", se = FALSE)


    plots <- lapply(seq_along(x$basehaz), function(i){
      geom_base <- ggplot(plotdat[[i]]) + geom_ylab + geom_xlab + ggplot2::theme_bw()
      geom_plot <- geom_base + do.call(ggplot2::geom_smooth, c(geom_maps, geom_args))
      if (limits == "ci") {
        lim_defs <- list(alpha = 0.3) # default plot args for ci
        lim_args <- c(defaults = list(lim_defs), ci_geom_args)
        lim_args <- do.call("set_geom_args", lim_args)
        lim_maps <- list(mapping = aes_string(x = "times", ymin = "lb", ymax = "ub"))
        lim_tmp  <- geom_base +
          ggplot2::stat_smooth(aes_string(x = "times", y = "lb"), method = "loess") +
          ggplot2::stat_smooth(aes_string(x = "times", y = "ub"), method = "loess")
        lim_build<- ggplot2::ggplot_build(lim_tmp)
        lim_data <- list(data = data.frame(times = lim_build$data[[1]]$x,
                                           lb    = lim_build$data[[1]]$y,
                                           ub    = lim_build$data[[2]]$y))
        lim_plot <- do.call(ggplot2::geom_ribbon, c(lim_data, lim_maps, lim_args))
      } else {
        lim_plot <- NULL
      }
      geom_plot + lim_plot
    })

    out <- cowplot::plot_grid(
      plotlist = plots,
      axis = "rlbt",
      labels = labels,
      label_size = 10)
    return(out)
  }
  NextMethod("plot")
}


# Set plotting defaults
set_geom_args <- function(defaults, ...) {
  dots <- list(...)
  if (!length(dots))
    return(defaults)
  dot_names <- names(dots)
  def_names <- names(defaults)
  for (j in seq_along(def_names)) {
    if (def_names[j] %in% dot_names)
      defaults[[j]] <- dots[[def_names[j]]]
  }
  extras <- setdiff(dot_names, def_names)
  if (length(extras)) {
    for (j in seq_along(extras))
      defaults[[extras[j]]] <- dots[[extras[j]]]
  }
  return(defaults)
}


# internal for plot.stanreg ----------------------------------------------

# Prepare argument list to pass to plotting function
#
# @param x stanreg object
# @param pars, regex_pars user specified pars and regex_pars arguments (can be
#   missing)
# @param ...  additional arguments to pass to the plotting function
# @param plotfun User's 'plotfun' argument
set_plotting_args <- function(x, pars = NULL, regex_pars = NULL, ...,
                              plotfun = character()) {

  plotfun <- mcmc_function_name(plotfun)
  if (!used.sampling(x))
    validate_plotfun_for_opt_or_vb(plotfun)

  .plotfun_is_type <- function(patt) {
    grepl(pattern = paste0("_", patt), x = plotfun, fixed = TRUE)
  }

  if (.plotfun_is_type("nuts")) {
    nuts_stuff <- list(x = bayesplot::nuts_params(x), ...)
    if (!.plotfun_is_type("energy"))
      nuts_stuff[["lp"]] <- bayesplot::log_posterior(x)
    return(nuts_stuff)
  }
  if (.plotfun_is_type("rhat")) {
    rhat <- bayesplot::rhat(x, pars = pars, regex_pars = regex_pars)
    return(list(rhat = rhat, ...))
  }
  if (.plotfun_is_type("neff")) {
    ratio <- bayesplot::neff_ratio(x, pars = pars, regex_pars = regex_pars)
    return(list(ratio = ratio, ...))
  }
  if (!is.null(pars) || !is.null(regex_pars)) {
    pars <- collect_pars(x, pars, regex_pars)
    pars <- allow_special_parnames(x, pars)
  }

  if (!used.sampling(x)) {
    if (!length(pars))
      pars <- NULL
    return(list(x = as.matrix(x, pars = pars), ...))
  }

  if (needs_chains(plotfun))
    list(x = as.array(x, pars = pars, regex_pars = regex_pars), ...)
  else
    list(x = as.matrix(x, pars = pars, regex_pars = regex_pars), ...)
}

mcmc_function_name <- function(fun) {
  # to keep backwards compatibility convert old function names
  if (fun == "scat") {
    fun <- "scatter"
  } else if (fun == "ess") {
    fun <- "neff"
  } else if (fun == "ac") {
    fun <- "acf"
  } else if (fun %in% c("diag", "stan_diag")) {
    stop(
      "For NUTS diagnostics, instead of 'stan_diag', ",
      "please specify the name of one of the functions listed at ",
      "help('NUTS', 'bayesplot')",
      call. = FALSE
    )
  }

  if (identical(substr(fun, 1, 4), "ppc_"))
    stop(
      "For 'ppc_' functions use the 'pp_check' ",
      "method instead of 'plot'.",
      call. = FALSE
    )

  if (!identical(substr(fun, 1, 5), "mcmc_"))
    fun <- paste0("mcmc_", fun)

  if (!fun %in% bayesplot::available_mcmc())
    stop(
      fun, " is not a valid MCMC function name.",
      " Use bayesplot::available_mcmc() for a list of available MCMC functions."
    )

  return(fun)
}

# check if a plotting function requires multiple chains
needs_chains <- function(x) {
  nms <- paste0("mcmc_",
                c(
                  "trace",
                  "trace_highlight",
                  "acf",
                  "acf_bar",
                  "hist_by_chain",
                  "dens_overlay",
                  "violin",
                  "combo"
                )
  )
  mcmc_function_name(x) %in% nms
}

# Select the correct plotting function
# @param plotfun user specified plotfun argument (can be missing)
set_plotting_fun <- function(plotfun = NULL) {
  if (is.null(plotfun))
    return("mcmc_intervals")
  if (!is.character(plotfun))
    stop("'plotfun' should be a string.", call. = FALSE)

  plotfun <- mcmc_function_name(plotfun)
  fun <- try(get(plotfun, pos = asNamespace("bayesplot"), mode = "function"),
             silent = TRUE)
  if (!inherits(fun, "try-error"))
    return(fun)

  stop(
    "Plotting function ",  plotfun, " not found. ",
    "A valid plotting function is any function from the ",
    "'bayesplot' package beginning with the prefix 'mcmc_'.",
    call. = FALSE
  )
}

# check if plotfun is ok to use with vb or optimization
validate_plotfun_for_opt_or_vb <- function(plotfun) {
  plotfun <- mcmc_function_name(plotfun)
  if (needs_chains(plotfun) ||
      grepl("_rhat|_neff|_nuts_", plotfun))
    STOP_sampling_only(plotfun)
}

# pairs method ------------------------------------------------------------

#' Pairs method for stanidm objects
#'
#' Interface to \pkg{bayesplot}'s \code{\link[bayesplot]{mcmc_pairs}} function
#' for use with \pkg{rstanarm} models. Be careful not to specify too
#' many parameters to include or the plot will be both hard to read and slow to
#' render.
#'
#' @method pairs stanidm
#' @export
#' @importFrom bayesplot pairs_style_np pairs_condition
pairs.stanidm <-
  function(x,
           pars = NULL,
           regex_pars = NULL,
           condition = pairs_condition(nuts = "accept_stat__"),
           ...) {

    if (!used.sampling(x))
      STOP_sampling_only("pairs")

    dots <- list(...)
    ignored_args <- c("np", "lp", "max_treedepth")
    specified <- ignored_args %in% names(dots)
    if (any(specified)) {
      warning(
        "The following arguments were ignored because they are ",
        "specified automatically by rstanarm: ",
        paste(sQuote(ignored_args[specified]), collapse = ", ")
      )
    }

    posterior <- as.array.stanreg(x, pars = pars, regex_pars = regex_pars)
    if (is.null(pars) && is.null(regex_pars)) {
      # include log-posterior by default
      lp_arr <- as.array.stanreg(x, pars = "log-posterior")
      dd <- dim(posterior)
      dn <- dimnames(posterior)
      dd[3] <- dd[3] + 1
      dn$parameters <- c(dn$parameters, "log-posterior")
      tmp <- array(NA, dim = dd, dimnames = dn)
      tmp[,, 1:(dd[3] - 1)] <- posterior
      tmp[,, dd[3]] <- lp_arr
      posterior <- tmp
    }
    posterior <- round(posterior, digits = 12)

    bayesplot::mcmc_pairs(
      x = posterior,
      np = bayesplot::nuts_params(x),
      lp = bayesplot::log_posterior(x),
      max_treedepth = .max_treedepth(x),
      condition = condition,
      ...
    )

  }


# pairs method ------------------------------------------------------------

#' Pairs method for stanreg objects
#'
#' Interface to \pkg{bayesplot}'s \code{\link[bayesplot]{mcmc_pairs}} function
#' for use with \pkg{rstanarm} models. Be careful not to specify too
#' many parameters to include or the plot will be both hard to read and slow to
#' render.
#'
#' @method pairs stanreg
#' @export
#' @importFrom bayesplot pairs_style_np pairs_condition
#' @export pairs_style_np pairs_condition
#' @aliases pairs_style_np pairs_condition
#'s
#' @templateVar stanregArg x
#' @param pars An optional character vetor of parameter names. All parameters
#'   are included by default, but for models with more than just a few
#'   parameters it may be far too many to visualize on a small computer screen
#'   and also may require substantial computing time.
#' @param condition Same as the \code{condition} argument to
#'   \code{\link[bayesplot]{mcmc_pairs}} except the \emph{default is different}
#'   for \pkg{rstanarm} models. By default, the \code{mcmc_pairs} function in
#'   the \pkg{bayesplot} package plots some of the Markov chains (half, in the
#'   case of an even number of chains) in the panels above the diagonal and the
#'   other half in the panels below the diagonal. However since we know that
#'   \pkg{rstanarm} models were fit using Stan (which \pkg{bayesplot} doesn't
#'   assume) we can make the default more useful by splitting the draws
#'   according to the \code{accept_stat__} diagnostic. The plots below the
#'   diagonal will contain realizations that are below the median
#'   \code{accept_stat__} and the plots above the diagonal will contain
#'   realizations that are above the median \code{accept_stat__}. To change this
#'   behavior see the documentation of the \code{condition} argument at
#'   \code{\link[bayesplot]{mcmc_pairs}}.
#' @param ... Optional arguments passed to \code{\link[bayesplot]{mcmc_pairs}}.
#'   The \code{np}, \code{lp}, and \code{max_treedepth} arguments to
#'   \code{mcmc_pairs} are handled automatically by \pkg{rstanarm} and do not
#'   need to be specified by the user in \code{...}. The arguments that can be
#'   specified in \code{...} include \code{transformations}, \code{diag_fun},
#'   \code{off_diag_fun}, \code{diag_args}, \code{off_diag_args},
#'   and \code{np_style}. These arguments are
#'   documented thoroughly on the help page for
#'   \code{\link[bayesplot]{mcmc_pairs}}.
#'
#'
#' @examples
#' \donttest{
#' if (!exists("example_model")) example(example_model)
#'
#' bayesplot::color_scheme_set("purple")
#'
#' # see 'condition' argument above for details on the plots below and
#' # above the diagonal. default is to split by accept_stat__.
#' pairs(example_model, pars = c("(Intercept)", "log-posterior"))
#'
#' pairs(
#'   example_model,
#'   regex_pars = "herd:[2,7,9]",
#'   diag_fun = "dens",
#'   off_diag_fun = "hex"
#' )
#' }
#'
#' \donttest{
#' # for demonstration purposes, intentionally fit a model that
#' # will (almost certainly) have some divergences
#' fit <- stan_glm(
#'   mpg ~ ., data = mtcars,
#'   iter = 1000,
#'   # this combo of prior and adapt_delta should lead to some divergences
#'   prior = hs(),
#'   adapt_delta = 0.9
#' )
#'
#' pairs(fit, pars = c("wt", "sigma", "log-posterior"))
#'
#' pairs(
#'   fit,
#'   pars = c("wt", "sigma", "log-posterior"),
#'   transformations = list(sigma = "log"), # show log(sigma) instead of sigma
#'   off_diag_fun = "hex" # use hexagonal heatmaps instead of scatterplots
#' )
#'
#'
#' bayesplot::color_scheme_set("brightblue")
#' pairs(
#'   fit,
#'   pars = c("(Intercept)", "wt", "sigma", "log-posterior"),
#'   transformations = list(sigma = "log"),
#'   off_diag_args = list(size = 3/4, alpha = 1/3), # size and transparency of scatterplot points
#'   np_style = pairs_style_np(div_color = "black", div_shape = 2) # color and shape of the divergences
#' )
#'
#' # Using the condition argument to show divergences above the diagonal
#' pairs(
#'   fit,
#'   pars = c("(Intercept)", "wt", "log-posterior"),
#'   condition = pairs_condition(nuts = "divergent__")
#' )
#'
#' }
#'
pairs.stanreg <-
  function(x,
           pars = NULL,
           regex_pars = NULL,
           condition = pairs_condition(nuts = "accept_stat__"),
           ...) {

    if (!used.sampling(x))
      STOP_sampling_only("pairs")

    dots <- list(...)
    ignored_args <- c("np", "lp", "max_treedepth")
    specified <- ignored_args %in% names(dots)
    if (any(specified)) {
      warning(
        "The following arguments were ignored because they are ",
        "specified automatically by rstanarm: ",
        paste(sQuote(ignored_args[specified]), collapse = ", ")
      )
    }

    posterior <- as.array.stanreg(x, pars = pars, regex_pars = regex_pars)
    if (is.null(pars) && is.null(regex_pars)) {
      # include log-posterior by default
      lp_arr <- as.array.stanreg(x, pars = "log-posterior")
      dd <- dim(posterior)
      dn <- dimnames(posterior)
      dd[3] <- dd[3] + 1
      dn$parameters <- c(dn$parameters, "log-posterior")
      tmp <- array(NA, dim = dd, dimnames = dn)
      tmp[,, 1:(dd[3] - 1)] <- posterior
      tmp[,, dd[3]] <- lp_arr
      posterior <- tmp
    }
    posterior <- round(posterior, digits = 12)

    bayesplot::mcmc_pairs(
      x = posterior,
      np = bayesplot::nuts_params(x),
      lp = bayesplot::log_posterior(x),
      max_treedepth = .max_treedepth(x),
      condition = condition,
      ...
    )

  }


# internal for pairs.stanreg ----------------------------------------------

# @param x stanreg object
.max_treedepth <- function(x) {
  control <- x$stanfit@stan_args[[1]]$control
  if (is.null(control)) {
    max_td <- 10
  } else {
    max_td <- control$max_treedepth
    if (is.null(max_td))
      max_td <- 10
  }
  return(max_td)
}


plot2.stanidm <- function(fit){

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



