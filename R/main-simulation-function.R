sim_spacetime <- function(nsims = 1){

}

#' Generate a spatiotemporal variable
#'
#' @param map.dims map dimensions, as in \code{generate_spcov()}.
#' @param ntime number of time points to simulate.
#' @param break.df a breakpoint dataframe, as created with \code{break_frame()}.
#' @param seasonality an optional named list of seasonality parameters, which
#' is passed to \code{season()}.
#' @param covar_FUN a covariance function, as in \code{generate_spcov()}.
#' @param covar.pars a list of covariance parameters, as in \code{generate_spcov()}.
#' @param sd the standard deviation around the underlying trend and seasonality.
#' @param arima.mod a list of ARiMA parameters, which is passed to the \code{model}
#' argument of \code{sparima_sim()}.
#' @param times an optional list of time points. When this argument is anything
#' other than \code{NA}, \code{ntime = length(times)}.
#'
#' @return a 'spacebreak' list object with the following elements
#'
#' \describe{
#'   \item{Y}{A matrix containing the spatiotemporal reponse variable, with
#'   rows corresponding to map pixels, and columns corresponding to time
#'   points. For each pixel, \code{Y = f.t + s.t + sptemp}.}
#'   \item{f.t}{a vector of length \code{ntime}, containing values of the
#'   piecewise timeseries function, evaluated at each time point.}
#'   \item{s.t}{an optional vector of length \code{ntime}, containing values
#'   evaluated from the \code{season()} function for each time point.}
#'   \item{time}{a vector of time points.}
#'   \item{space}{a list, as output by \code{generate_spcov()}.}
#'   \item{spatiotemp}{a matrix containing the spatiotemporal autocorrelation
#'   component of Y.}
#' }
#'
#' @export
#'
#' @examples
#'
#' ## No seasonality
#' bdf <- break_frame(nodes = c(0, 20, 21, 50), slopes = c(.1, -1, .1))
#' spacebreak <- generate_spacebreak(break.df = bdf, sd = 1)
#'
#' plot.spacebreak(x = spacebreak)
#'
#' ## With seasonality (every 4 time steps)
#' spacebreak2 <- generate_spacebreak(break.df = bdf, seasonality = list(amp = .5, per = 4), sd = 1)
#'
#' plot.spacebreak(x = spacebreak2)
#'
generate_spacebreak <- function(map.dims = c(3, 3), ntime = 50, break.df,
                               seasonality = NULL,
                               covar_FUN = remotePARTS::covar_exp,
                               covar.pars = list(range = .1),
                               sd = .2,
                               arima.mod = list(ar = .2),
                               times = NULL,
                               burn.in = 20){

  ## build time series, if not provided
  if(is.null(times)){
    times = seq(0, (ntime - 1), by = 1)
  } else {
    ntime = length(times)
  }

  # ## start with burn in, if needed
  # if(!is.na(burn.in) & burn.in > 0){
  #
  # }

  ## build piecewise function
  f.t = breakwise(times, break.df)
  ## create seasonality, if needed
  if(!is.null(seasonality)){
    s.t = do.call(season, append(list(times), seasonality))
  } else {
    s.t = 0
  }
  ## generate spatial backbone
  space <- generate_spcov(map.dims = map.dims, covar_FUN = covar_FUN,
                          covar.pars = covar.pars, sd = sd)
  ## generate spatiotemporal autocorrelation
  sptemp <- sparima_sim(ntime = ntime, covar = space$covar.mat,
                        model = arima.mod, burn.in = burn.in)


  ## Build response
  Y = t(apply(X = sptemp, MARGIN = 1, FUN = function(x)x + s.t + f.t))

  ## build output
  if(!is.null(seasonality)){
    seas = s.t
  } else {
    seas = s.t
  }


  out <- list(Y = Y, f.t = f.t, s.t = seas,
              time = times,
              space = space, spatiotemp = sptemp)
  class(out) <- append("spacebreak", class(out))

  ## return
  out
}


#' plot a spacebreak object
#'
#' @param x the spacebreak object
#' @export
#'
#' @return
plot.spacebreak <- function(x, type = "ts"){
  if(type == "ts"){
    matplot(y = t(x$Y), x = x$time, type = "l", xlab = "time", ylab = "Y")
    lines(x$f.t ~ x$time, col = "red", lty = 1, lwd = 2)
  } else if(type == "map"){
    cat("no 'map' method yet")
  }
}

# ## piecewise function
# breakframe = break_frame(nodes = c(0, 15, 19, 30),
#                          slopes = c(.1, -.25, .1),
#                          initial = 0)
# f.t <- breakwise(0:30, breakframe)
#
# ## seasonality
# s.t <- season(0:30, amp = .5, per = 4)
#
# ### visualize time trends
# plot(f.t + s.t, type = "l")
# lines(x = f.t, col = "red")
# abline(v = c(16, 20), lty = 2, col = "blue")
#
# ## spatial dimensions
# spcov <- generate_spcov(map.dims = c(3,3), covar.pars = list(range = .1),
#                         sd = .2)
#
# ## spatiotemporal autocorrelation
# sptemp <- sparima_sim(ntime = 31, covar = spcov$covar.mat,
#                       model = list(ar = .5))
#
# ### visualize temporal autocorrelation at each pixel
# matplot(t(sptemp), type = "l")
#
# ## Build the response
# Y.t <- t(apply(X = sptemp, MARGIN = 1, FUN = function(x){x + s.t + f.t}))
#
# ### Visualize the response and piecewise trend
# matplot(t(Y.t), type = "l");
# lines(x = f.t, col = "red", lty = 1, lwd = 2)
