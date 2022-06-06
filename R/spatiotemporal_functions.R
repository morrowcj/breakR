
#' Simulate an ARIMA model with spatially correlated errors
#'
#' @param ntime number of time points
#' @param covar matrix containing pairwise covariances among pixels
#' @param fast logical: should fast method be used? see \code{rnorm_spcor}
#' @param innov optional matrix of innovations to use in place of \code{ntime},
#' \code{covar}, and \code{fast}.
#' @param model a named list containing model specification for \code{arima.sim()}.
#' This list should NOT include an 'innov' element.
#' @param ... additional arguments passed to \code{arima.sim()}
#'
#' @return
#' @export
#'
#' @examples
#' sparima_sim(20, cor.mat)
#' sparima_sim(innov = rnorm_spcor(20, cor.mat))
#' ## specify burn-in period (from arima.sim)
#' sparima_sim(20, cor.mat, n.start = 20)
#' ## single time point (after 100)
#' sparima_sim(1, cor.mat, n.start = 100)
#'
#' ## compare observed covariance to true covariance
#' spcor = generate_spcov(c(10, 10), sd = 1, unit.map = TRUE)
#' D = spcov$D.mat
#' C = spcov$covar.mat
#' x <- sparima_sim(50, C, model = list(ar = .8), n.start = 100)
#' ## sample covariance vs. distance
#' plot(x = D, y = cor(t(x)), xlab = "distance", ylab = "covariance")
#' ## add line for covariance function
#' curve(remotePARTS::covar_exp(x, 1), from = 0, to = max(D),
#'       add = T, col = "red")
sparima_sim <- function(ntime, covar, fast = FALSE, model = list(ar = .2),
                        burn.in = 20, innov = NULL, ...){
  if(!is.null(innov)){
    ## use pre-calculat♣ed innovations
    innov = as.matrix(innov)
    ntime = ncol(innov)
  } else {
    ## calculate spatially-autocorrelated innovations
    innov = rnorm_spcor(ntime, covar, fast)
  }
  ## generate spatiotemporally autocorrelated random variable
  if(ntime > 1) {
    sptmp = t(apply(X = innov, MARGIN = 1, function(e){
      arima.sim(n = ntime, model = model, innov = e, n.start = burn.in, ...)
    }))
  } else if(ntime == 1){ # (for spatial-only case)
    sptmp = matrix(apply(X = innov, MARGIN = 1, function(e){
      arima.sim(n = ntime, model = model, innov = e,  n.start = burn.in, ...)
    }), ncol = ntime)
  }
  ## return
  sptmp
}
