# spatial_functions.R

#' Generate pixels with spatial correlations for a grid-based map
#'
#' @param map.dims a dimension vector of length 2, containing the rows and
#' columns, respectively of the map. \code{map.dims} are coerced to integers
#' with \code{as.integer(map.dims)}. The number of pixels (n) in the resulting
#' map is equal to \code{map.dims[1]*map.dims[2]}.
#' @param unit.map logical: should a unit map be used (range from 0 to 1)?
#' @param covar_FUN function to generate distance-based correlations from.
#' \code{cavar_FUN} should accept an nxn distance matrix (\code{D.mat}) and
#' return a nxn covariance matrix.
#' @param covar.pars a list of parameters used by \code{covar_FUN}.
#' @param sd the standard deviation of the system. The covariance matrix
#' that results from \code{covar_FUN(D.mat, ...)} is scaled by
#' \code{sqrt(sd)}.
#'
#' @return a list with the following elements
#'
#' \describe{
#'   \item{coord.df}{a data.frame with n rows and columns corresponding
#'   to pixel identifiers (\code{pixel}), the pixel's x-coordinate
#'   (\code{coords.x}) and y-coordinate (\code{coords.y})}
#'   \item{D.mat}{the nxn distance matrix, containing pairwise distances
#'   among all pixels}
#'   \item{covar.mat}{the nxn covariance matrix, containing pairwise
#'   covariances among all pixels}
#'
#' }
#'
#' @export
#'
#' @examples
#' results = generate_spcov(c(5, 5))
#' results$coord.df
#' results$covar.mat
#' results$D.mat
#'
#' ## smaller sd
#' generate_spcov(c(5, 5), sd = .5)$covar.mat
#'
#' ## unequal map dimensions
#' generate_spcov(c(2, 5))$coord.df
#'
#' ## plot distance-based covariance function
#' plot(results$covar.mat ~ results$D.mat, xlab = "distance",
#'      ylab = "covariance")
generate_spcov <- function(map.dims = c(5, 5), unit.map = TRUE,
                           covar_FUN = remotePARTS::covar_exp,
                           covar.pars = list(range = 1),
                           sd = 1){
  map.dims = as.integer(map.dims)
  npix <- unname(map.dims[1]*map.dims[2]) #total pixels
  coords <- as.matrix(expand.grid(x = 1:map.dims[2], y = 1:map.dims[1])) #coordinates

  if(unit.map){ #rescale coordinates to between 0 and 1
    coords = scales::rescale(coords, to = c(0, 1))
  }

  D.mat <- as.matrix(dist(coords)) #calculate distances
  covar.fun = match.fun(covar_FUN) #calculate covariance
  cov.mat <- sqrt(sd)*do.call(covar.fun, append(list(D.mat), covar.pars)) #spatial covariance
  ## return
  list(coord.df = data.frame(pixel = 1:npix, coords = coords),
       D.mat = D.mat,
       covar.mat = cov.mat)
}

#' Generate spatially-autocorrelated random variables
#'
#' @param n number of samples per location
#' @param covar a matrix containing pairwise covariances among
#' locations.
#' @param fast logical: should the faster method, whose realizations
#' observe lower fidelity to the covariance function than when
#' \code{fast = FALSE}.
#'
#' @return a matrix with rows corresponding to locations and columns
#' corresponding to random samples taken at that location
#'
#' @export
#'
#' @examples
#' cor.mat <- matrix(c(1.0, 0.3, 0.1,
#'                     0.3, 1.0, 0.0,
#'                     0.1, 0.0, 1.0), nrow = 3)
#' rnorm_spcor(5, cor.mat)
#'
#' ## plot sample covariance vs. true covariance
#' spcov = generate_spcov(c(10, 10), sd = 1, unit.map = TRUE)
#' D = spcov$D.mat
#' C = spcov$covar.mat
#' x <- rnorm_spcor(100, C)
#' ## sample covariance vs. distance
#' plot(x = D, y = cor(t(x)), xlab = "distance", ylab = "covariance")
#' ## add line for covariance function
#' curve(remotePARTS::covar_exp(x, 1), from = 0, to = max(D),
#'       add = TRUE, col = "red")
#'
rnorm_spcor <- function(n, covar, fast = FALSE){

  if(!fast){
    vals <- t(mvtnorm::rmvnorm(n, sigma = covar)) # more accurate but slower
  } else {
    if(n == 1){ #spatial errors
      vals <- covar %*% matrix(rnorm(1), nrow = nrow(covar))
    } else if(n > 1) {
      vals <- covar %*% t(replicate(nrow(covar), rnorm(n)))
    } else {
      stop("invalid n")
    }
  }
  ## return
  vals
}
