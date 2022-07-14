# example functions

#' Melt a spatiotemporal matrix into a long-form data frame
#'
#' @param Y a spatiotemporal matrix, as results from \code{generate_spacebreak()}.
#'
#' @return a data frame with three columns: indicating the \code{pixel} ID,
#' the \code{time} point, and the data \code{value}.
#'
#' @export
melt_spacebreak <- function(Y){
  reshape2::melt(Y, varnames = list(c("pixel"), c("time")))
}

#' Example of a fit_FUN, passed to simulate_spacebreak
#'
#' @param Y a matrix whose rows contain pixels and whose columns contain
#' timepoints.
#' @param ... additional arguments
#'
#' @export
example_fitfun <- function(Y, ...){
  # convert Y (matrix) into a long-form dataframe
  data = melt_spacebreak(Y)

  ## fit the models of interest
  model.lm <- stats::lm(value ~ time, data = data) # linear
  model.lmer <- lme4::lmer(value ~ (1|time), data = data) #random effects
  ## any other model that you want goes here ...

  # return the model outputs, even if they differ.
  output <- list(lm = model.lm, lmer = model.lmer)
}
