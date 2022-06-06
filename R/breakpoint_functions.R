# breakpoint_functions.R

#' Build a breakpoint dataframe, which describes the line segments of a linear piece-wise
#' time series function
#'
#' @param nodes a chronological vector of time points at which breakpoint nodes occur,
#' including the start and end of the full time series.
#' @param slopes a vector of slopes for the line segments between each node.
#' There needs to be one fewer slope than nodes.
#' @param initial the initial numeric value of the time series data at the first time
#' step (i.e., the first value in \code{nodes})
#' @param slopeseries an optional time series of linear slopes. if \code{slopeseries}
#' is not \code{NULL} (the default), then the arguments \code{nodes} and \code{slopes}
#' are ignored entirely.
#'
#' @return a dataframe object of class \code{break_df}, with the following
#' columns describing each line segment in the piece-wise time series function.
#'
#' \describe{
#'   \item{segment}{a numeric identifier of the line segment}
#'   \item{left}{the time at which the left node of the segment occurs}
#'   \item{right}{the time at which the right node of the segment occurs}
#'   \item{start}{the value at the left node of the segment}
#'   \item{end}{the value at the right node of the segment}
#'   \item{slope}{the slope of the segment}
#' }
#'
#' @export
#'
#' @examples
#' break_frame(nodes = c(0, 5, 9), slopes = c(.2, -.2))
#' break_frame(c(0, 5, 9), c(0.2, -.2), initial = 1)
#' break_frame(c(1, 6, 10), c(.2, -.2), 1)
#' break_frame(c(0, .2, .5, 1), c(1, 0, -1))
#' ## using a time series of slopes
#' break_frame(slopeseries = c(rep(.5, 2), rep(0, 2), rep(-.5, 2)))
break_frame <- function(nodes, slopes,
                        initial = 0,
                        slopeseries = NULL){

  if (!is.null(slopeseries)) {
    tmp = get_breaks(slopeseries)
    nodes = tmp$nodes
    slopes = tmp$slopes
  }

  ## check for appropriate dimensions
  stopifnot(length(slopes) == length(nodes) - 1)
  ## initialize
  start = initial
  break.df = data.frame(segment = 1:length(slopes), left = NA, right = NA,
                        start = NA, end = NA, slope = NA)
  for(i in 1:length(slopes)){
  ## set/update values in each segment
    left = nodes[i]
    right = nodes[i + 1]
    slope = slopes[i]
    end = start + slope * (right - left)
    break.df[i, ] = c(segment = i, left = left, right = right,
                      start = start, end = end, slope = slope)
    start = end
  }
  ## append a new class string
  class(break.df) <- append(class(break.df), 'break_df')

  ## return
  break.df
}

#' Get breakpoints and nodes from a time series of slopes
#'
#' @param slope.vec a vector of linear slopes for each time point of a
#' time series
#'
#' @return a list containing vectors of \code{nodes} and \code{slopes}
#' for each line segment
#'
#' @examples
#' (slope_vector = c(rep(0, 3), rep(.2, 4), rep(-.1, 1)))
#' (breakR:::get_breaks(slope.vec = slope_vector))
get_breaks <- function(slope.vec){
  ## starting values
  slp = slope.vec[1]
  nodes = 0
  slopes = slp
  ## loop
  for (i in seq_len(length(slope.vec))){
    ## append the values, if slope changed
    if(slope.vec[i] != slp){
      nodes = c(nodes, i)
      slopes = c(slopes, slope.vec[i])
    }

    ## add the last time point as a node
    if(i == length(slope.vec)){
      nodes = c(nodes, i)
    }

    ## update the slope
    slp = slope.vec[i]
  }

  ## return
  list(nodes = nodes, slopes = slopes)
}

#' Build a linear piece-wise time series function from a breakpoint dataframe
#'
#' @param t a vector of time points at which to evaluate the function
#' @param break.df a breakpoint dataframe, from which to build the function
#'
#' @return a vector of values, evaluated from the piece-wise function.
#' @export
#'
#' @examples
#' df = break_frame(nodes = c(0, 5, 9), slopes = c(.2, -.2), 1)
#' (bw.t = breakwise(0:9, df))
#' plot(0:9, bw.t, type = "l", xlab = "t", ylab = "breakwise(t)")
#'
#' # evaluate at non-integer
#' breakwise(1.32, df)
breakwise <- function(t, break.df){
  values = rep(NA, length(t)) # initialize output vector
  ## determine segment properties at each time point, and keep track of changes
  for(i in 1:length(t)){
    time = t[i] # time step
    ## which segment are we in?
    if(time == min(break.df$left)){
      seg = which(break.df$left == min(break.df$left))
    } else {
      seg =  which((break.df$left < time) & (time <= break.df$right))
    }
    ## get the value of the function at each breakpoint
    value = break.df[seg, "start"] + break.df[seg, "slope"]*(time - break.df[seg, "left"])
    values[i] = value
  }

  ## return
  values
}
