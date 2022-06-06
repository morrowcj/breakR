
#' Generate seasonality with a sinusoidal function
#'
#' @param t a vector of time points at which to evaluate the sinusoidal function
#' @param amp the amplitude of the sine wave
#' @param per periodicity of the wave: the number of time points over which one full
#' cycle is completed
#' @param noise.sd the standard deviation of the random normal variable, which
#' is added to the sinusoidal wave.
#' @param trend an optional linear trend coefficient to add to the sine wave.
#'
#' @details a sine wave is generated as \code{amp + sin(2*pi / per * t) + trend*t + e} where
#' \code{e} is a random normal variable with sd = \code{noise.sd}.
#'
#' @return a vector of values the same length as \code{t}
#'
#' @export
#'
#' @examples
#' season()
#' season(amp = .8)
#' season(t = seq(0, 1, .1), per = 1)
#' season(t = 15)
#' season(noise.sd = 1)
#'
#' # visualize
#' plot(season(), type = "l", ylim = c(-1, 1))
#' plot(season(amp = .5), col = "red", type = "l", ylim = c(-1, 1))
#' plot(season(amp = .5, noise.sd = .5), col = "blue", type = "l")
#' plot(season(0:100, trend = .1/12), type = "l")
season <- function(t = 0:24, amp = 1, per = 12, noise.sd = 0, trend = 0){
  ## sinusoidal wave
  sinus = amp * sin(2*pi / per * t)
  ## optional noise
  if(noise.sd > 0){
    out = sinus + rnorm(length(t), sd = noise.sd)
  } else {
    out = sinus
  }
  ## optional trend
  if(trend > 0){
    out = out + trend*t
  }
  ## return
  out
}
