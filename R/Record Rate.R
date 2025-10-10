#################### IID #####################
#' @title Record Rate for IID Case
#' @description Computes the probability of observing a record at time \eqn{t} for an i.i.d. sequence.
#' @param t Time index (integer).
#' @return Probability of a record occurring at time \eqn{t}, given by:
#' \deqn{ P(R_t) = \frac{1}{t} }
#' #' @examples
#' rec_rate_iid(10)  # Finite-time record probability at t = 10
#' 0.1
rec_rate_iid = function(t){
  1/t
}

#################### DTRW #####################
#' Record Rate at Finite Time for DTRW Process
#'
#' Computes the probability that an observation at time \code{t} is a record in a discrete-time random walk (DTRW) process.
#' This function uses the survival probability of a symmetric random walk.
#'
#' @param t Integer time step \code{t >= 1}.
#'
#' @return Numeric, probability that time \code{t} is a record.
#' @export
#'
#' @details
#' This is the exact finite-time record rate, based on the survival probability of the walk.
#' It is equal to the survival probability. See (\code{\link{Survival}})
#' @seealso \code{\link{P_DTRW}} for the asymptotic approximation.
#'
#' @examples
#' PT_DTRW(10)  # Finite-time record probability at t = 10
#' 0.1761971
rec_rate_DTRW = function(t){
  Survival(t)
}


#' Asymptotic Record Rate for DTRW Process
#'
#' Approximates the record rate at large times \code{t} for a discrete-time symmetric random walk (DTRW).
#'
#' @param t Numeric or integer time step \code{t > 0}.
#'
#' @return Numeric, the asymptotic record rate given by \eqn{1 / \sqrt{\pi t}}.
#' @export
#'
#' @details
#' For large \code{t}, the probability that time \code{t} is a record approaches \eqn{1 / \sqrt{\pi t}}.
#'
#' @seealso \code{\link{PT_DTRW}} for the exact finite-time expression.
#'
#' @examples
#' rec_rate_DTRW_approximate(1000)  # Asymptotic record probability at large t
#' 0.01784124
#' rec_rate_DTRW_approximate(10)
#' 0.1784124
rec_rate_DTRW_approximate = function(t){
  1/sqrt(pi*t)
}
#################### LDM #####################
#' @title Record Rate for Linear Drift Model (LDM)
#' @description Computes the probability of observing a record at time \eqn{t} in an LDM setting under \eqn{Gumbel} underlying distribution.
#' @param theta The drift parameter (\eqn{\theta > 0}).
#' @param t Time index (integer).
#' @param scale Scale parameter of the \eqn{Gumbel} distribution.
#' @details Probability of a record occurring at time \eqn{t}, given by:
#' \deqn{ P(t) = \frac{1 - e^{-\theta/\text{scale}}}{1 - e^{-\theta t/\text{scale}}} }
#' @return a probability value
#' @examples
#' rec_rate_LDM(theta=0.5, t=50, scale=1)
#' 0.3934693
#' rec_rate_LDM(theta=0.5, t=10, scale=1)
#' 0.3961385
rec_rate_LDM = function(theta, t, location=0, scale){
  (1 - exp(-theta/scale)) / (1 - exp(-theta*t/scale))
}

#' @title Asymptotic Record Rate for LDM
#' @description Computes the limiting record rate as \eqn{ t \to \infty } in the LDM setting.
#' @param theta The drift parameter (\eqn{\theta > 0}).
#' @param scale Scale parameter of the \eqn{Gumbel} distribution.
#' @details Limiting record probability:
#' \deqn{ P(t -> \infty) = 1 - e^{-\theta/\text{scale}} }
#' @return probability
#' @examples rec_rate_YNM_approximate(theta=0.5, scale=1)
#' 0.3934693
rec_rate_LDM_approximate = function(theta, location=0,scale){
  1 - exp(-theta / scale)
}

#################### YNM-Nevzorov #####################
#' @title Record Rate for YNM-Nevzorov Model
#' @description Computes the probability of observing a record at time \eqn{t} in the YNM-Nevzorov setting.
#' @param gamma The memory parameter (\(\gamma > 1\)).
#' @param t Time index (integer).
#' @details Probability of a record occurring at time \eqn{t}, given by:
#' \deqn{ P(R_t) = \frac{\gamma^t}{\gamma} \cdot \frac{1 - \gamma}{1 - \gamma^t} }
#' @return a probability
#' @examples
#' rec_rate_YNM(gamma=1.1, t=10)
#' 0.1479504
#' rec_rate_YNM(gamma=1.4, t=10)
#' 0.2959456
rec_rate_YNM = function(gamma, t){
  (gamma^t / gamma) * (1 - gamma) / (1 - gamma^t)
}

#' @title Asymptotic Record Rate for YNM-Nevzorov Model
#' @description Computes the limiting record rate as \( t \to \infty \) in the YNM-Nevzorov model.
#' @param gamma The memory parameter (\(\gamma > 1\)).
#' @details Limiting record probability:
#' \deqn{ P(R_\infty) = 1 - \frac{1}{\gamma} }
#' @return a probability
#' @examples
#' rec_rate_YNM_approximate(gamma=1.1)
#' [1] 0.09090909
#' rec_rate_YNM_approximate(gamma=1.4)
#' [1] 0.2857143
rec_rate_YNM_approximate = function(gamma){
  1 - (1/gamma)
}
