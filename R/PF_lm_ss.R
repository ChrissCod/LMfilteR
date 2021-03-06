#' Parameter Estimation Of Linear Regression Using Particle Filters With Simple Sampling

#' @description Estimation of the coefficients of a linear regression based on the particle filters algorithm. This function is similar to \code{PF_lm} except for the resampling method which in this case is the simple sampling.
#' As a result, the user can try higher number of particles.
#' @param Data1 matrix. The matrix containing the independent variables
#' @param Y numeric. The response variable
#' @param n integer. Number of particles, by default 500
#' @param sigma_l The variance of the normal likelihood, 1 by default
#' @param sigma_est logical. If \code{TRUE} takes the last row of \code{initDisPar} as prior estimation of the standard deviation, see more in \emph{Details}
#' @param initDisPar matrix. Values a, b of the uniform distribution (via \code{runif}) for each parameter to be estimated, see more in \emph{Details}
#' @param lbd numeric. A number to be added and substracted from the priors when \code{initDisPar} is not provided
#'
#' @details Estimation of the coefficients of a linear regression:
#'
#' \eqn{Y = \beta_0  + \beta_1*X_1 + ... + \beta_n*X_n  + \epsilon},  (\eqn{\epsilon \sim N(0,1)})
#'
#' using particle filter methods.
#' The state-space equations are:
#'
#' \eqn{(Eq. 1) X_{k} = a_0 + a_1 * X1_{k-1} + ... + a_n * Xn_{k-1}}
#'
#' \eqn{(Eq. 2) Y_{k} = X_{k} + \epsilon,  \epsilon \sim N(0,\sigma)}
#'
#' where, k = 2, ... , number of observations; \eqn{a_0, ... , a_n} are the parameters to be estimated (coefficients), and \eqn{\sigma = 1 }.
#'

#' The priors of the parameters are assumed uniformly distributed. \code{initDisPar} is a matrix for which the number of rows is the number of independent variables plus one when \code{sigma_est = FALSE},
#' (plus one since we also estimate the constant term of the regression), or plus two when \code{sigma_est = TRUE} (one for the constant term of the regression, and one for the estimation of sigma).
#' The first and second column of \code{initDisPar} are the corresponding arguments \code{a} and \code{b} of the uniform distribution (\code{stats::runif}) of each parameter prior.
#' The first row \code{initDisPar} is the prior guess of the constant term. The following rows are the prior guesses of the coefficients.
#' When sigma is estimated, i.e., if \code{sigma_est = TRUE} the last row of \code{initDisPar} corresponds to the prior guess for the standard deviation.
#' If \code{sigma_est = FALSE}, then the standard deviation of the likelihood is \code{sigma_l}.
#' If \code{sigma_est = TRUE}, the algorithm estimates the standard deviation of \eqn{\epsilon} together with the coefficients.
#' If \code{initDisPar} is missing, the initial priors are taken using \code{lm()} and \code{coeff()} plus-minus \code{lbd} as a reference.
#'
#' The resampling method used corresponds to the simple sampling, i.e., we take a sample of \code{n} particles with probability equals the likelihood computed on each iteration.
#' In addition, on each iteration white noise is added to avoid particles to degenerate.
#'
#' In case no \code{Data1} is provided, a synthetic data set is generated automatically taking three normal i.i.d. variables, and the dependent variable is computed as in
#' \eqn{Y = 2 + 1.25*X_1 + 2.6*X_2 - 0.7*X_3 + \epsilon}
#'
#' @return A list containing the following elemnts:
#'
#' \code{stateP_res}: A list of matrices with the PF estimation of the parameters on each observation; the number of rows is the number of observations in \code{Data1} and the number of columns is \code{n}.
#'
#' \code{Likel}: A matrix with the likelihood of each particle obtained on each observation.
#'
#' @author Christian Llano Robayo, Nazrul Shaikh.
#' @examples
#'
#' \dontrun{
#' #### Using default Data1, no sigma estimation ####
#' Res <- PF_lm_ss(n = 10000L, sigma_est = FALSE) #10 times more than in PF_lm
#' lapply(Res,class) # Structure of returning list.
#' ###Summary of estimated parameters
#' Res$summ
#' #Evolution of the estimated parameters
#' par(mfrow=c(2, 2))
#' for (i in 1:4){
#' plot(apply(Res$stateP_res[[i]],1,mean), main = colnames(Res$summ)[i], col="blue",
#'       xlab =  "", ylab = "",type="l")
#' }
#'
#'
#' #### Using default Data1, with sigma estimation ####
#' Res2 <- PF_lm_ss(n = 1000L, sigma_est = TRUE)
#' lapply(Res2,class) # Structure of returning list
#' ###Summary of the estimated parameters
#' Res2$summ
#' #Evolution of the estimated parameters
#' par(mfrow=c(2, 3))
#' for (i in 1:5){
#' plot(apply(Res2$stateP_res[[i]],1,mean), main = colnames(Res2$summ)[i], col="blue",
#'        xlab =  "", ylab = "",type="l")
#' }
#'
#'
#' #### Using default Data1, given initDisPar ####
#' b0 <- matrix(c(1.9, 2, # Prior of a_0
#'                1, 1.5, # Prior of a_1
#'                2, 3,   # Prior of a_2
#'                -1, 0),  # Prior of a_3
#'                ncol = 2, byrow = TRUE )
#' Res3 <- PF_lm_ss(n = 10000L, sigma_est = FALSE, initDisPar = b0)
#' lapply(Res3,class) # Structure of returning list.
#' ###Summary of the estimated parameters
#' Res3$summ
#' #Evolution of the estimated parameters
#' par(mfrow=c(2, 2))
#' for (i in 1:4){
#' plot(apply(Res3$stateP_res[[i]],1,mean), main = colnames(Res3$summ)[i], col="blue",
#'     xlab =  "", ylab = "",type="l")
#'     }
#' }
#'
#' @export
#' @references Ristic, B., Arulampalam, S., Gordon, N. (2004). Beyond the Kalman filter: particle filters for tracking applications. Boston, MA: Artech House. ISBN: 158053631X.
#' @references West, M., Harrison, J. (1997). Bayesian forecasting and dynamic models (2nd ed.). New York: Springer. ISBN: 0387947256.

PF_lm_ss <- function(Y, Data1, n = 500L, sigma_l = 1, sigma_est = FALSE, initDisPar, lbd = 2){

  if (missing(Data1) & missing(Y))
  {
    Data1 <- data.frame(MASS::mvrnorm(40, mu = rep(1,3), Sigma = diag(3), empirical = TRUE))
    eps <- stats::runif(40)
    Y <- 2 + 1.25*Data1[,1] + 2.6*Data1[,2] - 0.7*Data1[,3] + eps
  }

  if (!is.data.frame(Data1)& !is.matrix(Data1))
    stop("Argument 'Data1' must be a data frame or matrix", call. = FALSE)

  if (length(n)>1)
    stop("Argument 'n' with length > 1", call. = FALSE)

  if (!is.integer(n)) {
    n <- as.integer(round(n))
  }

  if (length(sigma_est)>1)
    stop("Argument 'sigma_est' with length > 1", call. = FALSE)

  if (!is.logical(sigma_est))
    stop("Argument 'sigma_est' must be a logic", call. = FALSE)

  if (nrow(Data1) != length(Y))
    stop("Data1 and Y with different dimensions", call. = FALSE)

  Data1 <- data.frame(Data1)
  obs <- nrow(Data1)
  nv <- ncol(Data1)

  if (sigma_est == FALSE) {
    if (length(sigma_l)!=1) stop("Argument 'sigma_l' must be a number", call. = FALSE)
    if (missing(initDisPar)){
      initDisPar <- rbind(cbind(stats::coef(stats::lm(Y ~., data = Data1)) -  lbd,
                                stats::coef(stats::lm(Y ~., data = Data1)) + lbd))
    }

    stateP <- mapply(stats::runif,
                     initDisPar[,1], initDisPar[,2], n = n, SIMPLIFY = FALSE) #Prior distributions/guess about parameters

    stateP_res <- stateP
    k <- 1
    likel <- numeric(n)
    A <- matrix(0, nrow = n, ncol = obs + 1)
    cdf_r <- matrix(0, nrow = n, ncol = obs + 1)

    while (k <= obs){
      k <- k+1
      Yh <- t(as.matrix(cbind(1, Data1[k-1, ])) %*% t(sapply(stateP, cbind)))
      likel <- -log(stats::dnorm(Y[k-1], Yh, sigma_l))
      likel <- likel / sum(likel)
      likel <- 1 / (1e-10 + likel)
      likel <- likel / sum(likel)
      A[,k] <- likel
      stateP <- lapply(1:(nv+1), function(j)
        sample(stateP[[j]], size = n, replace = TRUE, prob = likel)
        )
      stateP_res <- mapply(rbind, stateP_res, stateP, SIMPLIFY = FALSE)
    }

    # Summary of estimated parameters
    sumRes <- sapply(1:(nv+1), function(i)
      summary(apply(stateP_res[[i]],1,mean)))
    colnames(sumRes) <- paste0("b", 1:(nv+1)-1)
  }

  else {

    if (missing(initDisPar)){
      initDisPar <- rbind(cbind(stats::coef(stats::lm(Y~., data = Data1)) - lbd,
                                stats::coef(stats::lm(Y~., data = Data1)) + lbd),
                          s=c(0.95,1.1))
    }

    stateP <- mapply(stats::runif, initDisPar[,1], initDisPar[,2],
                     n = n, SIMPLIFY = FALSE)
    stateP_res <- stateP
    k <- 1
    likel <- numeric(n)
    cdf <- numeric(n)
    A <- matrix(0, n, obs+1) #likelihood on each step
    cdf_r <- matrix(0, n, obs+1) #cumulative distribution on each step

    while (k <= obs){
      #print(k)
      k <- k + 1
      Yh <- t(as.matrix(cbind(1, Data1[k-1, ])) %*% t(sapply(stateP[-(nv+2)], cbind)))
      likel <- -log(stats::dnorm(Y[k-1], Yh, stateP[[nv+2]]))
      likel <- likel / sum(likel)
      likel <- 1 / (1e-15 + likel)
      likel <- likel / sum(likel)
      A[,k] <-  likel
      stateP <- lapply(1:(nv+2), function(j)
        sample(stateP[[j]], size = n,
               replace = TRUE, prob = likel)
        )
      stateP_res <- mapply(rbind, stateP_res, stateP, SIMPLIFY = FALSE)
    }

    # Summary of estimated parameters
    sumRes <- sapply(1:(nv+2), function(i) summary(apply(stateP_res[[i]],1,mean)))
    colnames(sumRes) <- c(paste0("b", 1:(nv+1)-1), "s_hat")

  }

  list(stateP_res = stateP_res,
       Likel = A,
       CDF = cdf_r,
       summ = sumRes
  )
}
