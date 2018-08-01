#' Parameter Estimation Of Linear Regression Using Particle Filters

#' @description Estimation of the coefficients of a linear regression based on the particle filters algorithm. This function is similar to \code{PF_lm} except for the resampling method which in this case is the simple sampling.
#' As result, the user can try higher number of particles.
#' @param Data1 matrix. A matrix with the dependent and independent variables. The first column should contain the dependent variable.
#' @param n integer. Number of particles, by default 500.
#' @param sigma_l The variance of the normal likelihood, 1 by default.
#' @param sigma_est logical. If \code{TRUE} takes the last row of \code{initDisPar} as prior estimation of the standard deviation, see more in \emph{Details}.
#' @param initDisPar matrix. Values a, b of the uniform distribution (via \code{runif}) for each parameter to be estimated, see more in \emph{Details}.
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
#' The priors of the parameters are assumed uniformly distributed. \code{initDisPar} is a matrix for which the number of rows is the number of independient variables plus two when \code{sigma_est = FALSE}, (one for the state and one for the constant of the regression), or plus three when \code{sigma_est = TRUE} (one for the state,  one for the constant of the regression, and one for the estimation of sigma).
#' The first and second column of \code{initDisPar} are the corresponding arguments \code{a} and \code{b} of the uniform distribution (\code{stats::runif}) of each parameter prior.
#' The first row of \code{initDisPar} is the prior of the state equation.
#' The second row is the prior guess of the regression interecept term. The following rows are the prior guesses of the coefficients, or when is estimated, of sigma, i.e.,
#' if \code{sigma_est = TRUE} the last row of \code{initDisPar} is the uniform prior of the standard deviation.
#' If \code{sigma_est = FALSE}, then standard deviation of the likelihood is \code{sigma_l}.
#' If \code{sigma_est = TRUE}, the algorithm estimates the standard deviation of \eqn{\epsilon} together with the coefficients.
#' If \code{initDisPar} is missing, the initial priors are taken using \code{lm()} and \code{coeff()} plus-minus 1 as a reference.
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
#' Res <- PF_lm_ss(n=10000L, sigma_est = FALSE) #10 times more than in PF_lm
#' lapply(Res,class) # Structure of returning list.

#' sumRes <- sapply(2:5, function(i)
#' summary(apply(Res$stateP_res[[i]],1,mean))) # Summary of estimated parameters
#' colnames(sumRes) <-  c("a0", "a1", "a2", "a3")

#' sumRes

#' par(mfrow=c(2, 2)) #Evolution of the mean of PF estimation on each time
#' for (i in 2:5){
#' plot(apply(Res$stateP_res[[i]],1,mean), main = colnames(sumRes)[i-1], col="blue",
#'       xlab =  "", ylab = "",type="l")
#' }
#'
#'dev.off()

#' #### Using default Data1, with sigma estimation ####
#' Res2 <- PF_lm_ss(n = 1000L, sigma_est = TRUE)
#' lapply(Res2,class) # Structure of returning list
#' sumRes2 <- sapply(2:6, function(i)
#'   summary(apply(Res2$stateP_res[[i]],1,mean)))# Summary of estimated parameters
#' colnames(sumRes2) <-  c("a0", "a1", "a2", "a3", "s")
#' sumRes2
#'
#' par(mfrow=c(2, 3)) #Evolution of the mean of PF estimation
#' for (i in 2:6){
#' plot(apply(Res2$stateP_res[[i]],1,mean), main = colnames(sumRes2)[i-1], col="blue",
#'        xlab =  "", ylab = "",type="l")
#' }
#' dev.off()
#'
#' #### Using default Data1, given initDisPar ####
#' b0 <- matrix(c( -2, 11, # Prior of the state equation
#'                 1.9, 2, # Prior of a_0
#'                 1, 1.5, # Prior of a_1
#'                 2, 3,   # Prior of a_2
#'                 -1, 0)  # Prior of a_3
#'                 ,ncol = 2, byrow = TRUE )
#' Res3 <- PF_lm_ss(n=10000L, sigma_est = FALSE, initDisPar = b0)
#' lapply(Res3,class) # Structure of returning list.
#' sumRes3 <- sapply(2:5, function(i)
#'           summary(apply(Res3$stateP_res[[i]],1,mean))) # Summary of estimated parameters
#' colnames(sumRes3) <-  c("a0", "a1", "a2", "a3")
#' sumRes3
#'
#' par(mfrow=c(2, 2)) #Evolution of the mean of PF estimation
#' for (i in 2:5){
#' plot(apply(Res3$stateP_res[[i]],1,mean), main = colnames(sumRes3)[i-1], col="blue",
#'     xlab =  "", ylab = "",type="l")
#'     }
#' dev.off()
#' }
#'
#' @export
#' @references Ristic, B., Arulampalam, S., Gordon, N. (2004). Beyond the Kalman filter: particle filters for tracking applications. Boston, MA: Artech House. ISBN: 158053631X.
#' @references West, M., Harrison, J. (1997). Bayesian forecasting and dynamic models (2nd ed.). New York: Springer. ISBN: 0387947256.


PF_lm_ss <- function( Data1, n = 500L, sigma_l = 1, sigma_est = FALSE, initDisPar){

  if (missing(Data1))
  {
    Data1 <- MASS::mvrnorm(40, mu = rep(1,3), Sigma = diag(3), empirical = TRUE)
    eps <- stats::runif(40)
    Y <- 2 + 1.25*Data1[,1] + 2.6*Data1[,2] - 0.7*Data1[,3] + eps
    Data1 <- as.data.frame(cbind(Y,Data1))
  }

  if (!is.data.frame(Data1)) stop("Argument 'Data1' must be a data frame or matrix", call. = FALSE)
  if (length(n)>1)  stop("Argument 'n' with length > 1", call. = FALSE)
  if (!is.integer(n)) stop("Argument 'n' must be an integer number", call. = FALSE)
  if (length(sigma_est)>1)  stop("Argument 'sigma_est' with length > 1", call. = FALSE)
  if (!is.logical(sigma_est)) stop("Argument 'sigma_est' must be a logic", call. = FALSE)

  YY <- Data1[,1]
  obs <- dim(Data1)[1]
  nv <- dim(Data1)[2]-1

  if (sigma_est == FALSE) {
    if (length(sigma_l)!=1) stop("Argument 'sigma_l' must be a number", call. = FALSE)
    if (missing(initDisPar)){
      names(Data1) <- c("Y",paste("X",1:nv,sep = ""))
      initDisPar <- rbind(c(min(Data1[,1]), max(Data1[,1])),
                          cbind(stats::coef(stats::lm(Y ~., data = Data1)) -  1,
                                stats::coef(stats::lm(Y ~., data = Data1)) + 1))
    }

    p <- nv + 2

    stateP <- mapply(stats::runif,
                     initDisPar[,1], initDisPar[,2], n = n, SIMPLIFY = FALSE) #Prior guess

    stateP_res <- stateP

    k1 <- length(YY)
    k <- 1

    likel <- numeric(n)
    A <- matrix(0, n, obs+1)

    while (k <= k1){

      k <- k+1

      statePrdt <- stateP
      statePrdt[[1]] <- Reduce("+", Map('*', stateP[2:(nv+2)], c(1, Data1[k-1,2:(nv+1)]))) #Not considered sigma

      likel <- stats::dnorm(YY[k-1],statePrdt[[1]], 1)
      A[,k] <- likel

      stateP <- lapply(1:p,function(j) sample(statePrdt[[j]],
                        size = n, replace = TRUE, prob = likel/sum(likel) ))
      stateP_res <- mapply(rbind, stateP_res, stateP, SIMPLIFY = FALSE)
      stateP[2:(nv+2)] <-  mapply("+", stateP[2:(nv+2)], lapply(1:(nv+1),
                                   function(i) stats::rnorm(n)),SIMPLIFY = FALSE)

    }
  }

  else {

    if (missing(initDisPar)){
      initDisPar <- rbind(c(min(Data1[,1]), max(Data1[,1])),
                          cbind(stats::coef(stats::lm(Y~., data = Data1)) - 2,
                                stats::coef(stats::lm(Y~., data = Data1)) + 2),
                          s=c(0.99,1.5))
    }


    stateP <- mapply(stats::runif,
                     initDisPar[,1], initDisPar[,2], n = n, SIMPLIFY = FALSE)

    stateP_res <- stateP

    k1 <- length(YY)
    k <- 1

    likel <- numeric(n)
    cdf <- numeric(n)
    A <- matrix(0, n, obs+1) #likelihood on each step
    tt <- matrix(0, n, obs+1) #state prediction

    while (k <= k1){

      k <- k+1

      statePrdt <- stateP
      statePrdt[[1]] <- Reduce("+", Map('*', stateP[2:(nv+2)], c(1,Data1[k-1,2:(nv+1)])))
      tt[,k] <- statePrdt[[1]]

      likel <- stats::dnorm(YY[k-1], statePrdt[[1]], statePrdt[[nv+3]])
      A[,k] <-  likel

      stateP <- lapply(1:(nv+3), function(j) sample(statePrdt[[j]],
                        size = n, replace = TRUE, prob = likel/sum(likel) ))
      stateP_res <- mapply(rbind, stateP_res, stateP, SIMPLIFY = FALSE)
      stateP[2:(nv+2)] <-  mapply("+", stateP[2:(nv+2)], lapply(1:(nv+1),
                                    function(i) stats::rnorm(n)),SIMPLIFY = FALSE)

    }
  }

  list(stateP_res = stateP_res, Likel = A)

}


