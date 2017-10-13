#' Parameter Estimation Of Linear Regression Using Particle Filters

#' @description Estimation of the coefficients of a linear regression based on a particle filter algorithm. Given \code{Data1}, the first columnn is set
#' as the dependent variable, the remaining columns are the independent variables. Optional, the standard deviation of the error term can also be estimated.
#' Synthetic data is generated in case no data is provided.
#' @param Data1 matrix. Data set of dependent and independent variables. The first column corresponds to the dependent variable.
#' @param n integer. Number of particles, by default 500.
#' @param sigma_est logical. If \code{TRUE} takes the last row of \code{initDisPar} as prior estimation of the standard deviation, see more in \emph{Details}.
#' @param initDisPar matrix. Values a, b of the uniform distribution (via \code{runif}) for each parameter to be estimated, see more in \emph{Details}.
#'
#' @details We estimate the coefficients (in what follows: parameters) of a linear regression:
#'
#' \eqn{Y = \beta_0  + \beta_1*X_1 + ... + \beta_n*X_n  + \epsilon},  (\eqn{\epsilon \sim N(0,1)})
#'
#' using particle filter methods. The algorithm implementation follows the work of An, D., Choi, J. H., Kim, N. H. (2013) adapted for the case of linear models.
#'
#' The state-space equations are considered as follows:
#'
#' \eqn{(Eq. 1) X_{k} = a_0 + a_1 * X1_{k-1} + ... + a_n * Xn_{k-1}}
#'
#' \eqn{(Eq. 2) Y_{k} = X_{k} + \epsilon,  \epsilon \sim N(0,\sigma)}
#'
#' where, k = 2, ... , number of observations; \eqn{a_0, ... , a_n} are the parameters to be estimated (coefficients), and \eqn{\sigma =1 } by default.
#'
#' The priors of the parameters are considered uniformly distributed.
#' The minimum and maximum values for each prior should be placed in the first and second column of \code{initDisPar} respectively.
#' The first row in \code{initDisPar} is the prior of the state equation.
#' If \code{sigma_est = TRUE} the last row of \code{initDisPar} is the prior of the standard deviation.
#' The number of rows in \code{initDisPar} is the number of independient variables plus two when \code{sigma_est = FALSE} or plus three when \code{sigma_est = TRUE}.
#'
#' If \code{sigma_est = FALSE}, the algorithm assumes the standard deviation as one. If \code{sigma_est = TRUE}, the algorithm estimates the standard deviation of \eqn{\epsilon} together with the coefficients.
#' If \code{initDisPar} is missing, the initial priors are taken using \code{lm()} and \code{coeff()} plus-minus one as a reference.
#'
#' In case no \code{Data1} is provided, synthetic data set is generated automatically taking three normal i.i.d. variables, and the dependent variable is computed as in
#' \eqn{Y = 2 + 1.25*X_1 + 2.6*X_2 - 0.7*X_3 + \epsilon}
#'
#'
#'
#' @return A list containing the estimated parameter on each observation with \code{n} particles.
#' @author Christian Llano Robayo, Nazrul Shaikh.
#' @examples
#'
#' #### Using default Data1, no sigma estimation ####
#' Res <- PF_lm(n=1000L, sigma_est = FALSE)
#' lapply(Res,dim) # Structure of returning list.
#' sumRes <- sapply(1:5, function(i)
#'           summary(apply(Res[[i]],1,mean)))[,-1] # Summary of estimated parameters
#' colnames(sumRes) <-  c("a0", "a1", "a2", "a3")
#' sumRes
#' plot(apply(Res[['a2Resul']],1,mean),type="l")
#'
#' #### Using default Data1, with sigma estimation ####
#' Res2 <- PF_lm(n=1000L, sigma_est = TRUE)
#' lapply(Res2,dim) # Structure of returning list.
#' sumRes2 <- sapply(1:6, function(i)
#'            summary(apply(Res2[[i]],1,mean)))[,-1] # Summary of estimated parameters
#' colnames(sumRes2) <-  c("a0", "a1", "a2", "a3", "s")
#' sumRes2
#'
#' #### Using default Data1, given initDisPar ####
#' b0 <- matrix(c( -2, 11, # Prior of the state equation
#'                 1.9, 2, # Prior of a_0
#'                 1, 1.5, # Prior of a_1
#'                 2, 3,   # Prior of a_2
#'                 -1, 0)  # Prior of a_3
#'                 ,ncol = 2, byrow = TRUE )
#' Res3 <- PF_lm(n=1000L, sigma_est = FALSE, initDisPar = b0)
#' lapply(Res3,dim) # Structure of returning list.
#' sumRes3 <- sapply(1:5, function(i)
#'           summary(apply(Res3[[i]],1,mean)))[,-1] # Summary of estimated parameters
#' colnames(sumRes3) <-  c("a0", "a1", "a2", "a3")
#' sumRes3
#' plot(apply(Res3[['a2Resul']],1,mean),type="l")
#'
#' @export
#' @references An, D., Choi, J. H., Kim, N. H. (2013). Prognostics 101: A tutorial for particle filter-based prognostics algorithm using Matlab. Reliability Engineering & System Safety, 115, DOI: https://doi.org/10.1016/j.ress.2013.02.019
#' @references Ristic, B., Arulampalam, S., Gordon, N. (2004). Beyond the Kalman filter: particle filters for tracking applications. Boston, MA: Artech House. ISBN: 158053631X.
#' @references West, M., Harrison, J. (1997). Bayesian forecasting and dynamic models (2nd ed.). New York: Springer. ISBN: 0387947256.


PF_lm <- function( Data1, n = 500L, sigma_est = FALSE, initDisPar){

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

  gtemp <- function(i,g,a) {
    #Temporary function for resampling
    u <- stats::runif(1)
    loca <- (which(g>=u))
    a[loca[1]]
  }


  measuData <- Data1[,1]

  obs <- dim(Data1)[1]
  nv <- dim(Data1)[2]-1


  if (sigma_est == FALSE) {
    if (missing(initDisPar)){
      initDisPar <- rbind(c(min(Data1[,1]), max(Data1[,1])),
                          cbind(stats::coef(stats::lm(Y~., data = Data1)) -  1,
                                stats::coef(stats::lm(Y~., data = Data1)) + 1))
    }

    ParamName <- c('x', paste('a', 0:nv, sep = ''))
    p <- length(ParamName)

    param <- mapply(stats::runif, initDisPar[,1], initDisPar[,2], n = n, SIMPLIFY = FALSE)
    names(param) <- ParamName

    ParamResul <- param
    names(ParamResul) <- paste(ParamName, 'Resul', sep = '')

    k1 <- length(measuData)
    k <- 1

    likel <- numeric(0)
    cdf <- numeric(0)
    A <- numeric()

    while (k <= k1){

      k <- k+1

      paramPredi <- param
      paramPredi[[1]] <- Reduce("+", Map('*', param[2:(nv+2)], c(1,Data1[k-1,2:(nv+1)]))) # BE AWARE!! Not considered sigma
      likel <- stats::dnorm(measuData[k-1], paramPredi[[1]], 1)
      A <- cbind(A,likel)
      cdf <- cumsum(likel)
      cdf <- cdf/max(cdf)

      param <- lapply(1:p, function(j) sapply(1:100, function(i) gtemp(i, cdf, paramPredi[[j]]) ) )
      ParamResul <- mapply(rbind, ParamResul, param, SIMPLIFY = FALSE)
    }
  }

  else {

    if (missing(initDisPar)){
      initDisPar <- rbind(c(min(Data1[,1]), max(Data1[,1])),
                          cbind(stats::coef(stats::lm(Y~., data = Data1)) -3,
                                stats::coef(stats::lm(Y~., data = Data1)) +3),
                          c(0.5,1.1))
    }

    ParamName <- c('x', paste('a', 0:nv, sep = ''), 's')
    p <- length(ParamName)

    param <- mapply(stats::runif, initDisPar[,1], initDisPar[,2], n = n, SIMPLIFY = FALSE)
    names(param) <- ParamName

    ParamResul <- param
    names(ParamResul) <- paste(ParamName, 'Resul', sep = '')

    k1 <- length(measuData)
    k <- 1

    likel <- numeric(0)
    cdf <- numeric(0)
    A <- numeric()

    while (k <= k1){

      k <- k+1

      paramPredi <- param
      paramPredi[[1]] <- Reduce("+", Map('*', param[2:(nv+2)], c(1,Data1[k-1,2:(nv+1)])))
      likel <- stats::dnorm(measuData[k-1], paramPredi[[1]], paramPredi[[(nv+3)]])
      A <- cbind(A, likel)

      cdf <- cumsum(likel)
      cdf <- cdf/max(cdf)

      param <- lapply(1:p, function(j) sapply(1:100, function(i) gtemp(i, cdf, paramPredi[[j]]) ) )
      ParamResul <- mapply(rbind, ParamResul, param, SIMPLIFY = FALSE)

    }
  }

  #In our example, to plot the output user can make use of:
  #plot(apply(ParamResul[['a2Resul']],1,mean),type="l")
  #abline(h = 2.6, col='red')
  #A quick summary of the results:
  #apply(ParamResul[['a1Resul']],1,mean)
  #sapply(1:p, function(i) summary(apply(ParamResul[[i]],1,mean)))

  ParamResul
}

