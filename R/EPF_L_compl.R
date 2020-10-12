#' Parameter Estimation Of Linear Models Using Evolutionary Particle Filters (EPF) Method

#' @description Estimation of the parameters of a linear regression model using a particle filter method that includes two evolutionary algorithm-based steps. See \emph{Details}
#' @param Data matrix. The matrix containing the independent variables of the linear model
#' @param Y numeric. The dependent variable
#' @param nPart integer. The number of particles, by default 1000L
#' @param p_mut  numeric. The mutation probability in EPF, by default 0.01
#' @param p_cross numeric. The cross-over probability in EPF, by default 0.5
#' @param lmbd numeric. A number to be added and subtracted from the priors when \code{initDisPar} is not provided, by default 3
#' @param count numeric. The number of replications within EPF, by default 10
#' @param resample_met function. The resampling method used in EPF, by default \code{syst_rsmpl}, see \emph{Details}
#' @param sigma_l The value of the variance in the likelihood normal, 1 by default
#' @param initDisPar matrix. Values a, b of the uniform distribution (via \code{runif}) for each parameter to be estimated, see more in \emph{Details}
#'
#' @details Estimation of the coefficients of a linear regression:
#'
#' \eqn{Y = \beta_0  + \beta_1*X_1 + ... + \beta_n*X_n  + \epsilon},  (\eqn{\epsilon \sim N(0,1)})
#'
#' using Evolutionary Particle Filter. EPF includes additional steps within PF-base algorithm to avoid degeneracy of particles and prevent the curse of dimensionality.
#' One step is the inclusion of two evolutionary algorithm - based processes: mutation and cross-over. In mutation, \code{p_mut * nPart} particles are selected at random from the set of particles obtained in \code{k-1} PF-iteration.
#' This particles will be replaced by fresh new particles taken from a uniform distribution.
#' In cross-over, a pair of random particles is selected from the \code{k-1} propagated particles. The pair is combined into one particle (using the mean function) and the result replaces a particle selected at random from the \code{k-1} propagated particles.
#' This cross-over process is replicated \code{prob_cross * N} times on each iteration.
#' \code{EPF_L_compl} uses all the information to estimate the parameters, differs from function \code{EPF_L_MovH} by using \code{movH} equals \code{nrows(Data)}
#'
#' Similarly as in PF_lm and PF_lm_ss, the state-space equations of the linear model are:
#'
#' \eqn{(Eq. 1) X_{k} = a_0 + a_1 * X1_{k-1} + ... + a_n * Xn_{k-1}}
#'
#' \eqn{(Eq. 2) Y_{k} = X_{k} + \epsilon,  \epsilon \sim N(0,\sigma)}
#'
#' where, k = 2, ... , number of observations; \eqn{a_0, ... , a_n} are the parameters to be estimated (coefficients), and \eqn{\sigma = 1 }.
#'
#' The priors of the parameters are assumed uniformly distributed. \code{initDisPar} is a matrix; the number of rows is the number of independent variables plus one since the constant term of the regression is also estimated.
#' The first and second column of \code{initDisPar} are the corresponding arguments \code{a} and \code{b} of the uniform distribution (\code{stats::runif}) for each parameter prior.
#' The first row of \code{initDisPar} is the prior guess of the constant term. The following rows are the prior guesses of the coefficients.
#' If \code{initDisPar} is missing, the initial priors are taken using \code{lm()} and \code{coeff()} plus-minus \code{lmbd} as a reference.
#'
#' For \code{resample_met}, several resampling methods are used following Li, T., Bolic, M., & Djuric, P. M. (2015). Resampling methods for particle filtering: classification, implementation, and strategies. IEEE Signal processing magazine, 32(3), 70-86.
#' Currently available resampling methods are: \code{syst_rsmpl} - systematic resampling (default), \code{multin_rsmpl} - multinomial resampling,
#' \code{strat_rsmpl} - stratified resampling, and \code{simp_rsmpl} - simple resampling (similar to PF_lm_ss)
#'
#' @return A list with the following elements:
#' \code{smmry}: A summary matrix of the estimated parameters; first column is the OLS estimation, and second column is the EPF estimation.
#' \code{SSE}: The sum of squared errors of the EPF estimation.
#' \code{Yhat}: Estimation of the observation Y.
#' \code{sttp}: A matrix with the last iteration particles.
#' \code{estm}: A vector with the final parameter EPF estimation.
#' \code{c_time}: The time in seconds it takes to complete EPF.
#'
#' @author Christian Llano Robayo, Nazrul Shaikh.
#' @examples
#'
#' \dontrun{
#' #Example using linear model with 5 variables
#' #Using rcorrmatrix from cluterGeneration to generate a random correlation matrix
#' realval_coef <- c(2, -1.25, 2.6, -0.7, -1.8, 0.6)
#' Data1 <-  MASS::mvrnorm(100, mu = rep(0, 5),
#'                         Sigma = clusterGeneration::rcorrmatrix(5),
#'                         empirical = TRUE)
#' eps <- stats::rnorm(100) # error term
#' Y_ <- colSums(t(cbind(1,Data1)) * realval_coef) + eps
#' #run EPF
#' Res1 <- EPF_L_compl(Data = Data1, Y = Y_)
#' #Summary EPF estimation
#' Res1$smmry
#' }
#' @export


EPF_L_compl <- function(Data, Y, nPart = 1000L, p_mut = 0.01, p_cross = 0.5, sigma_l = 1, lmbd = 3, count = 10, initDisPar, resample_met = syst_rsmpl)
  {

  meth_rsmpl <- resample_met
  p_select <- 0
  p_avrg <- 1
  #Iteration 0
  nv <- ncol(Data)
  n_obs <- nrow(Data)
  movH <- n_obs
  Data1 <- Data[1:movH, ,drop = FALSE]
  names(Data1) <- paste("X", 1:nv, sep = "")
  Y1 <- Y[1:movH]
  D1 <- as.data.frame(Data1)
  pred <- matrix(NA,nrow(Data)-movH+1,nv+1)

  if (missing(initDisPar)) {
    initDisPar <-  cbind(stats::coef(stats::lm(Y1 ~ ., data = D1)) - lmbd,
                         stats::coef(stats::lm(Y1 ~ ., data = D1)) + lmbd)
  }
  st_temp <- mapply(stats::runif, initDisPar[, 1], #generating uniform distribution for each beta
                    initDisPar[,  2],
                    n = nPart, SIMPLIFY = TRUE)

  lsst<-list()
  j <- 1

  tick <- Sys.time()
  while(j < count+1){
    Yh <- t(cbind(1, Data1) %*% t(st_temp))
    tmplkl <- apply(Yh, 1, function(y) -log(stats::dnorm(x = y, mean = Y[1:movH], sd = sigma_l)))
    likl <- t(tmplkl)
    A <- apply(likl,1,sum)
    A_1 <- 1 / (0.00001+A)
    A_2 <- A_1/sum(A_1)

    # Select top p_select% as is
    rnk <- rank(A_2, ties.method = "random")
    bst_idx <- which( rnk > (nPart*(1-p_select)))
    st_temp_bst_idx <- st_temp[bst_idx,]
    w_bst_idx <- A_2[bst_idx]

    temp1 <- meth_rsmpl(A_2)
    st_temp <- st_temp[temp1,]

    #Cross-over
    n_st_tmp <- 1:nrow(st_temp)
    k <-1
    while(k < nPart*p_cross){
      a1 <- sample(n_st_tmp, 2)
      b1 <- colMeans(st_temp[a1,])
      a3 <- sample(n_st_tmp, 1)
      st_temp[a3,] <- b1
      k <- k+1
    }

    # Mutation
    if(j < count & p_mut>0){
      s1 <- sample(n_st_tmp, size =  p_mut*nPart)

      dif1 <- t(apply(st_temp, 2, function(x) {
        d1 <- mean(x)
        sd1 <- stats::sd(x)
        c(d1 - 2*sd1,d1 + 2*sd1)
      }))
      temp_part <- mapply(stats::runif, dif1[,1], #generating uniform distribution for each beta
                          dif1[,2], n = p_mut*nPart, SIMPLIFY = TRUE)
      st_temp[s1,] <- temp_part
    }

    st_temp <- rbind(st_temp,st_temp_bst_idx)

    j<-j+1
  }

  #Computing weights of last iteration
  Yh <- t(cbind(1, Data1) %*% t(st_temp))
  tmplkl <- apply(Yh, 1, function(y) -log(stats::dnorm(x = y, mean = Y[1:movH], sd = 1)))
  likl <- t(tmplkl)
  if(movH==1) likl <- t(likl)
  A <- apply(likl,1,sum)
  A_1 <- 1/(0.0001+A)
  A_2 <- A_1/sum(A_1)

  #Averaging top p_avrg
  bst_idx <- which(rank(A_2, ties.method = "random")>nPart*(1-p_avrg))
  st_temp_bst_idx <- st_temp[bst_idx,]
  w_bst_idx <- A_2[bst_idx]
  w_bst_idx <- w_bst_idx/sum(w_bst_idx)

  est_wh <- colMeans(st_temp_bst_idx)
  #Estimation
  pred[1,] <- apply(st_temp_bst_idx,2, function(x) sum(x*w_bst_idx))

  Yhat_ <- cbind(1,Data)%*%colMeans(st_temp_bst_idx)
  Yhat_w_ <- cbind(1,Data)%*%t(pred)
  tock <- Sys.time()-tick
  lmd1 <- stats::lm(Y~.,data = data.frame(Data))
  #6. Returning values
  list(smmry = data.frame(
    OLS_est = stats::coef(lmd1),
    PF_est = as.numeric(est_wh)),
    SSE = sum((Y - Yhat_)^2),
    Yhat = Yhat_,
    sttp = st_temp,
    estm = est_wh,
    c_time = as.numeric(tock)
  )
}


