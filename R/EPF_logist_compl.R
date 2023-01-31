#' Parameter Estimation Of Logistic Linear Models Using Evolutionary Particle Filters (EPF) Method

#' @description Estimation of the parameters of a logistic regression using a particle filter method that includes
#' two evolutionary algorithm-based steps. See \emph{Details}
#' @param Data, matrix. The matrix containing the independent variables of the linear model
#' @param Y numeric. The dependent variable
#' @param nPart  integer. The number of particles, by default 1000L
#' @param p_mut  numeric. The mutation probability in EPF, by default 0.01
#' @param p_cross numeric. The cross-over probability in EPF, by default 0.5
#' @param lmbd numeric. A number to be added and subtracted from the priors when \code{initDisPar} is not provided, by default 3
#' @param thr numeric. The threshold after which the estimated Y is classified as 1
#' @param count numeric. The number of replications within EPF, by default 10
#' @param resample_met function. The resampling method used in EPF, by default \code{syst_rsmpl}, see \emph{Details}
#' @param initDisPar matrix. Values a, b of the uniform distribution (via \code{runif}) for each parameter to be estimated, see more in \emph{Details}
#'
#' @details Estimation of the coefficients of a logistic linear regression using the  evolutionary particle filter.
#' EPF includes evolutionary algorithms (mutation and cross-over) within PF-base algorithm to avoid degeneracy of particles and prevent the curse of dimensionality.
#' In mutation, \code{p_mut * nPart} particles are selected at random from the set of particles obtained in \code{k-1} PF-iteration.
#' This particles will be replaced by fresh new particles taken from a uniform distribution.
#' In cross-over, a pair of random particles is selected from the \code{k-1} propagated particles. The pair is combined into one particle (using the mean) and the result replaces a particle selected at random from the \code{k-1} propagated particles.
#' This cross-over process is replicated \code{prob_cross * N} times on each iteration.
#' \code{EPF_logist_compl} uses all the information to estimate the parameters.
#'
#' The state-space equations of the logistic regression model are:
#'
#' \eqn{(Eq. 1) X_{k} = a_0 + a_1 * X1_{k-1} + ... + a_n * Xn_{k-1}}
#'
#' \eqn{(Eq. 2) Y_{k} =  1/(1 + exp(- X_{k}))}
#'
#' \eqn{(Eq. 3) Z_{k} \sim  Ber(Y_{k})}
#'
#' where, k = 1, ... , number of observations; \eqn{a_0, ... , a_n} are the parameters to be estimated (coefficients), and \eqn{Ber(\cdot)} is the Bernoulli distribution.
#'
#' The priors of the parameters are assumed uniformly distributed. In \code{initDisPar}, the number of rows is the number of independent variables plus one since the constant term of the regression is also estimated.
#' The first and second column of \code{initDisPar} are the corresponding arguments \code{a} and \code{b} of the uniform distribution (\code{stats::runif}) for each parameter prior.
#' The first row of \code{initDisPar} is the prior guess of the constant term. The following rows are the prior guesses of the coefficients.
#' If \code{initDisPar} is missing, the initial priors are taken using \code{glm()} and \code{coeff()} plus-minus \code{lmbd} as a reference.
#'
#' For \code{resample_met}, several resampling methods are used following Li, T., Bolic, M., & Djuric, P. M. (2015). Resampling methods for particle filtering: classification, implementation, and strategies. IEEE Signal processing magazine, 32(3), 70-86.
#' Currently available resampling methods are: \code{syst_rsmpl} - systematic resampling (default), \code{multin_rsmpl} - multinomial resampling,
#' \code{strat_rsmpl} - stratified resampling, and \code{simp_rsmpl} - simple resampling (similar to PF_lm_ss)
#'
#' @return A list with the following elements:
#' \code{smmry}: A summary matrix of the estimated parameters; first column is the GLM estimation, and second column is the EPF estimation.
#' \code{AIC}:  Akaike information criteria of the EPF estimation.
#' \code{Yhat}: Estimation of the observation Y.
#' \code{sttp}: A matrix with the last iteration parameters-particles.
#' \code{estm}: A vector with the final parameter EPF estimation.
#' \code{c_time}: The time in seconds it takes to complete EPF.
#'
#' @author Christian Llano Robayo, Nazrul Shaikh.
#' @examples
#'
#' \dontrun{
#' #Simulating logistic regression model
#' val_coef <- c(2, -1.25, 2.6, -0.7, -1.8)
#' Data1 <-  MASS::mvrnorm(100, mu = rep(0, 4),
#'                        Sigma = diag(4),
#'                        empirical = TRUE)
#' Y1 <- colSums(t(cbind(1,Data1))*val_coef)
#' prb <- 1 / (1 + exp(-Y1))
#' Y_ <- rbinom(100, size = 1, prob = prb) ##
#' #Run EPF
#' Res <- EPF_logist_compl(Data = Data1, Y = Y_)
#' #Summary EPF estimation
#' Res$smmry
#' }
#' @export

EPF_logist_compl <- function(Data, Y, nPart = 1000L,
                             p_mut = 0.01, p_cross = 0.5,
                             thr =.5, lmbd = 3, count = 10,
                             initDisPar, resample_met = syst_rsmpl)
  {
  meth_rsmpl <- resample_met
  p_select <- 0
  p_avrg <- 1

  #Iteration 0
  nv <- ncol(Data)
  n_obs <- nrow(Data)
  movH <- nrow(Data)
  Data1 <- Data[1:movH, ,drop = FALSE]
  names(Data1) <- paste("X", 1:nv, sep = "")
  Y1 <- Y[1:movH]
  D1 <- as.data.frame(Data1)
  pred <- matrix(NA, nrow(Data) - movH + 1, nv + 1)

  if (missing(initDisPar)) {
  initDisPar <-  cbind(stats::coef(stats::glm(Y1 ~ ., data = D1, family = "binomial")) - lmbd,
                       stats::coef(stats::glm(Y1 ~ ., data = D1, family = "binomial")) + lmbd)
  }
  st_temp <- mapply(stats::runif, initDisPar[, 1], #generating uniform distribution for each beta
                    initDisPar[,  2],
                    n = nPart, SIMPLIFY = TRUE)
  j <- 1

  tick <- Sys.time()
  while(j < count + 1) {
    Yh0 <- t(cbind(1, Data1) %*% t(st_temp))
    Yh <- logit(Yh0)
    tmplkl <- apply(Yh, 1, function(y) -log(stats::dbinom(Y[1:movH], 1,prob = y )))
    likl <- t(tmplkl)
    A <- apply(likl,1,sum)
    A_1 <- 1 / (0.0001+A)
    A_2 <- A_1 / sum(A_1)

    # Select top p_select% as is
    rnk <- rank(A_2, ties.method = "random")
    bst_idx <- which(rnk > (nPart*(1-p_select)))
    st_temp_bst_idx <- st_temp[bst_idx,]
    w_bst_idx <- A_2[bst_idx]

    temp1 <- meth_rsmpl(A_2)
    st_temp <- st_temp[temp1,]

    #Cross-over
    n_st_tmp <- 1:nrow(st_temp)
    k <-1

    while(k < nPart*p_cross){
      a1 <- sample(n_st_tmp, 2)
      b1<- colMeans(st_temp[a1,])
      a3 <- sample(n_st_tmp, 1)
      st_temp[a3,] <- b1
      k <- k + 1
    }

    # Mutation
    if(j < count & p_mut>0){
      s1 <- sample(n_st_tmp, size =  p_mut*nPart)
      temp_part <- mapply(stats::runif, initDisPar[, 1],
                          initDisPar[,  2],
                          n = p_mut*nPart, SIMPLIFY = TRUE)
      st_temp[s1,] <- temp_part
    }

    st_temp <- rbind(st_temp,st_temp_bst_idx)

    j<-j+1
  }

  #Computing weights of last iteration
  Yh0 <- t(cbind(1, Data1) %*% t(st_temp))
  Yh <- logit(Yh0)
  tmplkl <- apply(Yh, 1, function(y) -log(stats::dbinom(Y[1:movH], 1, prob = y) ))
  likl <- t(tmplkl)

  A <- apply(likl, 1, sum)
  A_1 <- 1 / (0.0001 + A)
  A_2 <- A_1 / sum(A_1)

  #Averaging top p_avrg
  bst_idx <- which(rank(A_2, ties.method = "random")>nPart*(1-p_avrg))
  st_temp_bst_idx <- st_temp[bst_idx,]
  w_bst_idx <- A_2[bst_idx]
  w_bst_idx <- w_bst_idx / sum(w_bst_idx)

  est_wh <- colMeans(st_temp_bst_idx)
  #Estimation
  pred[1,] <- apply(st_temp_bst_idx,2, function(x) sum(x*w_bst_idx))

  p1<- logit(cbind(1,Data) %*% colMeans(st_temp_bst_idx))
  LL1 <- -( sum(Y*log(p1)) + sum((1-Y)*log(1- p1)))
  McFadden_R2 <- sum(cbind(1,Data) %*% colMeans(st_temp_bst_idx))^2 / (2*LL1)
  AIC1 <- 2*LL1 +2*(nv+1)

  p1w <- logit(cbind(1,Data) %*% t(pred))
  LL1w <- -(sum(Y*log(p1w)) + sum((1-Y)*log(1- p1w)))
  McFadden_R2_w <- sum(cbind(1,Data)%*%t(pred))^2/(2*LL1w)
  AIC1w <- 2*LL1w + 2*(nv+1)

  Yhat_ <- ifelse(logit(cbind(1,Data)%*%colMeans(st_temp_bst_idx)) > thr, 1, 0)
  Yhat_w_ <- ifelse(logit(cbind(1,Data)%*%t(pred))>thr,1,0)
  tock <- Sys.time()-tick
  #Returning values
  list(smmry = data.frame(
    GLM_est = stats::coef(stats::glm(Y~Data,family = "binomial")), #OLS using all variables
    PF_est = as.numeric(est_wh)),
    AIC = AIC1,
    Yhat = Yhat_,
    sttp = st_temp,
    estm = est_wh,
    ll = LL1,
    McFadden_R2 = McFadden_R2,
    c_time = tock
  )
}
