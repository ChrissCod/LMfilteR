#Resampling using CDF method
gtemp <- function(g,a) {
   #Temporary function for resampling
   u <- stats::runif(1)
   loca <- which(g>=u)
   a[loca[1]]
}

##Systematic resampling
syst_rsmpl <- function(w) {
 w <- w / sum(w)
 cs <- cumsum(w)
 M <- length(w)
 y <- rep(NA, M)
 y[1] <- stats::runif(1, 0, 1 / M)
 y[1:M] <- ((1:M-1) + y[1]) / M
 idx <- numeric(M)
 k <- 1
 for (j in 1:M) {
   while(cs[k] < y[j]) {
      k <- k + 1
   }
   idx[j] <- k
   }
   idx
}

##Multinomial resampling
multin_rsmpl <- function(w) {
 w <- w / sum(w)
 M <- length(w)
 cs <- cumsum(w)
 idx <- numeric(M)
 i <- 1
 for (i in 1:M) {
   tempR <- stats::runif(1)
   j <- 1
   while (cs[j] < tempR) {
      j <- j + 1
      }
   idx[i] <- j
 }
 return(idx)
}

##Simple resampling
simp_rsmpl <- function(w) {
 w <- w / sum(w)
 idx <- sample(1:length(w), replace = TRUE, prob = w)
 idx
}

##Stratified resampling
strat_rsmpl <- function(w) {
 w <- w / sum(w)
 M <- length(w)
 cs <- cumsum(w)
 idx <- numeric(M)
 temp <- seq(from = 0, to = 1 - 1 / M, by = 1 / M) + stats::runif(M) / M
 i <- 1
 j <- 1
 while (i <= M & j <= M) {
   while (cs[j] < temp[i]) {
    j = j + 1
   }
    idx[i] <- j
    i <- i + 1
   }
   idx
}

#Logit function
logit <- function(x) {
   1 / (1 + exp(-x))
}
