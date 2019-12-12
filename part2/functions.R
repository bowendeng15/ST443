library(MASS) #mvrnorm
library(matrixcalc) #is.positive.definite
library(glmnet) # glmnet 3.0-1
library(glasso)

##' @title Generate inverse covariance matrix  and mutilvariate normal data 
##' @param p number of variables
##' @param n number of samples
##' @param prob probability of each entry of inverse covariance matrix being non-zero
##' @return a list object with attributes
##' @name Theta inverse convariance matrix
##' @name Sigma covariance matrix
##' @name X mutilvariate normal data
##' @name E_true set E (vector form)
generate <- function(p, n, prob=0.1){
  # generate positive definite Theta
  delta <- 3 # initial delta
  while (TRUE) {
    Theta <- matrix(0,p,p)
    num_edge <- p*(p-1)/2 # C(p,2)
    Theta[upper.tri(Theta)] <- 0.5*rbinom(num_edge, 1, prob)
    Theta <- Theta + t(Theta) # symmetric
    Theta <- (Theta + delta*diag(p))/delta
    # ensure Theta positive definite
    if ( is.positive.definite(Theta) ) break
    delta <- delta + 1
  }
  obj <- list()
  obj$Theta <- Theta
  obj$Sigma <- solve(Theta)
  obj$E_true <- Theta[upper.tri(Theta)]!=0 # derive set E
  obj$X <- mvrnorm(n=n, mu=rep(0, p), Sigma=obj$Sigma) # simulate mvnormal data
  return(obj)
}

##' @title Compute auc using vectors of TPRs and FPR
##' @param TPR vector of TPRs
##' @param FPR vector of FPRs
##' @return value of auc
auc <- function(TPR, FPR){
  TPR <- c(0, sort(TPR), 1)
  FPR <- c(0, sort(FPR), 1)
  dFPR <- c(diff(FPR), 0)
  dTPR <- c(diff(TPR), 0)
  return( sum(TPR * dFPR) + sum(dTPR * dFPR)/2 ) 
}

##' @title Performance of a prediction
##' @param predict vector of prediction
##' @param actual vector of actual outcome
##' @return a list object with attributes
##' @name fpr 
##' @name fnr 
##' @name error
performance <- function(predict, actual){
  neg <- sum(actual==0)
  pos <- sum(actual==1)
  obj <- list()
  obj$fpr = sum( predict==1 & actual==0 ) / neg
  obj$fnr = sum( predict==0 & actual==1 ) / pos
  obj$error = ( obj$fpr*neg+obj$fnr*pos )/( neg+pos )
  return(obj)
}

##' @title Predict the E set using nodewise lasso approach
##' @param X mutilvariate normal data
##' @param lambda parameter in lasso
##' @return a list object with attributes
##' @name Beta matrix of estimated beta
##' @name E_1 estimation of E using approach 1 (vector form)
##' @name E_2 estimation of E using approach 2 (vector form)
predict.nodewise <- function(X, lambda){
  Beta <- c() # estimate beta
  p <- ncol(X)
  for (j in 1:p){
    lasso <- glmnet(X[,-j], X[,j], lambda = lambda)
    beta_j <- coef(lasso)[2:p]
    Beta <- rbind(Beta, append(beta_j, NA, after=j-1) )  # NA just a placeholder
  }
  obj <- list()
  obj$Beta <- Beta
  obj$E_1 <- Beta[upper.tri(Beta)]!=0
  obj$E_2 <- Beta[upper.tri(Beta)]!=0 & t(Beta)[upper.tri(Beta)]!=0
  return(obj)
}

##' @title Compute the Perfomance with a grid of lambdas using nodewise lasso approach
##' @param X mutilvariate normal data
##' @param E_true the true E set (vector form)
##' @param grid vector of lambdas
##' @return a list object with attributes
##' @name tpr_1 vector of TPRs using approach 1
##' @name fpr_1 vector of FPRs using approach 1
##' @name error_1 vector of error rate using approach 1
##' @name auc_1 value of AUC using approach 1
##' @name tpr_2 vector of TPRs using approach 2
##' @name fpr_2 vector of FPRs using approach 2
##' @name error_2 vector of error rate using approach 2
##' @name auc_2 value of AUC using approach 2
performance.nodewise.grid <- function(X, E_true, grid){
  tpr_1 <- fpr_1 <- error_1 <- tpr_2 <- fpr_2 <- error_2 <- c()
  for (i in 1:length(grid)){
    pred.nodewise <- predict.nodewise(X, grid[i])
    perf_1 <- performance(pred.nodewise$E_1, E_true)
    perf_2 <- performance(pred.nodewise$E_2, E_true)
    # append measures
    tpr_1[i] <- 1-perf_1$fnr
    fpr_1[i] <- perf_1$fpr
    error_1[i] <- perf_1$error
    tpr_2[i] <- 1-perf_2$fnr
    fpr_2[i] <- perf_2$fpr
    error_2[i] <- perf_2$error
  }
  obj <- list(tpr_1=tpr_1, fpr_1=fpr_1, error_1=error_1,
              tpr_2=tpr_2, fpr_2=fpr_2, error_2=error_2)
  obj$auc_1 <- auc(tpr_1, fpr_1)
  obj$auc_2 <- auc(tpr_2, fpr_2)
  return(obj)
}

##' @title Predict the E set using graphical lasso approach
##' @param X mutilvariate normal data
##' @param lambda parameter in lasso
##' @return a list object with attributes
##' @name E estimation of E (vector form)
predict.glasso <- function(X, lambda){
  glasso <- glasso(cov(X), rho=lambda)
  obj <- list()
  obj$E <- glasso$wi[upper.tri(glasso$wi)]!=0
  return(obj)
}

##' @title Compute the Perfomance with a grid of lambdas using graphical lasso approach
##' @param X mutilvariate normal data
##' @param E_true the true E set (vector form)
##' @param grid vector of lambdas
##' @return a list object with attributes
##' @name tpr vector of TPRs 
##' @name fpr vector of FPRs 
##' @name error vector of error rate
##' @name auc value of AUC
performance.glasso.grid <- function(X, E_true, grid){
  tpr <- fpr <- error <- c()
  for (i in 1:length(grid)){
    pred.glasso <- predict.glasso(X, grid[i])
    perf <- performance(pred.glasso$E, E_true)
    # append measures
    tpr[i] <- 1-perf$fnr
    fpr[i] <- perf$fpr
    error[i] <- perf$error
  }
  obj <- list(tpr=tpr, fpr=fpr, error=error)
  obj$auc <- auc(tpr, fpr)
  return(obj)
}

##' @title Plot the ROC curve
##' @param TPR vector of TPRs
##' @param FPR vector of FPRs
##' @param title text of title
plot.roc <- function(TPR, FPR, ...){
  TPR <- c(0, sort(TPR), 1)
  FPR <- c(0, sort(FPR), 1)
  plot(FPR, TPR, "l", xlab="FPR", ylab="TPR", xlim=c(0,1), ylim=c(0,1), ...)
  abline(a=0, b=1, col="grey", lty=2)
}

##' @title Plot mean error rate against lambdas and return optimal lambda
##' @param Error matrix of error rates
##' @param grid vector of lambdas
##' @return a list object with attributes
##' @name lambda.min lambda generating lowest error
##' @name lambda.1se lambda generating lowest error 1 standard error away
plot.cv.error <- function(Error, grid, ...){
  mean <- apply(Error, 2, mean)
  sd <- apply(Error, 2, sd)
  id_min <- which.min(mean)
  id_1se <- which(mean<(mean+sd)[id_min])[1] # 1 standard error
  id_2se <- which(mean<(mean+2*sd)[id_min])[1] # 2 standard error
  id_3se <- which(mean<(mean+3*sd)[id_min])[1] # 3 standard error
  # plot
  plot(log10(grid), mean, "l", ...)
  abline(h=mean[id_min]+sd[id_min],lty=2, col=2)
  abline(v=log10(grid)[id_min],lty=2, col=2)
  abline(v=log10(grid)[id_1se],lty=2, col=2)
  lines(log10(grid), (mean+sd), col="grey", lty=1)
  lines(log10(grid), (mean-sd), col="grey", lty=1)
  return(list(lambda.min=grid[id_min], lambda.1se=grid[id_1se], lambda.2se=grid[id_2se], lambda.3se=grid[id_3se]))
}