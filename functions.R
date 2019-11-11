generate <- function(p, n, delta){
  obj <- list()
  obj$p <- p
  obj$n <- n
  obj$delta <- delta
  Theta <- matrix(0,p,p) + delta*diag(p)
  for (i in 2:p){
    for (j in 1:(i-1)){
      Theta[i,j] = 0.5*rbinom(1, 1, 0.1)
      Theta[j,i] = Theta[i,j]
    }
  }
  obj$Theta <- Theta/delta
  obj$Sigma <- solve(Theta)
  obj$X <- mvrnorm(n=n, mu=rep(0, p), Sigma=obj$Sigma)
  E_true <- matrix(nrow=0, ncol=3)
  for (i in 1:(p-1)){
    for (j in (i+1):p){
      if ( Theta[i,j] ){
        E_true <- rbind(E_true, c(i,j,1))
      }
      else{
        E_true <- rbind(E_true, c(i,j,0))
      }
    }
  }
  obj$E_true <- E_true
  return(obj)
}

auc <- function(TPR, FPR){
  TPR <- c(0, TPR, 1)
  FPR <- c(0, FPR, 1)
  dFPR <- c(diff(FPR), 0)
  dTPR <- c(diff(TPR), 0)
  return( sum(TPR * dFPR) + sum(dTPR * dFPR)/2 ) 
}

predict.nodewise <- function(X, lambda){
  p <- ncol(X)
  # estimate beta
  Beta <- matrix(nrow=0, ncol=p)
  for (j in 1:p){
    lasso <- glmnet(X[,-j], X[,j], lambda = lambda)
    beta_j <- coef(lasso)[2:p]
    Beta <- rbind(Beta
                  , append(beta_j, NA, after=j-1) # NA just a placeholder, meaningless
    ) 
  }
  # hence get E_1 and E_2
  E_1 <- matrix(nrow=0, ncol=3)
  E_2 <- matrix(nrow=0, ncol=3)
  for (i in 1:(p-1)){
    for (j in (i+1):p){
      if (Beta[i,j]){
        E_1 <- rbind(E_1, c(i,j,1))
      }
      else{
        E_1 <- rbind(E_1, c(i,j,0))
      }
      if (Beta[i,j] && Beta[j,i]){
        E_2 <- rbind(E_2, c(i,j,1))
      }
      else{
        E_2 <- rbind(E_2, c(i,j,0))
      }
    }
  }
  obj <- list()
  obj$Beta <- Beta
  obj$E_1 <- E_1
  obj$E_2 <- E_2
  return(obj)
}

performance.nodewise.grid <- function(X, E_true, grid){
  tpr_1 <- c()
  fpr_1 <- c()
  tpr_2 <- c()
  fpr_2 <- c()
  for (lambda in grid){
    pred.nodewise <- predict.nodewise(X, lambda)
    # compute TPR FPR
    Tab_1 <- table(pred.nodewise$E_1[,3], E_true[,3])
    if (!(0 %in% rownames(Tab_1))){
      Tab_1 <- rbind(c(0,0), Tab_1) # in case all predictions are 1
    }
    if (!(1 %in% rownames(Tab_1))){
      Tab_1 <- rbind(Tab_1, c(0,0)) # in case all predictions are 0
    }
    tpr_1 <- append(tpr_1, Tab_1[2,2]/(sum(Tab_1[,2])) )
    fpr_1 <- append(fpr_1, Tab_1[2,1]/(sum(Tab_1[,1])) ) 
    
    Tab_2 <- table(pred.nodewise$E_2[,3], E_true[,3])
    if (!(0 %in% rownames(Tab_2))){
      Tab_2 <- rbind(c(0,0), Tab_2)
    }
    if (!(1 %in% rownames(Tab_2))){
      Tab_2 <- rbind(Tab_2, c(0,0))
    }
    tpr_2 <- append(tpr_2, Tab_2[2,2]/(sum(Tab_2[,2])) )
    fpr_2 <- append(fpr_2, Tab_2[2,1]/(sum(Tab_2[,1])) )
  }
  neg <- sum(E_true[,3]==0)
  pos <- sum(E_true[,3]==1)
  fnr_1 <- 1-tpr_1
  fnr_2 <- 1-tpr_2
  obj <- list()
  obj$tpr_1 <- tpr_1
  obj$fpr_1 <- fpr_1
  obj$tpr_2 <- tpr_2
  obj$fpr_2 <- fpr_2
  obj$error_1 <- (fpr_1*pos+fnr_1*neg)/(pos+neg)
  obj$error_2 <- (fpr_2*pos+fnr_2*neg)/(pos+neg)
  obj$auc_1 <- auc(tpr_1, fpr_1)
  obj$auc_2 <- auc(tpr_2, fpr_2)
  return(obj)
}

predict.glasso <- function(X, lambda){
  glasso <- glasso(cov(X), rho=lambda)
  E <- matrix(nrow=0, ncol=3)
  for (i in 1:(p-1)){
    for (j in (i+1):p){
      if ( glasso$wi[i,j] ){
        E <- rbind(E, c(i,j,1))
      }
      else{
        E <- rbind(E, c(i,j,0))
      }
    }
  }
  obj <- list()
  obj$E <- E
  return(obj)
}

performance.glasso.grid <- function(X, E_true, grid){
  tpr <- c()
  fpr <- c()
  for (lambda in grid){
    pred.glasso <- predict.glasso(X, lambda)
    # compute TPR FPR
    Tab <- table(pred.glasso$E[,3], E_true[,3])
    if (!(0 %in% rownames(Tab))){
      Tab <- rbind(c(0,0), Tab) # in case all predictions are 1
    }
    if (!(1 %in% rownames(Tab))){
      Tab <- rbind(Tab, c(0,0)) # in case all predictions are 0
    }
    tpr <- append(tpr, Tab[2,2]/(sum(Tab[,2])) )
    fpr <- append(fpr, Tab[2,1]/(sum(Tab[,1])) ) 
  }
  neg <- sum(E_true[,3]==0)
  pos <- sum(E_true[,3]==1)
  fnr <- 1-tpr
  obj <- list()
  obj$tpr <- tpr
  obj$fpr <- fpr
  obj$error <- (fpr*pos+fnr*neg)/(pos+neg)
  obj$auc <- auc(tpr, fpr)
  return(obj)
}

plot.roc <- function(TPR, FPR, title=NULL){
  TPR <- c(0, TPR, 1)
  FPR <- c(0, FPR, 1)
  plot(FPR, TPR, "l"
       , xlab="FPR", ylab="TPR"
       , xlim=c(0,1), ylim=c(0,1)
       # , xaxs="i", yaxs="i"
  )
  abline(a=0, b=1, col="red")
  title(main = title)
}

plot.cv.error <- function(Error, grid, zoom=NULL, ...){
  mean <- apply(Error, 2, mean)
  sd <- apply(Error, 2, sd)
  id_min <- which.min(mean)
  # set zoom in area
  if ( is.null(zoom) ) id_zoom <- 1:length(grid)
  else id_zoom <- (id_min-zoom):(id_min+zoom)
  # plot
  plot(log10(grid)[id_zoom], mean[id_zoom]
       , "l", xlab="log10(lambda)", ylab="mean error rate"
       ,...)
  abline(h=mean[id_min]+sd[id_min],lty=2, ...) 
  abline(v=log10(grid)[id_min],lty=2, ...)
  lines(log10(grid)[id_zoom], (mean+sd)[id_zoom], col="grey", lty=3)
  lines(log10(grid)[id_zoom], (mean-sd)[id_zoom], col="grey", lty=3)
}