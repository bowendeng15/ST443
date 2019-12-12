setwd("./")
source("functions.R") #import own functions

### Data Simulation ########################################################
RNGkind(sample.kind = "Rounding")
set.seed(666)
p <- 50 # number of random variables
n <- 500 # number of samples generated
prob <- 0.1 # delta to construct Theta Matrix
data <- generate(p, n, prob)

table(data$E_true) 



### Apply Estimation Approaches #################################################
pred.nodewise <- predict.nodewise(data$X, lambda = 0.08)
pred.glasso <- predict.glasso(data$X, lambda = 0.08)

table(pred.nodewise$E_1, data$E_true) 
table(pred.nodewise$E_2, data$E_true) 
table(pred.glasso$E, data$E_true) 

performance(pred.nodewise$E_1, data$E_true)
performance(pred.nodewise$E_2, data$E_true)
performance(pred.glasso$E, data$E_true)



### ROC Curve and Overall Performance ###########################################
grid <- 10 ^ seq(-2, -0.8, length = 50)
perf.nodewise <- performance.nodewise.grid(data$X, data$E_true, grid)
perf.glasso <- performance.glasso.grid(data$X, data$E_true, grid)

par(mfrow=c(1,1))
plot.roc(perf.nodewise$tpr_2, perf.nodewise$fpr_2, col=1, lty=1, lwd=1)
par(new=T)
plot.roc(perf.nodewise$tpr_1, perf.nodewise$fpr_1, col=2, lty=1, lwd=1)
par(new=T)
plot.roc(perf.glasso$tpr, perf.glasso$fpr, col=4, lty=1, lwd=1)
legend("bottomright", legend=c("node1","node2","graphical"), col=c(2,1,4), lty=1)

perf.nodewise$auc_1 # 0.9941819
perf.nodewise$auc_2 # 0.9956785
perf.glasso$auc # 0.9936843



### Optimal Tuning Parameter #################################################
#### cross-validation
# compute Error Matrix
set.seed(0530)
grid <- 10 ^ seq(-0.5, -1, length = 50)
K <- 10 # number of folds
folds <- sample( rep(1:K, length = n) )
Error_1 <- Error_2 <- Error_3 <- c()
for (k in 1:K){
  cat(k, ". ")
  perf.nodewise <- performance.nodewise.grid(data$X[folds!=k,], data$E_true, grid)
  perf.glasso <- performance.glasso.grid(data$X[folds!=k,], data$E_true, grid)
  Error_1 <- rbind(Error_1, perf.nodewise$error_1)
  Error_2 <- rbind(Error_2, perf.nodewise$error_2)
  Error_3 <- rbind(Error_3, perf.glasso$error)
}

# cv plot
ylab <- "mean misclassification rate"
xlab <- expression(paste("log(",lambda,")"))
ylim <- c(0.01, 0.04)
xlim <- NULL
par(mfrow=c(1,3))
res1 <- plot.cv.error(Error_1, grid, main="Node-wise 1", xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab)
res2 <- plot.cv.error(Error_2, grid, main="Node-wise 2", xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab)
res3 <- plot.cv.error(Error_3, grid, main="Graphical", xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab)

# optimal lambdas using cv
res1$lambda.1se # 0.1490972
res2$lambda.1se # 0.1389495
res3$lambda.1se # 0.1716698

# optimise 3 approaches using optimal lambdas
pred.nodewise.1 <- predict.nodewise(data$X, lambda=res1$lambda.1se )
pred.nodewise.2 <- predict.nodewise(data$X, lambda=res2$lambda.1se )
pred.glasso.3 <- predict.glasso(data$X, lambda=res3$lambda.1se )

table(pred.nodewise.1$E_1, data$E_true) 
table(pred.nodewise.2$E_2, data$E_true) 
table(pred.glasso.3$E, data$E_true)

performance(pred.nodewise.1$E_1, data$E_true)
performance(pred.nodewise.2$E_2, data$E_true)
performance(pred.glasso.3$E, data$E_true)


# testing performance (Repetition)
set.seed(1)
fpr_1se_1 <- fnr_1se_1 <- error_1se_1 <- c()
fpr_1se_2 <- fnr_1se_2 <- error_1se_2 <- c()
fpr_1se_3 <- fnr_1se_3 <- error_1se_3 <- c()
t <- 1
while (t<=50){
  cat(t,". ")
  while (TRUE){
    test <- generate(p, n, prob)
    if ( sum(test$E_true)!=0 ) break
    }
  pred.nodewise.1se.1 <- predict.nodewise(test$X, lambda=res1$lambda.1se )
  pred.nodewise.1se.2 <- predict.nodewise(test$X, lambda=res2$lambda.1se )
  pred.glasso.1se <- predict.glasso(test$X, lambda=res3$lambda.1se )
  perf_1 <- performance(pred.nodewise.1se.1$E_1, test$E_true)
  perf_2 <- performance(pred.nodewise.1se.2$E_2, test$E_true)
  perf_3 <- performance(pred.glasso.1se$E, test$E_true)
  # append measures
  fpr_1se_1[t] <- perf_1$fpr
  fnr_1se_1[t] <- perf_1$fnr
  error_1se_1[t] <- perf_1$error
  fpr_1se_2[t] <- perf_2$fpr
  fnr_1se_2[t] <- perf_2$fnr
  error_1se_2[t] <- perf_2$error
  fpr_1se_3[t] <- perf_3$fpr
  fnr_1se_3[t] <- perf_3$fnr
  error_1se_3[t] <- perf_3$error
  t <- t+1
}
par(mfrow=c(1,3))
names <- c("node1","node2","graphical")
ylim <- NULL
boxplot(fnr_1se_1, fnr_1se_2, fnr_1se_3, main="FNR", names=rep("",3), ylim=ylim); axis(1, at=1:3, labels=names, las = 2)
boxplot(fpr_1se_1, fpr_1se_2, fpr_1se_3, main="FPR", names=rep("",3), ylim=ylim); axis(1, at=1:3, labels=names, las = 2)
boxplot(error_1se_1, error_1se_2, error_1se_3, main="Misclassification\nRate", names=rep("",3), ylim=ylim); axis(1, at=1:3, labels=names, las = 2)

sd(fpr_1se_1)
sd(fpr_1se_2)
sd(fpr_1se_3)

sd(fnr_1se_1)
sd(fnr_1se_2)
sd(fnr_1se_3)

sd(error_1se_1)
sd(error_1se_2)
sd(error_1se_3)



### Repetition ###################################################################
# it takes about 10 min to run this code chunk
set.seed(0530)
grid <- 10 ^ seq(-2, 0, length = 50)
auc_1 <- auc_2 <- auc_3 <- c()
Error_rep_1 <- Error_rep_2 <- Error_rep_3 <- c()
t <- 1
while (t<=50){
  cat(t, ". ")
  while (TRUE){ 
    data <- generate(p, n, prob) 
    if ( sum(data$E_true)!=0 ) break
    }
  perf.nodewise <- performance.nodewise.grid(data$X, data$E_true, grid)
  perf.glasso <- performance.glasso.grid(data$X, data$E_true, grid)
  auc_1[t] <- perf.nodewise$auc_1
  auc_2[t] <- perf.nodewise$auc_2
  auc_3[t] <- perf.glasso$auc
  Error_rep_1 <- rbind(Error_rep_1, perf.nodewise$error_1)
  Error_rep_2 <- rbind(Error_rep_2, perf.nodewise$error_2)
  Error_rep_3 <- rbind(Error_rep_3, perf.glasso$error)
  t <- t+1
}
min_error_1 <- apply(Error_rep_1, 1, min)
min_error_2 <- apply(Error_rep_2, 1, min)
min_error_3 <- apply(Error_rep_3, 1, min)

# boxplot
par(mfrow=c(1,2))
ylim=NULL
boxplot(auc_1, auc_2, auc_3, main="AUC", names=rep("",3), ylim=ylim); axis(1, at=1:3, labels=names, las = 2)
boxplot(min_error_1, min_error_2, min_error_3, main="Minimum \nMisclassification Rate", names=rep("",3), ylim=ylim); axis(1, at=1:3, labels=names, las = 2)



### Different Simulation Settings #################################################
# compute AUC against n
set.seed(666)
p <- 50
prob <- 0.1
grid <- 10 ^ seq(-2, 0, length = 50)
# vector of n
ns <- c( seq(10,130,length=6), seq(200,1100,length=5) ) 
auc_n_1 <- auc_n_2 <- auc_n_3 <- c()
for (n in ns){
  n <- round(n)
  cat(n, ". ")
  while (TRUE){
    data <- generate(p, n, prob) 
    if (sum(data$E_true)!=0) break
    }
  perf.nodewise <- performance.nodewise.grid(data$X, data$E_true, grid)
  perf.glasso <- performance.glasso.grid(data$X, data$E_true, grid)
  auc_n_1 <- append(auc_n_1, perf.nodewise$auc_1)
  auc_n_2 <- append(auc_n_2, perf.nodewise$auc_2)
  auc_n_3 <- append(auc_n_3, perf.glasso$auc)
}
# plot AUC against n
par(mfrow=c(1,1))
plot(ns, auc_n_1, "l", col=1, lty=1
     , xlab = "n", ylab = "AUC"
)
lines(ns, auc_n_2, col=2, lty=1)
lines(ns, auc_n_3, col=3, lty=1)
title("p=50, prob=0.1")
legend("bottomright",legend=c("node-wise 1","node-wise 2", "graphical"), col=c(1,2,3), lty=1, cex = 1)


# compute AUC against prob
set.seed(666)
n <- 500
p <- 50
grid <- 10 ^ seq(-2, 0, length = 50)
# vector of prob
probs <- 10 ^ seq(-6, -0.1, length = 20)
auc_prob_1 <- auc_prob_2 <- auc_prob_3 <- c()
for (prob in probs){
  cat(prob, ". ")
  data <- generate(p, n, prob)
  while (TRUE){ 
    data <- generate(p, n, prob) 
    if ( sum(data$E_true)!=0 ) break
    }
  perf.nodewise <- performance.nodewise.grid(data$X, data$E_true, grid)
  perf.glasso <- performance.glasso.grid(data$X, data$E_true, grid)
  auc_prob_1 <- append(auc_prob_1, perf.nodewise$auc_1)
  auc_prob_2 <- append(auc_prob_2, perf.nodewise$auc_2)
  auc_prob_3 <- append(auc_prob_3, perf.glasso$auc)
}
# plot AUC against prob
plot(log10(probs), auc_prob_1, "l", col=1, lty=1
     , ylim = NULL #c(0.7, 0.995)
     , xlab = "log(prob)", ylab = "AUC"
)
lines(log10(probs), auc_prob_2, col=2, lty=1)
lines(log10(probs), auc_prob_3, col=3, lty=1)
title("p=50, n=500")
legend("bottomleft",legend=c("node-wise 1","node-wise 2", "graphical"), col=c(1,2,3), lty=1, cex = 1)