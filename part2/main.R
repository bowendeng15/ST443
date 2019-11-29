setwd("/Users/Bowen.Deng/Desktop/LSE/ST443/group_project/ST443/part2")
source("functions.R") #import own functions

### Data Simulation ########################################################
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

neg <- sum(data$E_true==0)
pos <- sum(data$E_true==1)
(fpr_1 = sum( pred.nodewise$E_1==1 & data$E_true==0 ) / neg)
(fnr_1 = sum( pred.nodewise$E_1==0 & data$E_true==1 ) / pos)
(error_1 = ( fpr_1*neg+fnr_1*pos )/( neg+pos ))
(fpr_2 = sum( pred.nodewise$E_2==1 & data$E_true==0 ) / neg)
(fnr_2 = sum( pred.nodewise$E_2==0 & data$E_true==1 ) / pos)
(error_2 = ( fpr_2*neg+fnr_2*pos )/( neg+pos ))
(fpr_3 = sum( pred.glasso$E==1 & data$E_true==0 ) / neg)
(fnr_3 = sum( pred.glasso$E==0 & data$E_true==1 ) / pos)
(error_3 = ( fpr_3*neg+fnr_3*pos )/( neg+pos ))


### ROC Curve and Overall Performance ###########################################
grid <- 10 ^ seq(-2, 0, length = 50)
perf.nodewise <- performance.nodewise.grid(data$X, data$E_true, grid)
perf.glasso <- performance.glasso.grid(data$X, data$E_true, grid)

par(mfrow=c(1,3))
plot.roc(perf.nodewise$tpr_1, perf.nodewise$fpr_1, "ROC of \nNode-wise 1")
plot.roc(perf.nodewise$tpr_2, perf.nodewise$fpr_2, "ROC of \nNode-wise 2")
plot.roc(perf.glasso$tpr, perf.glasso$fpr, "ROC of \nGraphical Lasso")

perf.nodewise$auc_1
perf.nodewise$auc_2
perf.glasso$auc

par(mfrow=c(1,1))
plot(log10(grid), perf.glasso$error, type="l", ylim=c(0, 0.1))
abline(v=log10(grid)[which.min(perf.glasso$error)], lty=3)
lines(log10(grid), perf.glasso$fpr, col=2)
lines(log10(grid), 1-perf.glasso$tpr, col=8)

### Optimal Tuning Parameter #################################################
#### cross-validation
# compute Error Matrix
set.seed(0530)
grid <- 10 ^ seq(-0.8, -1.8, length = 50)
K <- 10 # number of folds
folds <- sample( rep(1:K, length = n) )
Error_1 <- c()
Error_2 <- c()
Error_3 <- c()
for (k in 1:K){
  cat(k, ". ")
  perf.nodewise <- performance.nodewise.grid(data$X[folds!=k,], data$E_true, grid)
  perf.glasso <- performance.glasso.grid(data$X[folds!=k,], data$E_true, grid)
  Error_1 <- rbind(Error_1, perf.nodewise$error_1)
  Error_2 <- rbind(Error_2, perf.nodewise$error_2)
  Error_3 <- rbind(Error_3, perf.glasso$error)
}

# cv plot
ylab = "mean misclassification rate"
xlab = expression(paste("log(",lambda,")"))
ylim = c(0, 0.2)
xlim = NULL
par(mfrow=c(1,3))
res1 = plot.cv.error(Error_1, grid, main="Node-wise 1", xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab)
res2 = plot.cv.error(Error_2, grid, main="Node-wise 2", xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab)
res3 = plot.cv.error(Error_3, grid, main="Graphical", xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab)

# optimal lambdas using cv
res1$lambda.1se
res2$lambda.1se
res3$lambda.1se

# performance using these lambda
pred.nodewise.1 <- predict.nodewise(data$X, lambda=round(res1$lambda.1se,4) )
pred.nodewise.2 <- predict.nodewise(data$X, lambda=round(res2$lambda.1se,4) )
pred.glasso <- predict.glasso(data$X, lambda=round(res3$lambda.1se,4) )

table(pred.nodewise.1$E_1, data$E_true) 
table(pred.nodewise.2$E_2, data$E_true) 
table(pred.glasso$E, data$E_true)

neg <- sum(data$E_true==0)
pos <- sum(data$E_true==1)
(fpr_1 = sum( pred.nodewise.1$E_1==1 & data$E_true==0 ) / neg)
(fnr_1 = sum( pred.nodewise.1$E_1==0 & data$E_true==1 ) / pos)
(error_1 = ( fpr_1*neg+fnr_1*pos )/( neg+pos ))
(fpr_2 = sum( pred.nodewise.2$E_2==1 & data$E_true==0 ) / neg)
(fnr_2 = sum( pred.nodewise.2$E_2==0 & data$E_true==1 ) / pos)
(error_2 = ( fpr_2*neg+fnr_2*pos )/( neg+pos ))
(fpr_3 = sum( pred.glasso$E==1 & data$E_true==0 ) / neg)
(fnr_3 = sum( pred.glasso$E==0 & data$E_true==1 ) / pos)
(error_3 = ( fpr_3*neg+fnr_3*pos )/( neg+pos ))


# testing performance (Repetition)
set.seed(0530)
fpr_1se_1=c(); fnr_1se_1=c(); fpr_1se_2 =c(); fnr_1se_2=c(); fpr_1se_3=c(); fnr_1se_3=c()
error_1se_1=c(); error_1se_2=c(); error_1se_3=c()
t = 1
while (t<=50){
  cat(t,". ")
  test <- generate(p, n, prob)
  while ( sum(test$E_true)==0 ){ test <- generate(p, n, prob) }
  pred.nodewise.1se.1 <- predict.nodewise(test$X, lambda=round(res1$lambda.1se,4) )
  pred.nodewise.1se.2 <- predict.nodewise(test$X, lambda=round(res2$lambda.1se,4) )
  pred.glasso.1se <- predict.glasso(test$X, lambda=round(res3$lambda.1se,4) )
  neg <- sum(test$E_true==0)
  pos <- sum(test$E_true==1)
  fpr_1se_1[t] = sum( pred.nodewise.1se.1$E_1==1 & test$E_true==0 )/neg
  fnr_1se_1[t] = sum( pred.nodewise.1se.1$E_1==0 & test$E_true==1 )/pos
  fpr_1se_2[t] = sum( pred.nodewise.1se.2$E_2==1 & test$E_true==0 )/neg
  fnr_1se_2[t] = sum( pred.nodewise.1se.2$E_2==0 & test$E_true==1 )/pos
  fpr_1se_3[t] = sum( pred.glasso.1se$E==1 & test$E_true==0 )/neg
  fnr_1se_3[t] = sum( pred.glasso.1se$E==0 & test$E_true==1 )/pos
  error_1se_1[t] = ( fpr_1se_1[t]*neg+fnr_1se_1[t]*pos )/( neg+pos )
  error_1se_2[t] = ( fpr_1se_2[t]*neg+fnr_1se_2[t]*pos )/( neg+pos )
  error_1se_3[t] = ( fpr_1se_3[t]*neg+fnr_1se_3[t]*pos )/( neg+pos )
  t <- t+1
}
par(mfrow=c(1,3))
ylim = c(0, 0.1)
names = c("node1","node2","graphical")
boxplot(fpr_1se_1, fpr_1se_2, fpr_1se_3, main="FPR", names=rep("",3), ylim=ylim); axis(1, at=1:3, labels=names, las = 2)
boxplot(fnr_1se_1, fnr_1se_2, fnr_1se_3, main="FNR", names=rep("",3), ylim=ylim); axis(1, at=1:3, labels=names, las = 2)
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
auc_1 <- c(); auc_2 <- c(); auc_3 <- c()
Error_1 <- c(); Error_2 <- c(); Error_3 <- c()
t = 1
while (t<=50){
  cat(t, ". ")
  data <- generate(p, n, prob)
  while ( sum(data$E_true)==0 ){ data <- generate(p, n, prob) }
  perf.nodewise <- performance.nodewise.grid(data$X, data$E_true, grid)
  perf.glasso <- performance.glasso.grid(data$X, data$E_true, grid)
  auc_1 <- append(auc_1, perf.nodewise$auc_1)
  auc_2 <- append(auc_2, perf.nodewise$auc_2)
  auc_3 <- append(auc_3, perf.glasso$auc)
  Error_1 <- rbind(Error_1, perf.nodewise$error_1)
  Error_2 <- rbind(Error_2, perf.nodewise$error_2)
  Error_3 <- rbind(Error_3, perf.glasso$error)
  t <- t+1
}
min_error_1 <- apply(Error_1, 1, min)
min_error_2 <- apply(Error_2, 1, min)
min_error_3 <- apply(Error_3, 1, min)

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
grid <- 10 ^ seq(0.5,-2.5, length = 20)
# vector of n
ns <- c(seq(10,130,length = 6), seq(200,1100,length = 5)) 
auc_n_1 <- c()
auc_n_2 <- c()
auc_n_3 <- c()
for (n in ns){
  n <- round(n)
  cat(n, ", ")
  data <- generate(p, n, prob)
  while ( sum(data$E_true)==0 ){ data <- generate(p, n, prob) }
  perf.nodewise <- performance.nodewise.grid(data$X, data$E_true, grid)
  perf.glasso <- performance.glasso.grid(data$X, data$E_true, grid)
  auc_n_1 <- append(auc_n_1, perf.nodewise$auc_1)
  auc_n_2 <- append(auc_n_2, perf.nodewise$auc_2)
  auc_n_3 <- append(auc_n_3, perf.glasso$auc)
}
# plot AUC against n
plot(ns, auc_n_1, "l", col=1, lty=1
     , xlab = "n", ylab = "AUC"
)
lines(ns, auc_n_2, col=2, lty=4)
lines(ns, auc_n_3, col=3, lty=4)
title("p=50, prob=0.1")
legend("bottomright",legend=c("node-wise 1","node-wise 2", "graphical"), col=c(1,2,3), pch="-", cex = 1)

# compute AUC against prob
set.seed(666)
n <- 500
p <- 50
grid <- 10 ^ seq(0.5,-2.5, length = 20)
# vector of prob
probs <- 10 ^ seq(-6,-3, length = 10)
auc_prob_1 <- c()
auc_prob_2 <- c()
auc_prob_3 <- c()
for (prob in probs){
  cat(prob, ", ")
  data <- generate(p, n, prob)
  while ( sum(data$E_true)==0 ){ data <- generate(p, n, prob) }
  perf.nodewise <- performance.nodewise.grid(data$X, data$E_true, grid)
  perf.glasso <- performance.glasso.grid(data$X, data$E_true, grid)
  auc_prob_1 <- append(auc_prob_1, perf.nodewise$auc_1)
  auc_prob_2 <- append(auc_prob_2, perf.nodewise$auc_2)
  auc_prob_3 <- append(auc_prob_3, perf.glasso$auc)
}
# plot AUC against prob
plot(log10(probs), auc_prob_1, "l", col=1, lty=1
     , ylim = c(0.7, 0.995)
     , xlab = "log10(prob)", ylab = "AUC"
)
lines(log10(probs), auc_prob_2, col=2, lty=1)
lines(log10(probs), auc_prob_3, col=3, lty=1)
title("p=50, n=500")
legend("bottomleft",legend=c("node-wise 1","node-wise 2", "graphical"), col=c(1,2,3), pch="-", cex = 1)