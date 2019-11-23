setwd("/Users/Bowen.Deng/Desktop/LSE/ST443/group_project/ST443/part2")
source("functions.R") #import own functions

### Data Simulation ########################################################
set.seed(666)
p <- 50 # number of random variables
n <- 1000 # number of samples generated
prob <- 0.1 # delta to construct Theta Matrix
data <- generate(p, n, prob)



### Apply Estimation Approaches #################################################
pred.nodewise <- predict.nodewise(data$X, lambda = 0.05)
table(pred.nodewise$E_1, data$E_true)
table(pred.nodewise$E_2, data$E_true) 

pred.glasso <- predict.glasso(data$X, lambda = 0.05)
table(pred.glasso$E, data$E_true) 



### ROC Curve and Overall Performance ###########################################
grid <- 10 ^ seq(-2.5, 0.5, length = 50)
perf.nodewise <- performance.nodewise.grid(data$X, data$E_true, grid)
perf.glasso <- performance.glasso.grid(data$X, data$E_true, grid)

par(mfrow=c(1,3))
plot.roc(perf.nodewise$tpr_1, perf.nodewise$fpr_1, "ROC of Node-wise 1")
plot.roc(perf.nodewise$tpr_2, perf.nodewise$fpr_2, "ROC of Node-wise 2")
plot.roc(perf.glasso$tpr, perf.glasso$fpr, "ROC of Graphical Lasso")

perf.nodewise$auc_1
perf.nodewise$auc_2
perf.glasso$auc



### Optimal Tuning Parameter #################################################
#### no repetition 
# plot
plot(log10(grid), perf.nodewise$error_1, type="l", col=2
     , xlab = "log10(lambda)", ylab = "error rate"
)
lines(log10(grid), perf.nodewise$error_2, col=3)
lines(log10(grid), perf.glasso$error, col=4)
legend("bottomright",legend=c("node-wise 1","node-wise 2", "graphical"), col=c(2,3,4), pch="-")
lambda_1 <- grid[which.min(perf.nodewise$error_1)]
lambda_2 <- grid[which.min(perf.nodewise$error_2)]
lambda_3 <- grid[which.min(perf.glasso$error)]
abline(v=log10(lambda_1), col=2, lty=2)
abline(v=log10(lambda_2), col=3, lty=4)
abline(v=log10(lambda_3), col=4, lty=2)

# optimal lambdas
lambda_1
lambda_2
lambda_3

#### cross-validation
# compute Error Matrix
set.seed(666)
grid <- 10 ^ seq(-0.7, -1.3, length = 50)
K <- 10 # number of folds
folds <- sample( rep(1:K, length = n) )
Error_1 <- c()
Error_2 <- c()
Error_3 <- c()
for (k in 1:K){
  perf.nodewise <- performance.nodewise.grid(data$X[folds!=k,], data$E_true, grid)
  perf.glasso <- performance.glasso.grid(data$X[folds!=k,], data$E_true, grid)
  Error_1 <- rbind(Error_1, perf.nodewise$error_1)
  Error_2 <- rbind(Error_2, perf.nodewise$error_2)
  Error_3 <- rbind(Error_3, perf.glasso$error)
}

# cv plot
plot.cv.error(Error_1, grid, col=2, ylim=c(0, 0.014)
              , xlab=expression(paste("log(",lambda,")")), ylab="mean error rate")
par(new=T)
plot.cv.error(Error_2, grid, col=3, ylim=c(0, 0.014)
              , xlab=expression(paste("log(",lambda,")")), ylab="mean error rate")
par(new=T)
plot.cv.error(Error_3, grid, col=4, ylim=c(0, 0.014)
              , xlab=expression(paste("log(",lambda,")")), ylab="mean error rate")
legend("topright",legend=c("node-wise 1","node-wise 2", "graphical"), col=c(2,3,4), pch="-", cex = 0.75)

# zoom-in plot
par(mfrow=c(1,3))
plot.cv.error(Error_1, grid, zoom = 2
              , xlab=expression(paste("log(",lambda,")")), ylab="mean error rate")
title(main = "node-wise 1")
plot.cv.error(Error_2, grid, zoom = 2
              , xlab=expression(paste("log(",lambda,")")), ylab="mean error rate")
title(main = "node-wise 2")
plot.cv.error(Error_3, grid, zoom = 6
              , xlab=expression(paste("log(",lambda,")")), ylab="mean error rate")
title(main = "graphical")

# optimal lambdas using cv
10^-1.045
10^-1.05
10^-0.96



### Repetition ###################################################################
# it takes about 10 min to run this code chunk
set.seed(666)
p <- 50
n <- 1000
prob <- 0.1
grid <- 10 ^ seq(0.5,-2.5, length = 50)
auc_1 <- c()
auc_2 <- c()
auc_3 <- c()
Error_1 <- c()
Error_2 <- c()
Error_3 <- c()
t = 1
while (t<=50){
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
boxplot(auc_1, auc_2, auc_3)
title("AUC")
boxplot(min_error_1, min_error_2, min_error_3)
title("Minumun Error Rate")



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
probs <- 10 ^ seq(-4,-0.3, length = 10)
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