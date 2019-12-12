### IMPORT PACKAGES NEEDED ############################################################################################
library(ROCR) # performance prediction
library(MASS) # lda stepAIC
library(class) # knn
library(glmnet) 
library(glmnetUtils) # glmnet.formula
library(leaps) # regsubset
library(gbm) # gbm
library(caret) # trainControl
library(e1071)
library(randomForest)
library(boot)

#' @title: general performance and plot ROC
perf = function(prob, y, ...){
  pred_test <- rep("no", length(y))
  pred_test[prob > .5] = "yes"
  mat <- table(pred_test, y_test)
  cat("\nconfusion matrix----------------------\n")
  print( mat )
  cat("\naccuracy----------------------\n")
  print( mean(pred_test==y_test) )
  cat("\nrecall----------------------\n")
  print( mat[2,2]/sum(mat[,2]) )
  pred_glm <- prediction(prob, y)
  cat("\nAUC----------------------\n")
  print( performance(pred_glm, "auc")@y.values[[1]] )
  #ROC
  perf_glm <- performance(pred_glm, measure = "tpr", x.measure = "fpr")   
  plot(perf_glm, ...)    
  abline(0, 1, col="grey", lty=2)
}



### 1. Data ############################################################################################
setwd("./")
bank_df<-read.csv("bank-additional-full.csv",sep=";",dec=".")
bank_df$pdays = 1/(1+bank_df$pdays) # transform to values in [0,1]
bank_df<- bank_df[,which(colnames(bank_df)!= "duration")]



### 2. Preprocessing ############################################################################################
RNGkind(sample.kind = "Rounding")
set.seed(1)
train<- sample(x=dim(bank_df)[1], size=dim(bank_df)[1]*2/3)

# standardize the numerical variables
numBank<- bank_df[,c(1,11:13,15:19)] #numercial column
catBank<- bank_df[,-c(1,11:13,15:19,20)] #categorcial column

scale<-apply(numBank[train,], 2, sd)
center<-apply(numBank[train,], 2, mean)

center["campaign"]<- 1
center["previous"]<- 0

tmp <- scale(numBank,center = center,scale = scale) #scale
tmp <- cbind(tmp, catBank, y=bank_df$y) # bind y
tmp <- model.matrix(y~., data = tmp)[,-1] # design matrix

cor(tmp[,c("loanunknown", "housingunknown")]) # cor=1, collinear
bank <- cbind(data.frame(tmp[,which(colnames(tmp)!="loanunknown")]), y=bank_df$y)

# train-test split
bank_test = bank[-train,] 
bank_train = bank[train,]
y_test = bank$y[-train]

# oversampling
bank_train_yes<- bank_train[bank_train$y=="yes",]
set.seed(2)
resample_idx<- sample(nrow(bank_train_yes),sum(bank_train$y=="no"),replace = T)
bank_train_res<- rbind(bank_train[bank_train$y=="no",],bank_train_yes[resample_idx,]) #resampled training data


### 3. Try LR, LDA, KNN and RF ############################################################

## a) Building logistic model 
glm_fit <- glm(y~., data = bank_train, family = "binomial")

prob_test <- predict(glm_fit, bank_test, type = "response")
perf(prob_test, bank_test$y)

#Using resampled dataset
glm_fit_res <- glm(y~., data = bank_train_res, family = "binomial")

prob_test <- predict(glm_fit_res, bank_test, type = "response")
perf(prob_test, bank_test$y)


## b) Linear Discriminant Analysis
lda_fit = lda(y~., data = bank_train)

prob_test = predict(lda_fit, bank_test)$posterior[,2]
perf(prob_test, bank_test$y)

#Using resampled dataset
lda_fit_res = lda(y~., data = bank_train_res)

prob_test = predict(lda_fit_res, bank_test)$posterior[,2]
perf(prob_test, bank_test$y)


## c) K-Nearest Neighbors ( k=15 )
knn_pred <- knn(bank_train_res[,-ncol(bank_train_res)], bank_test[,-ncol(bank_test)], bank_train_res$y, k = 15, prob = T)

knn_prob <- attr(knn_pred,"prob")
prob_test <- ifelse(knn_pred == "no", 1-knn_prob, knn_prob)  
perf(prob_test, bank_test$y)


## d) Random Forest
set.seed(135)
rf.bank <-randomForest(y~., data=bank_train_res, importance=TRUE, ntree = 500)
rf.bank$importance

prob_test <- predict(rf.bank, bank_test,type= "prob")[,2] 
perf(prob_test, bank_test$y)

### 4. LR and Tuning ############################################################################################

# Detecting outliers by calculating cook's distance
cooksd<-cooks.distance(glm_fit_res)
tmp <- cbind(bank_train_res,cooksd)
bank_train_df<- tmp[tmp$cooksd< 8/(dim(bank_train_res)[1]-2*53),-dim(tmp)[2]] #exclude the outliers

nrow(bank_train_res)-nrow(bank_train_df)

## a) Relaxed Lasso ( 5-fold cross validation )
set.seed(2)
cv.lasso <-cv.glmnet(as.matrix(bank_train_df[,-ncol(bank_train_res)]), bank_train_df$y
                     , nfolds = 5, family = "binomial", type.measure = "class")
plot(cv.lasso, ylim=c(0.23, 0.27))
col_lasso = rownames(coef(cv.lasso))[which(coef(cv.lasso)!=0)][-1]
length(col_lasso)
col_lasso = append(col_lasso, "y")

# fit a new LR model
lr_lasso <- glm(y ~ ., data = bank_train_df[,col_lasso], family = "binomial")

prob_test <- predict(lr_lasso, newdata = bank_test[,col_lasso], type = "response")
perf(prob_test, bank_test$y)

## b) Stepwise Selection
## Forward Stepwise Selection by regsubset()
regfit_fwd = regsubsets(y ~ .
                        , data = bank_train_df
                        , really.big = T
                        , nvmax = 51
                        , method = "forward"
)
reg_summary = summary(regfit_fwd)


plot(reg_summary$bic, xlab="Number of Variables", ylab="BIC", type="l")
best_model_bic = which.min(reg_summary$bic)
print(best_model_bic) #23
points(best_model_bic, reg_summary$bic[best_model_bic], col="red", cex=2, pch=20)

col_fwd = names(coef(regfit_fwd, id=which.min(reg_summary$bic)))[-1]
col_fwd = append(col_fwd, "y")
lr_fwd <- glm(y ~ ., data = bank_train_df[,col_fwd], family = "binomial")

prob_test <- predict(lr_fwd, newdata = bank_test[,col_fwd], type = "response")
perf(prob_test, bank_test$y)

## Backward Stepwise Selection by regsubset()
regfit_bwd = regsubsets(y ~ .
                        , data = bank_train_df
                        , really.big = T
                        , nvmax = 51
                        , method = "backward"
)
bwd_summary = summary(regfit_bwd)

plot(bwd_summary$bic, xlab="Number of Variables", ylab="BIC", type="l")
best_model_bic_bwd = which.min(bwd_summary$bic)
print(best_model_bic_bwd)#24
points(best_model_bic_bwd, bwd_summary$bic[best_model_bic_bwd], col="red", cex=2, pch=20)

col_bwd = names(coef(regfit_bwd, id=which.min(bwd_summary$bic)))[-1]
col_bwd = append(col_bwd, "y")
lr_bwd <- glm(y ~ ., data = bank_train_df[,col_bwd], family = "binomial")

prob_test <- predict(lr_bwd, newdata = bank_test[,col_bwd], type = "response")
perf(prob_test, bank_test$y)

## Backward Stepwise Selection by stepAIC() DONT RUN IT!!!!!!!!!!!!!!!!!!!!!TAKE LOOOOONG TIME
lr_full <- glm(
  y~.,
  data = bank_train_df,
  family = binomial,
)
summary(lr_full)
step1<- stepAIC(lr_full,direction = "backward",k=log(nrow(bank_train_df)),trace = F)
## Setting k to make BIC play the role
step1$anova
##Final Model:
## y ~ campaign + pdays + previous + emp.var.rate + cons.price.idx + 
##  euribor3m + jobretired + jobstudent + maritalsingle + educationuniversity.degree + 
##  defaultunknown + contacttelephone + monthaug + monthdec + 
##  monthjul + monthjun + monthmar + monthmay + monthnov + monthoct + 
##  monthsep + day_of_weekmon + day_of_weektue + poutcomenonexistent + 
##  poutcomesuccess

col_bwd_BIC<- c("campaign" ,"pdays", "previous", "emp.var.rate" , "cons.price.idx" , 
  "euribor3m" ,"jobretired", "jobstudent", "maritalsingle",  "educationuniversity.degree",
  "defaultunknown" , "contacttelephone" , "monthaug" ,"monthdec","monthjul", "monthjun" ,
  "monthmar" , "monthmay" , "monthnov" ,"monthoct", "monthsep", "day_of_weekmon" , 
  "day_of_weektue", "poutcomenonexistent" , "poutcomesuccess","y")
length(col_bwd_BIC)

lr_bwd_BIC<- glm(y ~ ., data = bank_train_df[,col_bwd_BIC], family = "binomial")

prob_test <- predict(lr_bwd_BIC, newdata = bank_test[,col_bwd_BIC], type = "response")
perf(prob_test, bank_test$y)


# ROC plot for all Four
par(mfrow=c(1,1))
prob_test <- predict(lr_lasso, newdata = bank_test[,col_lasso], type = "response")
perf(prob_test, bank_test$y, col="red")
par(new=T)
prob_test <- predict(lr_fwd, newdata = bank_test[,col_fwd], type = "response")
perf(prob_test, bank_test$y, col="blue")
par(new=T)
prob_test <- predict(lr_bwd, newdata = bank_test[,col_bwd], type = "response")
perf(prob_test, bank_test$y, col="orange")
par(new=T)
prob_test <- predict(lr_bwd_BIC, newdata = bank_test[,col_bwd_BIC], type = "response")
perf(prob_test, bank_test$y, col="pink")
text(x=0.8,y=0.3,"AUC_lasso=0.794",cex =0.8,col="red")
text(x=0.8,y=0.2,"AUC_reg_fwd=0.796",cex =0.8,col="blue")
text(x=0.8,y=0.1,"AUC_reg_bwd=0.797",cex =0.8,col="orange")
text(x=0.8,y=0.0,"AUC_stepAIC=0.798",cex =0.8,col="pink")


## Polynomial Regression and lasso
set.seed(1)
cv.poly <-cv.glmnet(y ~ .
                     + I(age^2)
                     + I(campaign^2)
                     + I(pdays^2)
                     + I(previous^2)
                     + I(emp.var.rate^2)
                     + I(cons.price.idx^2)
                     + I(cons.conf.idx^2)
                     + I(euribor3m^2)
                     + I(nr.employed^2)
                    , data = bank_train_df, nfolds = 5
                    ,family = "binomial", type.measure = "class")
plot(cv.poly, ylim=c(0.24, 0.26))
col_poly_lasso = rownames(coef(cv.poly))[which(coef(cv.poly)!=0)][-1]
length(col_poly_lasso)

( formula = paste("y~", paste(col_poly_lasso, collapse = "+") ) )
lr_poly_lasso <- glm(formula, data = bank_train_df, family = "binomial")

prob_test <- predict(lr_poly_lasso, newdata = bank_test, type = "response")
perf(prob_test, bank_test$y, col="red")


### 5. Adaboost and Tuning #############################################################################

#### Cross Validation for Tuning Parameters 
set.seed (123)
fitControl = trainControl(method = "repeatedcv", number = 5, repeats = 1)
gbmGrid = expand.grid(
  interaction.depth = c(1, 2, 3, 4)
  , n.trees = 500
  , shrinkage = 0.1
  , n.minobsinnode = c(5, 10, 15, 20)
)
# it takes 30 mins to run
gbmFit = train(y~., data=bank_train_res, distribution = 'adaboost', verbose=FALSE
               , method = 'gbm', trControl=fitControl, tuneGrid=gbmGrid)
gbmFit$results[order(gbmFit$results$Accuracy, decreasing=TRUE),]
gbmFit$bestTune
##n.trees interaction.depth shrinkage n.minobsinnode
##15     500              4       0.1             10
# prediction
prob_test <- predict(gbmFit, bank_test,type= "prob")[,2]
perf(prob_test, bank_test$y, col="red")



### 6. Final Model ############################################################
C = summary(lr_poly_lasso)$coefficient
nrow(C)-1
