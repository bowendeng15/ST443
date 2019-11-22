################################################################
### INSTALL PACKAGES NEEDED
################################################################
install.packages("Hmisc")
install.packages("car")
install.packages("leaps")

library("Hmisc")
library("car")
library(plyr)
###1. read and clean data
################################################################
bank0<-read.csv("bank-additional-full.csv",sep=";",dec=".")
table(is.na(bank0))
bank0$pdays[bank0$pdays==999]=-1 #999 means no previous contact, would make other value meaningless
head(bank0)
summary(bank0)

################################################################
###2. checking for collinearity
################################################################
M <- model.matrix(y~., data = bank0)[,-1]
dim(M)
corrM<- round(cor(M),3)
corrP<- rcorr(M)$P
# cormat : matrix of the correlation coefficients
# pmat : matrix of the correlation p-values
listedCorrMatrix <- function(corrmat, pmat) {
  ut <- upper.tri(corrmat)
  data.frame(
    row = rownames(corrmat)[row(corrmat)[ut]],
    column = rownames(corrmat)[col(corrmat)[ut]],
    corr  =(corrmat)[ut],
    p = pmat[ut]
  )
}
corrplist<- as.data.frame(listedCorrMatrix(corrM,corrP))
ord_corrplist<- arrange(corrplist, abs(corrplist$corr), decreasing=T)
## "loanunknown" and " housingunknown " have the exact correlation of 1 
nloan<-bank0[loan=="unknown",]
nhouse<- bank0[housing=="unknown",]

table(nhouse==nloan) ##"all  true" indicates that the clients whose condition of 
## personal loan is unknown also have no information related to 
## house loan

testloglm<- glm(y~., data=bank0, family = "binomial")
summary(testloglm)  

## In order to ensure little or no multicollinearity, 
## we first exclude the "loanunknown" form the set of model covariates
bank <- cbind(data.frame(M[,which(colnames(M)!="loanunknown")]), y)
head(bank)
testloglm1<- glm(y~., data=bank, family = "binomial")
summary(testloglm1) 

vif1<- as.data.frame(vif(testloglm1))
vif1_1<- as.data.frame(cbind(variable=row.names(vif1), vif1))
vif_ind<- arrange(vif1_1,vif(testloglm1), decreasing = T)


testloglm2<- update(testloglm1,.~.-nr.employed)
summary(testloglm2)
vif2<- as.data.frame(vif(testloglm2))


testloglm3<- update(testloglm2,.~.-emp.var.rate)
summary(testloglm3)
vif3<- as.data.frame(vif(testloglm3))
print(vif3)
table(vif3>10)

## We need to be aware of excluding the "loanunknown", "nr.employed" and "emp.var.rate" 
## from the covariates when building the logistic regression model to satisfy 
## the assumption of no or litte multicollinearity.


###############################################################################################
###TESTING LOGISTICS REGRESSION, LDA, QDA AND KNN ALL BASED ON THE VALIDATION SET APPROACH
################################################################################################
## 1. Generate the training and testing datasets
set.seed(1)
train<- sample(x=dim(bank)[1], size=dim(bank)[1]*2/3)

# Create testing and training data set
bank_test = bank[-train,] 
bank_train = bank[train,]
# Create testing data for Y
y_test = bank$y[-train]

## 2. logistic regression on the training data set
glm_fit0 = glm(
  y~.-nr.employed-emp.var.rate,
  data = bank_train,
  family = binomial,
)
summary(glm_fit0)

# generate the training dataset without outliers
cooksd<-cooks.distance(glm_fit0)
bank_cd<-cbind(bank_train,cooksd)
bank_train_df<- bank_cd[bank_cd$cooksd< 8/(dim(bank_train)[1]-2*53),-dim(bank_cd)[2]] #exclude the outliers
head(bank_train_df)

plot(cooksd, pch="*", cex=2, main="Influential Obs by Cook's distance")  # plot cook's distance
abline(h = 8/(dim(bank[train,])[1]-2*53), col="red")  # add cutoff line
# add labels
text(x=1:length(cooksd)+1, y=cooksd, labels=ifelse(cooksd>0.002,names(cooksd),""), col="red",cex=0.7)  
dev.copy(jpeg,"Influential Obs by Cook's Distance.jpg",height=8*72,width=8*72)
dev.off()

# building logistic model with clean training dataset
glm_fit <- glm(
  y~.-nr.employed-emp.var.rate,
  data = bank_train_df,
  family = binomial,
)
summary(glm_fit)

#Predicted probabilities for the testing data set 
glm_probs <- predict(glm_fit, bank_test, type = "response")

# Sample size for the testing data
dim(bank_test)
# Checking for default
contrasts(y)
# For predicted probabilities greater than 0.5, assign Y to be "yes"; otherwise assign Y to be "no"
glm_pred = rep("no", dim(bank_test)[1])
glm_pred[glm_probs > .5] = "yes"

# Confusion matrix
print(table(glm_pred, y_test))
# Proportation of make correct classification
mean(glm_pred==y_test)
# Misclassfication error rate
print(mean(glm_pred!=y_test))
                  # glm_pred is the predicted Y for testing data 
                  # and y_test is the true Y for testing data

## 3. Linear Discriminant Analysis
library(MASS) 
# Perform LDA on the traning data set 
lda_fit = lda(y~., data = bank, subset = train)
names(predict(lda_fit, bank_test))

lda_pred_posterior = predict(lda_fit, bank_test)$posterior
head(lda_pred_posterior)

lda_pred = predict(lda_fit, bank_test)$class
head(lda_pred)

# Confusion matrix
table(lda_pred, y_test)
# Misclassfication error rate
mean(lda_pred != y_test) ## sightly higher than the rate from logistic regression


## 4. Quadratic Discriminant Analysis
# QDA cannot be performed 

## 5. K-Nearest Neighbors
# Perform K-nearest neighbours on the traning data set
library(class)

#standardize the numerical variables
numBank<- bank0[,c(1,11:14,16:20)]
catBank<- bank0[,-c(1,11:14,16:20,21)]
sd(numBank$age)
stdnumBank<- scale(numBank)
stdBank<- cbind(stdnumBank,catBank,y)
stdBankM<- model.matrix(y~., data = stdBank)[,-1]
std_df<- cbind(data.frame(stdBankM),y)
head(std_df)

# decide the number of K to use 
set.seed(123)
ksample<- sample(x=dim(std_df)[1], size = 5000)
kdataset<- std_df[ksample,]
ktrain<- sample(5000, 5000*2/3, replace = F)
x<- seq(from=3,to=79,by=2)
y<- c()
for( k in x ){
  knn_pred<- knn(train = kdataset[ktrain, -ncol(kdataset)],
                 test = kdataset[-ktrain, -ncol(kdataset)],
                 cl = kdataset$y[ktrain],
                 k=k
                 )
  y = append(y, mean(knn_pred==kdataset$y[-ktrain]))
}
print(which.max(y))
plot(x,y,"b")
points(which.max(y), y[which.max(y)], col="red", cex=2, pch=20)

# Create training dataset
train_X <- std_df[train,-ncol(std_df)]
# Create testing data 
test_X <- std_df[-train,-ncol(std_df)]
# Create training data for Y
train_y <- std_df$y[train]

# Set k=13, which is selected above
knn_pred13 <- knn(train_X, test_X, train_y, k = 13)
#confusion matrix
table(knn_pred13, std_df$y[-train])
#missclassification rate
print(mean(knn_pred13 != std_df$y[-train])) #which is also higher than using logistic regression 
                                            #based on the validation set approach
############################################################################################
###Model Selection
############################################################################################
##1. Best Subset Selection
library(leaps)


#regfit_full <- regsubsets(y ~ ., data = bank,method = "exhaustive", nvmax = 20,really.big=T)
#cannot be performed 


##2. Forward Stepwise Selection
regfit_fwd <- regsubsets(y ~. , #to exempt linear dependency
                        data = bank,
                        nvmax = 52,  # 52 dummy variables in total
                        method = "forward") 

reg_summary <- summary(regfit_fwd)
head(reg_summary)
names(reg_summary)

par(mfrow=c(2,2),mar=c(4,4,1,1))
plot(reg_summary$rss, xlab="Number of Variables", ylab="RSS")

# Plot adjusted R2 vs Number of Variables
plot(reg_summary$adjr2, xlab="Number of Variables", ylab="Adjusted RSq", type = "l")

# which_max() function is used to identify the location of the maximum point of a vector
best_model_adjr2 <- which.max(reg_summary$adjr2)
print(best_model_adjr2)

# Plot a red dot to indicate the model with the largest adjusted R2
points(best_model_adjr2, reg_summary$adjr2[best_model_adjr2], col="red", cex=2, pch=20)

## In a similar fashion, we can plot Cp and BIC
plot(reg_summary$cp, xlab="Number of Variables", ylab="Cp", type="l")
best_model_cp = which.min(reg_summary$cp)
points(best_model_cp, reg_summary$cp[best_model_cp], col="red", cex=2, pch=20)
print(best_model_cp)

plot(reg_summary$bic, xlab="Number of Variables", ylab="BIC", type="l")
best_model_bic = which.min(reg_summary$bic)
points(best_model_bic, reg_summary$bic[best_model_bic], col="red", cex=2, pch=20)
print(best_model_bic)

## Check the coefficient estimates associated with models of size 21, indicated by using BIC.

print(names(coef(regfit_fwd, 21))[-1])




