library(MASS) # needed for studres
library(leaps)
library(car) # for vif
setwd("C:/Users/yuqi/Documents/R/wd")


bodyfatmendata=read.csv("bodyfatmen.csv")

lm.fit <- lm(density ~ age + weight + height + neck + chest + abdomen
             + hip + thigh + knee + ankle + biceps + forearm + wrist, data = bodyfatmendata)
#residual analysis
summary(lm.fit)
par(mfrow = c(2, 2))
qqnorm(rstandard(lm.fit))
qqnorm(lm.fit$residuals)

sres <- studres(lm.fit)
data.frame(sres)
plot(lm.fit$fitted.values, sres, main = "Studentized residuals")
# R-student
rstud <- rstudent(lm.fit)
data.frame(rstud)

plot(lm.fit$fitted.values, rstud, main = "R-student residuals")


X1 <- model.matrix(lm.fit) # get the X matrix
X <- X1[, -1] # remove the intercept column (column 1)

p <- length(lm.fit$coefficients)
n <- nrow(bodyfatmendata)
#Leverage, outliers & influential points

leverage.cutoff <- 2*p/n # Montgomery p. 213
cooks.cutoff <- qf(0.5, p, n - p, lower.tail = FALSE) # Montgomery p. 215
studres.cutoff <- qt(0.05/2, n - p, lower.tail = FALSE) # Montgomery p. 135
leverage.cutoff
cooks.cutoff
studres.cutoff
cooksd <- cooks.distance(lm.fit)
plot(cooksd, pch="*", main="Cooks distance measure of influence")
abline(h=cooks.cutoff)
plot(sres, pch="*", main="Studentized residual measure of outlier")
abline(h=studres.cutoff)
hatdiag <- hatvalues(lm.fit)
plot(hatdiag, pch="*", main="Hat diagonal measure of leverage")
abline(h=leverage.cutoff)

#remove 39th data in bodyfatmendata
bodyfatmendata <- bodyfatmendata[-c(39), ]
lm.fit <- lm(density ~ age + weight + height + neck + chest + abdomen
             + hip + thigh + knee + ankle + biceps + forearm + wrist, data = bodyfatmendata)
summary(lm.fit)
p <- length(lm.fit$coefficients)
n <- nrow(bodyfatmendata)

#multicollinearity
WtW <- cor(X)
WtW

# VIF
vif(lm.fit)
# or diagonals in
solve(WtW)

# Condition number
eval <- eigen(WtW)
eval
condition_number <- max(eval$values) / min(eval$values)
condition_number # (100<k<1000 so moderate to strong multicollinearity  see p. 298 Montgomery)

#all possible regression
predict.regsubsets =function (object ,newdata ,id ,...){
 form=as.formula (object$call[[2]])
 mat=model.matrix (form, newdata )
 coefi =coef(object, id=id)
 xvars =names (coefi )
 mat[,xvars ]%*%coefi
}


n <- nrow(bodyfatmendata)
set.seed(2)
train_valid <- sample(1:n, 85*n/100)
test <- (-train_valid)
bodyfatmendata.train_valid <- bodyfatmendata[train_valid, ]
bodyfatmendata.test <- bodyfatmendata[test, ]
#valid <- sample(1:nrow(bodyfatmendata.train_valid), 15*nrow(bodyfatmendata.train_valid)/85)
#bodyfatmendata.valid <- bodyfatmendata[valid, ]
#bodyfatmendata.train <- bodyfatmendata[-valid, ]

#all possible regression on data
regfit.full = regsubsets (density~.,data=bodyfatmendata.train_valid ,nvmax =13)
reg.summary = summary(regfit.full)
par(mfrow =c(2,2))
plot(reg.summary$adjr2 , col=ifelse(reg.summary$adjr2==max(reg.summary$adjr2), "red", "black"), xlab =" Number of Variables ",
     ylab="Adjusted RSq",type="p")
plot(reg.summary$cp , col=ifelse(reg.summary$cp==min(reg.summary$cp), "red", "black"), xlab =" Number of Variables ",ylab="Cp",
     type="p")
plot(reg.summary$bic , col=ifelse(reg.summary$bic==min(reg.summary$bic), "red", "black"), xlab=" Number of Variables ",ylab="BIC",
     type="p")
plot(reg.summary$rss ,xlab=" Number of Variables ",ylab="RSS",
     type="p")

par(mfrow =c(1,2))
mse <- replicate(13, 0)
r2.adj <- replicate(13, 0)
for(i in 1:13) {
   pred=predict (regfit.full, bodyfatmendata.test, id=i, return = TRUE)
   mse[i]=mean( (bodyfatmendata.test$density-pred)^2)
   p <- i+1
   n <- nrow(bodyfatmendata.test)
   ss.res <- sum( (bodyfatmendata.test$density - pred)^2)
   ss.tot <- sum( (bodyfatmendata.test$density 
                   - mean(bodyfatmendata.test$density)) ^ 2)
   r2.adj[i] <- 1 - (ss.res / (n - p)) / (ss.tot / (n - 1))
}
plot(mse, col=ifelse(mse==min(mse), "red", "black"), xlab=" Number of Variables ",ylab="MSE")
plot(r2.adj, col=ifelse(r2.adj==max(r2.adj), "red", "black"), xlab=" Number of Variables ",ylab="Adjusted RSq")

#K-folds CV repeating N times
N=50
k=5
sum.mean.cv.errors <- replicate(13, 0)
sum.mean.cv.r2.adj <- replicate(13, 0)
sum.mean.cv.cp <- replicate(13, 0)
sum.mean.cv.bic <- replicate(13, 0)

for(aa in 1:N){
set.seed (aa)
folds=sample (1:k,nrow(bodyfatmendata.train_valid),replace =TRUE)
cv.mse =matrix (NA ,k,13, dimnames =list(NULL , paste (1:13) ))
cv.bic =matrix (NA ,k,13, dimnames =list(NULL , paste (1:13) ))
cv.cp =matrix (NA ,k,13, dimnames =list(NULL , paste (1:13) ))
cv.r2.adj =matrix (NA ,k,13, dimnames =list(NULL , paste (1:13) ))

for(j in 1:k){
   best.fit =regsubsets (density~.,data=bodyfatmendata.train_valid [folds !=j,],
                         nvmax =13)
   reg_summary = summary(best.fit)
   for(i in 1:13) {
      pred=predict (best.fit, bodyfatmendata.train_valid [folds ==j,], id=i, return = TRUE)
      cv.mse[j,i]=mean( (bodyfatmendata.train_valid$density[folds ==j]-pred)^2)
      p <- i+1
      n <- nrow(bodyfatmendata.train_valid[folds ==j,])
      ss.res <- sum( (bodyfatmendata.train_valid$density[folds ==j] - pred)^2)
      ss.tot <- sum( (bodyfatmendata.train_valid$density[folds ==j] 
                     - mean(bodyfatmendata.train_valid$density[folds ==j])) ^ 2)
      cv.r2.adj[j,i] <- 1 - (ss.res / (n - p)) / (ss.tot / (n - 1))
      cv.bic[j,i] <- reg_summary$bic[i]
      cv.cp[j,i] <- reg_summary$cp[i]
   }
}

mean.cv.errors <- apply(cv.mse ,2, mean)
mean.cv.r2.adj <- apply(cv.r2.adj ,2, mean)
mean.cv.cp <- apply(cv.cp ,2, mean)
mean.cv.bic <- apply(cv.bic ,2, mean)

sum.mean.cv.errors = sum.mean.cv.errors + mean.cv.errors
sum.mean.cv.r2.adj = sum.mean.cv.r2.adj + mean.cv.r2.adj
sum.mean.cv.cp = sum.mean.cv.cp + mean.cv.cp
sum.mean.cv.bic = sum.mean.cv.bic + mean.cv.bic
}
par(mfrow =c(2,2))
plot(sum.mean.cv.errors/N, main = "MSE on validation data", xlab = "No. of Variables", ylab = "MSE")
plot(sum.mean.cv.r2.adj/N, main = "AdjR2 on validation data", xlab = "No. of Variables", ylab = "AdjR2")
plot(sum.mean.cv.cp/N, main = "Cp on training data", xlab = "No. of Variables", ylab = "Cp")
plot(sum.mean.cv.bic/N, main = "BIC on training data", xlab = "No. of Variables", ylab = "BIC")

reg.best=regsubsets (density~.,data=bodyfatmendata.train_valid , nvmax =4)
coef(reg.best ,4)
coef(reg.best ,3)

#evaluate on test data
pred=predict (reg.best, bodyfatmendata.test, id=3, return = TRUE)
test.mse=mean( (bodyfatmendata.test$density-pred)^2)
test.mse
p <- 3+1
n <- nrow(bodyfatmendata.test)
ss.res <- sum( (bodyfatmendata.test$density - pred)^2)
ss.tot <- sum( (bodyfatmendata.test$density 
                - mean(bodyfatmendata.test$density)) ^ 2)
test.r2.adj <- 1 - (ss.res / (n - p)) / (ss.tot / (n - 1))
test.r2.adj

###################
#### Bootstrap ####
###################
# Resample data 3000 times. The intercept value behaves as a normal distribution, as expected.

# 4 variables model
lm.fit4 <- lm(density ~ height + chest + abdomen + wrist, data = bodyfatmendata.train_valid)
system.time(bo1<-Boot(lm.fit4,R=3000))
summary(bo1)
confint(bo1)
plot(bo1,main="Values of intercept")
inf.fit= influence.measures(lm.fit4)
summary(inf.fit)
# 3 variables model
lm.fit3 <- lm(density ~ height + abdomen + wrist, data = bodyfatmendata.train_valid)
system.time(bo2<-Boot(lm.fit3,R=3000))
summary(bo2)
confint(bo2)
plot(bo2,main="Values of intercept")
inf.fit= influence.measures(lm.fit3)
summary(inf.fit)
