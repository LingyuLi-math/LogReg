##  clear
rm(list = ls())

################################################################################
# Download R package 
################################################################################
library(tidyverse)
library(caret)        # ten-fold cross-validation 
library(glmnet)
library("l0ara")
library(ggplot2)
library(caTools)
library(magrittr)
library(ROCR)
library(pROC)
library(glmnetUtils)
library(ncpen)
library(stargazer)    # transfer into latex 
library(broom)        # save regression results
library(ncvreg)
library("l0ara")
# install.packages("caret")
library(plyr)
library(pROC)


################################################################################
# Load data  
################################################################################
x = read.table(".\\GSE59491_15\\Data\\GSE59491_scale_DE.txt", header = T, check.names = FALSE)
data <- data.frame(t(x))


## Data split: training data + testing data
set.seed(1234)
training.samples <- data$Lable %>% createDataPartition(p = 0.7, list = FALSE)
train.data  <- data[training.samples, ]
test.data <- data[-training.samples, ]


x <- model.matrix(Lable ~., train.data)[,-1]   # delete Lable
y <- train.data$Lable                          
x.test <- model.matrix(Lable ~., test.data)[,-1] 
y.test <- test.data$Lable
# sum(y == 0)


################################################################################
# pathways 
################################################################################
setwd(".\\GSE59491_15\\Data37")
set.seed(123)

ridge.fit <- glmnet(x, y, family="binomial", alpha = 0, lambda = NULL) 
jpeg(file = "ridge_fit.jpg")
# postscript("ridge_fit.eps")
plot(ridge.fit, xvar = "lambda")
dev.off()
cv.ridge = cv.glmnet(x, y, alpha = 0, family = "binomial", nfolds = 10, type.measure = "class")
# jpeg(file = "cv_ridge.jpg")
postscript("cv_ridge.eps")
plot(cv.ridge)
dev.off()

lasso.fit <- glmnet(x, y, family="binomial", alpha = 1, lambda = NULL) 
jpeg(file = "lasso_fit.jpg")
# postscript("lasso_fit.eps")
plot(lasso.fit, xvar = "lambda")
dev.off()
cv.lasso = cv.glmnet(x, y, alpha = 1, family = "binomial", nfolds = 10, type.measure = "class")
jpeg(file = "cv_lasso.jpg")
# postscript("cv_lasso.eps")
plot(cv.lasso)
dev.off()

elastic.fit <- glmnet(x, y, family="binomial", alpha = 0.5, lambda = NULL) 
jpeg(file = "elastic_fit.jpg")
# postscript("elastic_fit.eps")
plot(elastic.fit, xvar = "lambda")
dev.off()
cv.elastic = cv.glmnet(x, y, alpha = 0.5, family = "binomial", nfolds = 10, type.measure = "class")
jpeg(file = "cv_elastic.jpg")
# postscript("cv_elastic.eps")
plot(cv.elastic)
dev.off()

cv.Bridge <- cv.ncpen(y.vec=y, x.mat=x, family="binomial", penalty="mbridge")
jpeg(file = "cv_Bridge.jpg")
# postscript("cv_Bridge.eps")
plot(cv.Bridge, type = "rmse", log.scale = T)
dev.off()
Bridge.model <- ncpen(y.vec = y, x.mat = x, family="binomial", penalty="mbridge")
jpeg(file = "Bridge_fit.jpg")
# postscript("Bridge_fit.eps")
plot(Bridge.model)
dev.off()

lam <- seq(1,0.05,-0.05)
cv.l0 <- cv.l0ara(x, y, family="logit", lam, measure = "class")
lambda.min <- cv.l0$lam.min
jpeg(file = "cv_l0.jpg")
# postscript("cv_l0.eps")
plot(cv.l0, col = 2)
dev.off()

l0.fit <- l0ara(x, y, family = "logit", 0.4)
jpeg(file = "l0_fit.jpg")
# postscript("l0_fit.eps")
plot(l0.fit, auc = F, split = F, col = 4) # Non-locally convex regions of solution path
dev.off()



## using 'X', 'y', not 'x', 'y' 
X <- model.matrix(Lable ~., train.data)[,-1]   # delete Lable
y <- train.data$Lable                          
x.test <- model.matrix(Lable ~., test.data)[,-1] 
y.test <- test.data$Lable

cv.SCAD <- cv.ncvreg(X, y, family ="binomial", penalty="SCAD") 
jpeg(file = "cv_SCAD.jpg")
# postscript("cv_SCAD.eps")
plot(cv.SCAD)
dev.off()
SCAD.fit <- ncvreg(X, y, family ="binomial", penalty="SCAD")
jpeg(file = "SCAD_fit.jpg")
# postscript("SCAD_fit.eps")
plot(SCAD.fit)
dev.off()

cv.MCP <- cv.ncvreg(X, y, family ="binomial", penalty="MCP")
jpeg(file = "cv_MCP.jpg")
# postscript("cv_MCP.eps")
plot(cv.MCP)
dev.off()
MCP.fit <- ncvreg(X, y, family ="binomial", penalty="MCP")
jpeg(file = "MCP_fit.jpg")
# postscript("MCP_fit.eps")
plot(MCP.fit)
dev.off()


################################################################################
# glmnet: Ridge penalty 
################################################################################
## save prediction results
coef_ridge <- matrix()       # save coefficients
pred_ridge <- matrix()       # save prediction


for(i in 1:30){
  
  ## set seeds
  set.seed(i)
  
  ## Cross-validation to obtain model parameters
  cv.ridge = cv.glmnet(x, y, alpha = 0, family = "binomial", nfolds = 10, type.measure = "class")
  lambda.min <- cv.ridge$lambda.min
  lambda.1se <- cv.ridge$lambda.1se
  print("***lambda.min、lambda.1se***")
  print(lambda.min)
  print(lambda.1se)
  
  ## fit
  ridge.model <- glmnet(x, y, alpha = 0, family = "binomial", lambda = cv.ridge$lambda.1se)
  coef <- as.matrix(coef(ridge.model))
  tcross <- rep(i, length(coef))                # i is the number of times the cycle crosses, totally K times.
  step_ridge <- data.frame(cbind(coef, tcross))
  coef_ridge <- cbind(coef_ridge, step_ridge)   # temp is merged with pred by row
  
  ## predict
  p <- predict(ridge.model, newx = x.test, type = "response")
  kcross <- rep(i, length(p)) 
  temp_ridge <- data.frame(cbind(y.test, p, kcross))
  pred_ridge <- cbind(pred_ridge,temp_ridge)    # temp is merged with pred by row
  
  print(paste("The：",i)) 
}

## set save oathway
setwd(".\\GSE59491_15\\Data37")
# write.csv(coef_ridge, file = "ridge\\coef_ridge.csv")
# write.csv(pred_ridge, file = "ridge\\pred_ridge.csv")

## define 'my_pred' function, Input: p_ridge
my_pred <- function(x){
  p_ridge1 <- x[,-2]
  p_ridge2 <- matrix(data=0, nrow = dim(p_ridge1)[1], ncol = 1, byrow = FALSE, dimnames=list(c(as.character(x[,1])),c("prob")))
  for (j in 1:10){
    p_ridge2 <- p_ridge2[] + p_ridge1[,3*j]
  }
  p_ridge3 <- p_ridge2/10
  return(p_ridge3)
}


## compute average 
p_ridge <- read.table("ridge\\pred_ridge.csv", header=TRUE, sep = ',')
pred_ridge <- cbind(y.test,my_pred(p_ridge))
# write.csv(pred_ridge, file = "ridge\\pred_ridge0.csv")

## ROC curve
# jpeg(file = "ridge\\pAUC_ridge.jpg")
plot.roc(pred_ridge[,1], pred_ridge[,2], print.auc=T, main="pAUC")
# dev.off()

## Performance index
predict = ifelse(pred_ridge[,2] > 0.5, 1, 0)
predict_value = predict 
true_value = pred_ridge[,1]
error = predict_value-true_value

## Calculate model accuracy, precision, recall rate (Recall), sensitivity, true positive rate (TPR), F-measure and confusion matrix.
## 'Precision' calculates the proportion of "items that should be retrieved (TP)" among all retrieved items (TP+FP);
## 'Recall' calculates the ratio of all retrieved items (TP) to all "items that should be retrieved (TP+FN)"
## 'Precision' is how many of the retrieved items are accurate. 'Recall' is how many of all accurate items have been retrieved.
accuracy = (nrow(data)-sum(abs(error)))/nrow(data) 
precision = sum(true_value & predict_value)/sum(predict_value)  # The true value and predicted value are all 1, and the predicted value is all 1.
recall = sum(predict_value & true_value)/sum(true_value)        # The true value and the predicted value are all 1 / the true value is all 1.
## 'precision' and 'recall' sometimes appear to be contradictory, so need to be considered comprehensively. The most common method is 'F-measure'.
F_measure= 2*precision*recall/(precision+recall)                # 'F-measure' is the weighted harmonic average of 'Precision' and 'Recall'.
specificity = ((nrow(data)-sum(abs(error)))-sum(predict_value & true_value))/(nrow(data)-sum(true_value)) 
### AUC
pred <- prediction(pred_ridge[,2], pred_ridge[,1])
auc <- performance(pred,'auc')
auc <- unlist(slot(auc,'y.values'))

## save Confusion matrix
# table(true_value, predict_value) %>% as.matrix(table) %>% write.csv(file = "ridge\\table.csv")     


## Output results
result <- matrix(0, 7, 6)
colnames(result) <- c("accuracy", "precision", "recall", "F_measure", "specificity", "AUC")
rownames(result) <- c("ridge", "lasso", "elatic net", "L1/2", "L0", "SCAD", "MCP")

i <- 1
result[i,1] <- accuracy
result[i,2] <- precision
result[i,3] <- recall
result[i,4] <- F_measure
result[i,5] <- specificity
result[i,6] <- auc
# result %>% write.csv(file = "result\\result.csv") 


################################################################################
# glmnet: Lasso penalty
################################################################################
## save prediction results
coef_lasso <- matrix()       # save coefficients
pred_lasso <- matrix()       # save prediction

for(i in 1:30){
  
  ## set seeds
  set.seed(i)
  
  ## Cross-validation to obtain model parameters
  cv.lasso = cv.glmnet(x, y, alpha = 1, family = "binomial", nfolds = 10, type.measure = "class" )
  lambda.min <- cv.lasso$lambda.min
  lambda.1se <- cv.lasso$lambda.1se
  print("***lambda.min、lambda.1se***")
  print(lambda.min)
  print(lambda.1se)
  
  ## fit
  lasso.model <- glmnet(x, y, alpha = 1, family = "binomial", lambda = cv.lasso$lambda.1se)
  coef <- as.matrix(coef(lasso.model))
  tcross <- rep(i, length(coef))                 
  step_lasso <- data.frame(cbind(coef, tcross))
  coef_lasso <- cbind(coef_lasso, step_lasso)    
  
  ## predict
  p <- predict(lasso.model, newx = x.test, type = "response")
  kcross <- rep(i, length(p)) 
  temp_lasso <- data.frame(cbind(y.test, p, kcross))  
  pred_lasso <- cbind(pred_lasso,temp_lasso)    
  
  print(paste("The：",i)) 
}

## set save pathway
# write.csv(coef_lasso, file = "lasso\\coef_lasso.csv")
# write.csv(pred_lasso, file = "lasso\\pred_lasso.csv")

## average value
p_lasso <- read.table("lasso\\pred_lasso.csv", header=TRUE, sep = ',')
pred_lasso <- cbind(y.test,my_pred(p_lasso))
# write.csv(pred_lasso, file = "lasso\\pred_lasso0.csv")

## ROC curve
# jpeg(file = "lasso\\pAUC_lasso.jpg")
plot.roc(pred_lasso[,1], pred_lasso[,2], print.auc=T, main="pAUC")
# dev.off()

## Performance index
predict = ifelse(pred_lasso[,2] > 0.5, 1, 0)
predict_value = predict 
true_value = pred_lasso[,1]
error = predict_value-true_value

accuracy = (nrow(data)-sum(abs(error)))/nrow(data) 
precision = sum(true_value & predict_value)/sum(predict_value)   
recall = sum(predict_value & true_value)/sum(true_value)         
F_measure= 2*precision*recall/(precision+recall)     
specificity = ((nrow(data)-sum(abs(error)))-sum(predict_value & true_value))/(nrow(data)-sum(true_value)) 
## AUC
pred <- prediction(pred_lasso[,2], pred_lasso[,1])
auc <- performance(pred,'auc')
auc <- unlist(slot(auc,'y.values'))

## Confusion matrix
# table(true_value, predict_value) %>% as.matrix(table) %>% write.csv(file = "lasso\\table.csv")     

## Output results
i <- 2
result[i,1] <- accuracy
result[i,2] <- precision
result[i,3] <- recall
result[i,4] <- F_measure
result[i,5] <- specificity
result[i,6] <- auc
# result %>% write.csv(file = "result\\result.csv") 


################################################################################
# glmnet: Elastic Net penalty 
################################################################################
## save prediction results
coef_elastic <- matrix()       # save coefficients
pred_elastic <- matrix()       # save prediction

for(i in 1:30){
  
  ## set seed
  set.seed(i)
  
  ## Cross-validation to obtain model parameters
  cv.elastic = cv.glmnet(x, y, alpha = 0.5, family = "binomial", nfolds = 10, type.measure = "class")
  lambda.min <- cv.elastic$lambda.min
  lambda.1se <- cv.elastic$lambda.1se
  print("***lambda.min、lambda.1se***")
  print(lambda.min)
  print(lambda.1se)
  
  ## fit model
  elastic.model <- glmnet(x, y, alpha = 0.5, family = "binomial", lambda = cv.elastic$lambda.1se)
  coef <- as.matrix(coef(elastic.model))
  tcross <- rep(i, length(coef))              
  step_elastic <- data.frame(cbind(coef, tcross))
  coef_elastic <- cbind(coef_elastic, step_elastic)   
  
  ## predict
  p <- predict(elastic.model, newx = x.test, type = "response")
  kcross <- rep(i, length(p)) 
  temp_elastic <- data.frame(cbind(y.test, p, kcross))
  
  pred_elastic <- cbind(pred_elastic,temp_elastic)   
  
  print(paste("The：",i)) 
}

## set save pathway
# write.csv(coef_elastic, file = "elastic_net\\coef_elastic.csv")
# write.csv(pred_elastic, file = "elastic_net\\pred_elastic.csv")

## average value
p_elastic <- read.table("elastic_net\\pred_elastic.csv", header=TRUE, sep = ',')
pred_elastic <- cbind(y.test,my_pred(p_elastic))
# write.csv(pred_elastic, file = "elastic_net\\pred_elastic0.csv")

## ROC curve
# jpeg(file = "elastic_net\\pAUC_elastic.jpg")
plot.roc(pred_elastic[,1], pred_elastic[,2], print.auc=T, main="pAUC")
# dev.off()

## Performance index
predict = ifelse(pred_elastic[,2] > 0.5, 1, 0)
predict_value = predict 
true_value = pred_elastic[,1]
error = predict_value-true_value

accuracy = (nrow(data)-sum(abs(error)))/nrow(data) 
precision = sum(true_value & predict_value)/sum(predict_value)   
recall = sum(predict_value & true_value)/sum(true_value)         
F_measure= 2*precision*recall/(precision+recall)     
specificity = ((nrow(data)-sum(abs(error)))-sum(predict_value & true_value))/(nrow(data)-sum(true_value)) 
## AUC
pred <- prediction(pred_elastic[,2], pred_elastic[,1])
auc <- performance(pred,'auc')
auc <- unlist(slot(auc,'y.values'))

## Confusion matrix 
# table(true_value, predict_value) %>% as.matrix(table) %>% write.csv(file = "elastic_net\\table.csv")     

## Output results
i <- 3
result[i,1] <- accuracy
result[i,2] <- precision
result[i,3] <- recall
result[i,4] <- F_measure
result[i,5] <- specificity
result[i,6] <- auc
# result %>% write.csv(file = "result\\result.csv") 


################################################################################
# ncpen: Bridge 1/2 penalty 
################################################################################

##  save prediction results
coef_Bridge <- matrix()       # save coefficients
pred_Bridge <- matrix()   # save prediction

for(i in 1:30){
  
  ## set seed
  set.seed(i)
  
  ## Cross-validation to obtain model parameters
  # cv.Bridge <- cv.ncpen(y.vec=y, x.mat=x, family="binomial", penalty="mbridge")
  # plot(cv.Bridge, type = "rmse", log.scale = T)
  # coef(cv.Bridge)$lambda
  
  ## fit model
  Bridge.model <- ncpen(y.vec = y, x.mat = x, family="binomial", penalty="mbridge")
  opt.lambda <- gic.ncpen(Bridge.model, pch="*", type="b")$opt.lambda
  
  print("*** optional lambda ***")
  print(opt.lambda)
  
  ## Extract coefficients
  coef <- coef(Bridge.model)
  # coef_Bridge <- as.matrix(coef[, dim(coef)[2]])
  coef <- as.matrix(coef[, dim(coef)[2]])
  tcross <- rep(i, length(coef))              
  step_Bridge <- data.frame(cbind(coef, tcross))
  coef_Bridge <- cbind(coef_Bridge, step_Bridge)   
  ## predict
  p <- predict(Bridge.model, "prob", new.x.mat = x.test)
  p <- p[,dim(p)[2]]
  kcross <- rep(i, length(p)) 
  temp_Bridge <- data.frame(cbind(y.test, p, kcross))
  pred_Bridge <- cbind(pred_Bridge,temp_Bridge)   
  
  print(paste("The：",i)) 
}


## set save pathway
# write.csv(coef_Bridge, file = "L0.5\\coef_Bridge.csv")
# write.csv(pred_Bridge, file = "L0.5\\pred_Bridge.csv")
## average value
p_Bridge <- read.table("L0.5\\pred_Bridge.csv", header=TRUE, sep = ',')
pred_Bridge <- cbind(y.test,my_pred(p_Bridge))
# write.csv(pred_Bridge, file = "L0.5\\pred_Bridge0.csv")

## ROC curve
# jpeg(file = "L0.5\\pAUC_Bridge.jpg")
plot.roc(pred_Bridge[,1], pred_Bridge[,2], print.auc=T, main="pAUC")
# dev.off()

## Performance index
predict = ifelse(pred_Bridge[,2] > 0.5, 1, 0)
predict_value = predict 
true_value = pred_Bridge[,1]
error = predict_value-true_value

accuracy = (nrow(data)-sum(abs(error)))/nrow(data) 
precision = sum(true_value & predict_value)/sum(predict_value) 
recall = sum(predict_value & true_value)/sum(true_value)       
F_measure= 2*precision*recall/(precision+recall)    
specificity = ((nrow(data)-sum(abs(error)))-sum(predict_value & true_value))/(nrow(data)-sum(true_value)) 
## AUC
pred <- prediction(pred_Bridge[,2], pred_Bridge[,1])
auc <- performance(pred,'auc')
auc <- unlist(slot(auc,'y.values'))

## Confusion matrix
# table(true_value, predict_value) %>% as.matrix(table) %>% write.csv(file = "L0.5\\table.csv")     

## Output results
i <- 4
result[i,1] <- accuracy
result[i,2] <- precision
result[i,3] <- recall
result[i,4] <- F_measure
result[i,5] <- specificity
result[i,6] <- auc
# result %>% write.csv(file = "result\\result.csv") 


################################################################################
# ncvreg: SCAD penalty  
################################################################################

X <- model.matrix(Lable ~., train.data)[,-1]   
y <- train.data$Lable 
x.test <- model.matrix(Lable ~., test.data)[,-1] 
y.test <- test.data$Lable


library(ncvreg)

##  save prediction results
coef_SCAD <- matrix()       # save coefficients
pred_SCAD <- matrix()       # save prediction

for(i in 1:30){
  
  set.seed(i)
  ## Cross-validation to obtain model parameters
  cv.SCAD <- cv.ncvreg(X, y, family ="binomial", penalty="SCAD") 
  lambda.min <- cv.SCAD$lambda.min 
  
  print("*** lambda.min ***")
  print(lambda.min)
  
  ## fit model
  SCAD.model <- ncvreg(X, y, family ="binomial", lambda = cv.SCAD$lambda.min, penalty="SCAD")
  coef <- as.matrix(coef(SCAD.model))
  tcross <- rep(i, length(coef))              
  step_SCAD <- data.frame(cbind(coef, tcross))
  coef_SCAD <- cbind(coef_SCAD, step_SCAD)   
  
  ## predict
  p <- predict(SCAD.model, x.test, type = "response")
  kcross <- rep(i, length(p)) 
  temp_SCAD <- data.frame(cbind(y.test, p, kcross))
  pred_SCAD <- cbind(pred_SCAD,temp_SCAD)   
  
  print(paste("The：",i)) 
  
}


## Store the entire y.test, p and fold
# write.csv(coef_SCAD, file = "SCAD\\coef_SCAD.csv")
# write.csv(pred_SCAD, file = "SCAD\\pred_SCAD.csv")
## average value
p_SCAD <- read.table("SCAD\\pred_SCAD.csv", header=TRUE, sep = ',')
pred_SCAD <- cbind(y.test,my_pred(p_SCAD))
# write.csv(pred_SCAD, file = "SCAD\\pred_SCAD0.csv")

## ROC curve
# jpeg(file = "SCAD\\pAUC_SCAD.jpg")
plot.roc(pred_SCAD[,1], pred_SCAD[,2], print.auc=T, main="pAUC")
# dev.off()

predict = ifelse(pred_SCAD[,2] > 0.5, 1, 0)
predict_value = predict 
true_value = pred_SCAD[,1]
error = predict_value-true_value

accuracy = (nrow(data)-sum(abs(error)))/nrow(data) 
precision = sum(true_value & predict_value)/sum(predict_value)  
recall = sum(predict_value & true_value)/sum(true_value)        
F_measure= 2*precision*recall/(precision+recall)   
specificity = ((nrow(data)-sum(abs(error)))-sum(predict_value & true_value))/(nrow(data)-sum(true_value)) 
## AUC
pred <- prediction(pred_SCAD[,2], pred_SCAD[,1])
auc <- performance(pred,'auc')
auc <- unlist(slot(auc,'y.values'))

## Confusion matrix
# table(true_value, predict_value) %>% as.matrix(table) %>% write.csv(file = "SCAD\\table.csv")     

## Output results
i <- 6
result[i,1] <- accuracy
result[i,2] <- precision
result[i,3] <- recall
result[i,4] <- F_measure
result[i,5] <- specificity
result[i,6] <- auc
result %>% write.csv(file = "result\\result.csv")


################################################################################
# ncvreg: MCP penalty 
################################################################################

library(ncvreg)
##  save prediction results
coef_MCP <- matrix()       # save coefficients
pred_MCP <- matrix()       # save prediction

for(i in 1:30){
  
  ## set seed
  set.seed(i)
  ## Cross-validation to obtain model parameters
  cv.MCP <- cv.ncvreg(X, y, family ="binomial", penalty="MCP") 
  lambda.min <- cv.MCP$lambda.min 
  
  print("*** lambda.min ***")
  print(lambda.min)
  
  ## fit model
  MCP.model <- ncvreg(X, y, family ="binomial", lambda = cv.MCP$lambda.min, penalty="MCP")
  coef <- as.matrix(coef(MCP.model))
  tcross <- rep(i, length(coef))              
  step_MCP <- data.frame(cbind(coef, tcross))
  coef_MCP <- cbind(coef_MCP, step_MCP)   
  
  ## predict
  p <- predict(MCP.model, x.test, type = "response")
  kcross <- rep(i, length(p)) 
  temp_MCP <- data.frame(cbind(y.test, p, kcross))  
  pred_MCP <- cbind(pred_MCP,temp_MCP)   
  
  print(paste("The：",i)) 
}

## set save pathway
# write.csv(coef_MCP, file = "MCP\\coef_MCP.csv")
# write.csv(pred_MCP, file = "MCP\\pred_MCP.csv")
## average value
p_MCP <- read.table("MCP\\pred_MCP.csv", header=TRUE, sep = ',')
pred_MCP <- cbind(y.test,my_pred(p_MCP))
# write.csv(pred_MCP, file = "MCP\\pred_MCP0.csv")

## ROC curve
# jpeg(file = "MCP\\pAUC_MCP.jpg")
plot.roc(pred_MCP[,1], pred_MCP[,2], print.auc=T, main="pAUC")
# dev.off()

## Performance index
predict = ifelse(pred_MCP[,2] > 0.5, 1, 0)
predict_value = predict 
true_value = pred_MCP[,1]
error = predict_value-true_value

accuracy = (nrow(data)-sum(abs(error)))/nrow(data) 
precision = sum(true_value & predict_value)/sum(predict_value)  
recall = sum(predict_value & true_value)/sum(true_value)        
F_measure= 2*precision*recall/(precision+recall)    
specificity = ((nrow(data)-sum(abs(error)))-sum(predict_value & true_value))/(nrow(data)-sum(true_value)) 
## AUC
pred <- prediction(pred_MCP[,2], pred_MCP[,1])
auc <- performance(pred,'auc')
auc <- unlist(slot(auc,'y.values'))

## Confusion matrix
# table(true_value, predict_value) %>% as.matrix(table) %>% write.csv(file = "MCP\\table.csv")     

## Output results
i <- 7
result[i,1] <- accuracy
result[i,2] <- precision
result[i,3] <- recall
result[i,4] <- F_measure
result[i,5] <- specificity
result[i,6] <- auc
# result %>% write.csv(file = "result\\result.csv")



stargazer(result)


################################################################################
# l0ara: Logistic regression + L0 penalty  
################################################################################
## Construct 'for loop' to get the AUC value predicted by ten-fold cross-validation. 
## And record the group with the largest value as the optimal training set and test set division.

# install.packages("l0ara")
library("l0ara")

coef_L0 <- matrix()   # save coefficients
pred_L0 <- matrix()   # save prediction

for(i in 1:30){
  # i <- 1
  ## set seed
  set.seed(i)
  ## Cross-validation to obtain model parameters
  # lam <- seq(1,0.05,-0.05)
  lam <- seq(1,0.5,-0.5)
  cv.l0 <- cv.l0ara(x, y, family="logit", lam, measure = "mse")
  lambda.min <- cv.l0$lam.min
  
  print("*** lambda.min ***")
  print(lambda.min)
  
  ## fit model
  l0.model <- l0ara(x, y, family = "logit", lambda.min)
  coef_l0 <- coef(l0.model)
  coef_l0 = as.matrix(coef_l0)
  
  new_coef_l0 <- matrix()
  for (j in 2:nrow(coef_l0)){
    new_coef_l0[j+1] <- coef_l0[j]
  }
  for (k in 1:2){
    new_coef_l0[k] <- coef_l0[k]
  }
  
  coef <- as.matrix(new_coef_l0)
  tcross <- rep(i, length(coef))              
  step_L0 <- data.frame(cbind(coef, tcross))
  coef_L0 <- cbind(coef_L0, step_L0)   
  
  ## predict
  p <- predict(l0.model, newx = x.test, type = "response")
  kcross <- rep(i, length(p)) 
  temp_L0 <- data.frame(cbind(y.test, p, kcross)) 
  pred_L0 <- cbind(pred_L0,temp_L0)   
  
  print(paste("The：",i)) 
}


## Store the entire y.test, p and fold
# write.csv(coef_L0, file = "L0\\coef_L0.csv")
# write.csv(pred_L0, file = "L0\\pred_L0.csv")

## average value
p_L0 <- read.table("L0\\pred_L0.csv", header=TRUE, sep = ',')

my_pred <- function(x){
  p_ridge1 <- x[,-2]
  p_ridge2 <- matrix(data=0, nrow = dim(p_ridge1)[1], ncol = 1, byrow = FALSE, dimnames=list(c(as.character(x[,1])),c("prob")))
  for (j in 1:10){
    p_ridge2 <- p_ridge2[] + p_ridge1[,3*j]
  }
  p_ridge3 <- p_ridge2/10
  return(p_ridge3)
}
pred_L0 <- cbind(y.test,my_pred(p_L0))
# pred_L0 <- cbind(y,my_pred(p_L0))
# write.csv(pred_L0, file = "L0\\pred_L0.csv")

## ROC curve
# jpeg(file = "L0\\pAUC_L01.jpg")
plot.roc(pred_L0[,1], pred_L0[,2], print.auc=T, main="pAUC")
# dev.off()

## Performance index
predict = ifelse(pred_L0[,2] > 0.5, 1, 0)
predict_value = predict
true_value = pred_L0[,1]
error = predict_value-true_value

accuracy = (nrow(data)-sum(abs(error)))/nrow(data)
precision = sum(true_value & predict_value)/sum(predict_value)  
recall = sum(predict_value & true_value)/sum(true_value)        
F_measure = 2*precision*recall/(precision+recall)    
specificity = ((nrow(data)-sum(abs(error)))-sum(predict_value & true_value))/(nrow(data)-sum(true_value))
## AUC
pred <- prediction(pred_L0[,2], pred_L0[,1])
auc <- performance(pred,'auc')
auc <- unlist(slot(auc,'y.values'))

## Confusion matrix
# table(true_value, predict_value) %>% as.matrix(table) %>% write.csv(file = "L0\\table.csv")     

## Output results
i <- 5
result[i,1] <- accuracy
result[i,2] <- precision
result[i,3] <- recall
result[i,4] <- F_measure
result[i,5] <- specificity
result[i,6] <- auc
# result %>% write.csv(file = "result\\result.csv")


################################################################################
# gene integration 
################################################################################

## clear
rm(list = ls())


## load data (i.e., the coefficients of above seven methods)
setwd(".\\GSE59491_15\\Data37")
Coef_ridge <- read.table("ridge\\coef_ridge.csv", header=TRUE, sep = ',')
Coef_lasso <- read.table("lasso\\coef_lasso.csv", header=TRUE, sep = ',')
Coef_Elastic <- read.table("elastic_net\\coef_elastic.csv", header=TRUE, sep = ',')
Coef_Bridge <- read.table("l0.5\\coef_Bridge.csv", header=TRUE, sep = ',')
Coef_l0 <- read.table("L0\\coef_L0name.csv", header=TRUE, sep = ',')
Coef_SCAD <- read.table("SCAD\\coef_SCAD.csv", header=TRUE, sep = ',')
Coef_MCP <- read.table("MCP\\coef_MCP.csv", header=TRUE, sep = ',')


## extract coefficients of 30 runs
## define 'my_cbind' function to extract 30-times coefficients -----------------
my_cbind <- function(x){
  x1 <- matrix()
  x1 <- as.character(x[,1])
  for (i in 0:29){
    x1 <- cbind(x1, x[, 3+2*i])
  }
  return(x1)
}


## compute ---------------------------------------------------------------------
coef_ridge1 <- my_cbind(Coef_ridge)
coef_lasso1 <- my_cbind(Coef_lasso)
coef_Elastic1 <- my_cbind(Coef_Elastic) 
coef_Bridge1 <- my_cbind(Coef_Bridge)
coef_l01 <- my_cbind(Coef_l0)  
coef_SCAD1 <- my_cbind(Coef_SCAD) 
coef_MCP1 <- my_cbind(Coef_MCP)  


## extract genes ---------------------------------------------------------------
## define 'my_union_coef', getunion of 30-times non-zeros coefficients ---------
my_union_coef <- function(x){
  x1 <- matrix(data=NA)
  j <- 1
  for (i in 2:dim(x)[1]){
    # i <- 2
    if (x[i,2] !=0 | x[i,3] !=0 | x[i,4] !=0 | x[i,5] !=0
        | x[i,6] !=0 | x[i,7] !=0 | x[i,8] !=0 | x[i,9] !=0
        | x[i,10] !=0 | x[i,11] !=0 | x[i,12] !=0 | x[i,13] !=0
        | x[i,14] !=0 | x[i,15] !=0 | x[i,16] !=0 | x[i,17] !=0
        | x[i,18] !=0 | x[i,19] !=0 | x[i,20] !=0 | x[i,21] !=0
        | x[i,22] !=0 | x[i,23] !=0 | x[i,24] !=0 | x[i,25] !=0
        | x[i,26] !=0 | x[i,27] !=0 | x[i,28] !=0 | x[i,29] !=0
        | x[i,30] !=0 | x[i,31] !=0) {
      x1[j] <- as.character(x[i,1])
      j <- j+1
    }else{
      print("Wrong!")
    }
  }
  return(as.matrix(x1))
}  


## Selected gene ---------------------------------------------------------------
coef_ridge <- my_union_coef(coef_ridge1)
coef_lasso <- my_union_coef(coef_lasso1)
coef_Elastic <- my_union_coef(coef_Elastic1)
coef_Bridge <- my_union_coef(coef_Bridge1)
coef_l0 <- my_union_coef(coef_l01)
coef_SCAD <- my_union_coef(coef_SCAD1)
coef_MCP <- my_union_coef(coef_MCP1)


## compute 2-2 Overlap ---------------------------------------------------------
gene_rl <- intersect(coef_ridge, coef_lasso)
gene_re <- intersect(coef_ridge, coef_Elastic)
gene_r1 <- intersect(coef_ridge, coef_l0)
gene_rb <- intersect(coef_ridge, coef_Bridge)
gene_rs <- intersect(coef_ridge, coef_SCAD)
gene_rm <- intersect(coef_ridge, coef_MCP)

gene_le <- intersect(coef_lasso, coef_Elastic)
gene_l1 <- intersect(coef_lasso, coef_l0)
gene_lb <- intersect(coef_lasso, coef_Bridge)
gene_ls <- intersect(coef_lasso, coef_SCAD)
gene_lm <- intersect(coef_lasso, coef_MCP)

gene_e1 <- intersect(coef_Elastic, coef_l0)
gene_eb <- intersect(coef_Elastic, coef_Bridge)
gene_es <- intersect(coef_Elastic, coef_SCAD)
gene_em <- intersect(coef_Elastic, coef_MCP)

gene_1b <- intersect(coef_l0, coef_Bridge)
gene_1s <- intersect(coef_l0, coef_SCAD)
gene_1m <- intersect(coef_l0, coef_MCP)

gene_bs <- intersect(coef_Bridge, coef_SCAD)
gene_bm <- intersect(coef_Bridge, coef_MCP)

gene_sm <- intersect(coef_SCAD, coef_MCP)

# View(gene_rl)  
# View(gene_re)  
# View(gene_r1)  
# View(gene_rb)  
# View(gene_rs)
# View(gene_rm)  
# 
# View(gene_le)  
# View(gene_l1)  
# View(gene_lb)  
# View(gene_ls)  
# View(gene_lm)  
# 
# View(gene_e1)
# View(gene_eb) 
# View(gene_es) 
# View(gene_em) 
# 
# View(gene_1b)
# View(gene_1s)
# View(gene_1m)
# 
# View(gene_bs)
# View(gene_bm)
# 
# View(gene_sm)



## compute number of features  -------------------------------------------------
sum(coef_ridge != 0)
sum(coef_lasso != 0)
sum(coef_Elastic != 0)
sum(coef_l0 != 0)
sum(coef_Bridge != 0)
sum(coef_SCAD != 0)
sum(coef_MCP != 0)


## save data  ------------------------------------------------------------------
setwd(".\\GSE59491_15\\Data37\\result")
# write.csv(coef_ridge, "coef_ridge1.csv")
# write.csv(coef_lasso, "coef_lasso1.csv")
# write.csv(coef_Elastic, "coef_Elastic1.csv")
# write.csv(coef_l0, "coef_l01.csv")
# write.csv(coef_Bridge, "coef_Bridge1.csv")
# write.csv(coef_SCAD, "coef_SCAD1.csv")
# write.csv(coef_MCP, "coef_MCP1.csv")

## save biomarlers
gene <- intersect(intersect(intersect(coef_Bridge, coef_lasso), coef_Elastic), coef_SCAD)
# write.csv(gene, "gene_overlap_scale20.csv")
stargazer(gene)  # Use stargazer to generate descriptive tables in LaTeX.


################################################################################
# plot 'ROC curve' of seven methods
################################################################################

## clear
rm(list = ls())


## load data
setwd(".\\GSE59491_15\\Data37")
coef_ridge <- read.table("ridge\\pred_ridge0.csv", header=TRUE, sep = ',')
coef_lasso <- read.table("lasso\\pred_lasso0.csv", header=TRUE, sep = ',')
coef_Elastic <- read.table("elastic_net\\pred_elastic0.csv", header=TRUE, sep = ',')
coef_Bridge <- read.table("l0.5\\pred_Bridge0.csv", header=TRUE, sep = ',')
coef_l0 <- read.table("l0\\pred_L00.csv", header=TRUE, sep = ',')
coef_SCAD <- read.table("SCAD\\pred_SCAD0.csv", header=TRUE, sep = ',')
coef_MCP <- read.table("MCP\\pred_MCP0.csv", header=TRUE, sep = ',')


library(pROC)
# pdf(file = "result\\pAUCall7.pdf")
roc_Elastic <- plot.roc( coef_Elastic[,2], coef_Elastic[,3], col="Salmon" )
roc_lasso <- lines.roc( coef_lasso[,2], coef_lasso[,3], col="Aquamarine" )
roc_SCAD <- lines.roc( coef_SCAD[,2], coef_SCAD[,3], col="Magenta" )# VioletRed
roc_Bridge <- lines.roc( coef_Bridge[,2], coef_Bridge[,3], col="Green" )
roc_l0 <- lines.roc( coef_l0[,2], coef_l0[,3], col="Tan" )
roc_MCP <- lines.roc( coef_MCP[,2], coef_MCP[,3], col="Cyan" )
roc_ridge <- lines.roc( coef_ridge[,2], coef_ridge[,3], col="Orchid" )  # Purple
legend("bottomright", legend=c("Elastic Net", "Lasso", "SCAD", "L1/2", "L0", "MCP", "Ridge"), 
       col=c("Salmon", "Aquamarine", "Magenta", "Green", "Tan", "Cyan", "Orchid"), lwd=2)
# dev.off()


## load data
bar <- read.table("result\\bar.csv", header=TRUE, sep = ',')
colnames(bar) <- c("AUC","Penalty")
bar$AUC <- round(bar$AUC,3)         # Keep three decimal places
library(ggplot2)
# pdf(file = "result\\bar_fig.pdf")
ggplot(data = bar, aes(x = Penalty, y = AUC, fill = Penalty))+
  geom_bar(stat = 'identity')+
  geom_text(aes(label = AUC))
# dev.off()


################################################################################
# the end
################################################################################
