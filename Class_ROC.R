# 2020.6.29

## clear
rm(list = ls())


## package
library(pROC)


## load data
setwd("D:\\E\\博士\\R_程序\\GSE59491_15")
Data  = read.table("Data\\GSE73685_17.txt", header = T, check.names = FALSE)
Data2 = read.table("Data\\GSE59491_17.txt", header = T, check.names = FALSE)


# split data -------------------------------------------------------------------
x.train <- data.frame(t(Data2)[,-1])
y.train <- t(Data2)[,1]
x.test <- data.frame(t(Data)[,-1])
y.test <- t(Data)[,1]


# fit ---------------------------------------------------------------------
glm.fit <- glm(y.train~., data = x.train, family = binomial)
# glm.fit <- glm(y.train~., data = x.train, family = binomial, control = list(maxit = 100))
summary(glm.fit)


# train --------------------------------------------------------------------
p_train <- predict(glm.fit, x.train, type = "response")
p_train = as.matrix(p_train)
A_train <- data.frame(p_train, y.train)
names(A_train)<- c("p", "outcome")
plot.roc(A_train$outcome, A_train$p, print.auc=T, main="pAUC")


# test --------------------------------------------------------------------
p_test <- predict(glm.fit, x.test, type = "response")
pred_glm <- cbind(p_test, y.test)
colnames(pred_glm) <- c('y.test', 'Lable')
p <- pred_glm[,1]
p_glm <- cbind(log(p/(1-p)), pred_glm[,2])
colnames(p_glm) <- c('y.test', 'Lable')
# write.table(p_glm,"Data37\\DE\\pre59491_46510_7cvlog_zf_new.txt",quote=F,sep="\t")  # 2020.6.30


library(pROC)
p_test = as.matrix(p_test)
A_test <- data.frame(p_test, y.test)
names(A_test)<- c("p", "outcome")
# pdf(file = "Data37\\DE\\pAUC_test_glm_cv20_zf_new.pdf",width = 4,height = 4)  # 2020.6.29后
plot.roc(A_test$outcome, A_test$p, print.auc=T)
# plot.roc(A_test$outcome, A_test$p)
# legend("bottomright", legend=c("Acc = 0.941", "Pre = 1.000 ", "Sn = 0.800", "F-measure = 0.889", "Sp = 1.000", "AUC = 0.933"))
# dev.off()


## performance 
predict = ifelse(pred_glm[,1] > 0.65, 1, 0)  # 2020.6.29 later
predict_value = predict
true_value = pred_glm[,2]
error = predict_value-true_value

data <- t(Data)
# 计算模型准确性（accuracy）、精确度（Precision），召回率（Recall）-- 敏感性（sensitivity）-- 真阳性率（TPR）和 F测度（F-measure）和混淆矩阵
accuracy = (nrow(data)-sum(abs(error)))/nrow(data)
precision = sum(true_value & predict_value)/sum(predict_value)  
recall = sum(predict_value & true_value)/sum(true_value)        
F_measure = 2*precision*recall/(precision+recall)    
specificity = ((nrow(data)-sum(abs(error)))-sum(predict_value & true_value))/(nrow(data)-sum(true_value))

accuracy
precision
recall
F_measure
specificity

## 混淆矩阵，显示结果依次为TP、FN、FP、TN
table(true_value, predict_value) 

