## 用于对 series_marix 添加标签后的数据，进行 T 检验
## 确定 T 检验的 p 值，并将 p/pdf <0.05 的gene挑出来


## clear
rm(list = ls())

## package
library(BiocGenerics)
library(parallel)
library(Biobase)
library(dplyr)       # ％>％ 管道函数的调用，传参
library(tidyr)
library(tidyverse)   # tibble 的调用
library(fdrtool)     # fdr校正


## load function
source("C:\\Users\\LiLingyu\\Desktop\\LogReg\\R\\Ttest.R")


## load data
setwd("D:\\E\\博士\\R_程序\\GSE59491_15")
x1 = read.table("Data\GSE59491_outcome.txt", header = T, sep='\t', fill=TRUE, strip.white = T, check.names = F)
GSE59491 <- t(x1)
dim(x1)    # 24479   326


############################  GSE59491 差异基因表达（t 检验） #############################
# p < 0.05 的   2,809 
test1 <- as.matrix(my_test_1(GSE59491))
test0 <- as.matrix(my_test_0(GSE59491))
p59491 <- my_p(GSE59491)                    
# write.csv(p59491,"p59491.csv")

p59491_BH_fdr <- my_BH_fdr(p59491)
# write.csv(p59491_BH_fdr,"p59491_BH_fdr.csv")

# p59491_fdr <- my_fdr(p59491,p59491)
# View(p59491_fdr)
# write.csv(p59491_fdr,"p59491_fdr.csv")

GSE59491_T <- my_p_data(p59491_BH_fdr, GSE59491)
# GSE59491_T <- my_p_data(p59491, GSE59491)   # GSE73685
dim(GSE59491_T)    # 326 360
# write.table(GSE59491_T,"GSE59491_DE.txt",quote=F,sep="\t") 

## 标准化,scale的数据时data.fram
data59491 <- my_scale(data.frame(t(GSE59491_T)))
dim(data59491)    # 326 360
View(data59491[,1:10])
# write.table(data59491,"Data/GSE59491_scale_DE.txt",quote=F,sep="\t") 

