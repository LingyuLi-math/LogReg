## Used to perform T test on the data after labeling series_marix
## Determine the p value of the T test and select genes with p/pdf <0.05


## clear
rm(list = ls())

## package
library(BiocGenerics)
library(parallel)
library(Biobase)
library(dplyr)        
library(tidyr)
library(tidyverse)    
library(fdrtool)      


## load function
source("C:\\Users\\LiLingyu\\Desktop\\LogReg\\R\\Ttest.R")   # change the pathway using your storaged pathway of "Ttest.R"


## load data
setwd(".\\GSE59491_15")
x1 = read.table("Data/GSE59491_outcome.txt", header = T, sep='\t', fill=TRUE, strip.white = T, check.names = F)
GSE59491 <- t(x1)
dim(x1)    # 24479   326

##########################################################################################
# GSE59491 Differential gene expression (t-test) 
##########################################################################################
## genes that  p < 0.05 
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

## Scale
data59491 <- my_scale(data.frame(t(GSE59491_T)))
dim(data59491)    # 326 360
View(data59491[,1:10])
# write.table(data59491,"GSE59491_scale_DE.txt",quote=F,sep="\t") 

##########################################################################################
# End the analysis for GSE59491
##########################################################################################
