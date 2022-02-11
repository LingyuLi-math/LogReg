library(tidyverse)
library(plyr)
library(stargazer) #转换latex 


# 计算超几何分布检验 ---------------------------------------------------------------

## 真实
N <- 359  # 总数
M <- 359  # 方法一
n <- 67   # 方法二
k <- 67   # overlap

# my_phyper 函数 ---------------------------------------------------------------

my_phyper <- function(N,M,n,k){
  q <- k
  m <- M
  p <- N-M
  l <- n
  p_value <- phyper(q-1, m, p, l, lower.tail=F)
  return(p_value)
}


# 结果 ----------------------------------------------------------------------

p_rl <- my_phyper(359,359,67,67)
p_re <- my_phyper(359,359,78,78)
p_r1 <- my_phyper(359,359,31,31)
p_rb <- my_phyper(359,359,25,25)
p_rs <- my_phyper(359,359,49,49)
p_rm <- my_phyper(359,359,27,27)

p_le <- my_phyper(359,67,78,67)
p_l1 <- my_phyper(359,67,31,28)
p_lb <- my_phyper(359,67,25,21)
p_ls <- my_phyper(359,67,49,48)
p_lm <- my_phyper(359,67,27,27)

p_e1 <- my_phyper(359,78,31,29)
p_eb <- my_phyper(359,78,25,21)
p_es <- my_phyper(359,78,49,48)
p_em <- my_phyper(359,78,27,27)

p_1b <- my_phyper(359,31,25,13)
p_1s <- my_phyper(359,31,49,19)
p_1m <- my_phyper(359,31,27,24)

p_bs <- my_phyper(359,25,49,20)
p_bm <- my_phyper(359,25,27,13)

p_sm <- my_phyper(359,49,27,27)


# 输出结果 --------------------------------------------------------------------


result <- matrix(0, 6, 6)
colnames(result) <- c( "lasso", "elatic net", "L1/2", "L0", "SCAD", "MCP")
rownames(result) <- c("ridge", "lasso", "elatic net", "L1/2", "L0", "SCAD")

i <- 1
result[i,1] <- signif(p_rl, digits = 3)
result[i,2] <- signif(p_re, digits = 3)
result[i,3] <- signif(p_r1, digits = 3)
result[i,4] <- signif(p_rb, digits = 3)
result[i,5] <- signif(p_rs, digits = 3)
result[i,6] <- signif(p_rm, digits = 3)

result[i+1,2] <- signif(p_le, digits = 3)
result[i+1,3] <- signif(p_l1, digits = 3)
result[i+1,4] <- signif(p_lb, digits = 3)
result[i+1,5] <- signif(p_ls, digits = 3)
result[i+1,6] <- signif(p_lm, digits = 3)

result[i+2,3] <- signif(p_e1, digits = 3) # format(p_e1, scientific = T)
result[i+2,4] <- signif(p_eb, digits = 3) 
result[i+2,5] <- signif(p_es, digits = 3)
result[i+2,6] <- signif(p_em, digits = 3)

result[i+3,4] <- signif(p_1b, digits = 3)
result[i+3,5] <- signif(p_1s, digits = 3)
result[i+3,6] <- signif(p_1m, digits = 3)

result[i+4,5] <- signif(p_bs, digits = 3)
result[i+4,6] <- signif(p_bm, digits = 3)

result[i+5,6] <- signif(p_sm, digits = 3)

setwd("D:\\E\\博士\\R_程序\\GSE59491_15")
# result %>% write.csv(file = "Data37\\result\\phyper.csv") 

# 表格 ----------------------------------------------------------------------
stargazer(result) 

# 百分比填色 -------------------------------------------------------------------
result <- as.matrix(result)

p <- matrix(data = NA, nrow = 21, ncol = 1,byrow = F, dimnames = NULL)

for (j in 1:6){
  for (i in j:6){
    p[6*(j-1)+i] <- result[j,i]
  }
}

View(p)


# 删除缺失值 -------------------------------------------------------------------
p <- data.frame(p)
p_new <- na.omit(p)
View(p_new)



# normalized --------------------------------------------------------------
p_matrix <- as.matrix(p_new)[7:21]
View(p_matrix)
p_sort <- sort(p_matrix)
View(p_sort)
p_scale <- scale(p_matrix)
View(p_scale)


# sort --------------------------------------------------------------------
sort1 <- matrix(data = NA, nrow = 15, ncol = 1,byrow = F, dimnames = NULL)

for (k in 1:15){
  sort1[k] <- p_sort[k]/p_sort[15]
}
View(sort1)
