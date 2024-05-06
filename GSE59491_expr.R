
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


## load data
setwd(".\\GSE59491_15")
expSet <- read.table("Data/GSE59491_series_matrix.txt", header=T, sep='\t', fill=TRUE, strip.white = T)
expSet$ID_REF <- as.character(expSet$ID_REF)  # Convert all ID_REF columns into symbols, and merge with anno2


## load annotation
anno <- read.table("Data/GPL18964.txt",header=T,sep='\t',fill=TRUE,strip.white = T,quote = "")  # 探针和基因ID的对应文件
anno2 <- anno[,c('ID','ENTREZ_GENE_ID')]     # Extract these two columns of labels
colnames(anno2) <- c('ID_REF','EntrezID')    # Replace these two columns's labels, with 'GeneSymbol''Gene.Symbol'
anno2$ID_REF <- as.character(anno2$ID_REF)   # ConvertID_REF columns to symbols, in order to merge with expSet


## Correspond gene expression data to probe names (in chip annotation files)
expset2 <- expSet %>%                      # ％>％, pipe function from dplyr package,
  inner_join(anno2,by='ID_REF') %>%        # Directly pass result of previous step to the function of next step, omitting intermediate step.
  select(ID_REF,EntrezID, everything())    # Rearrange
# View(expset2[,1:3])                      # Compared with expset, expset2 has an extra gene number in column 2


## Organize the chip annotation files and separate the multiple genes corresponding to one probe.
expset3 <- expset2
a <- tibble(expset3[,1:2])                 # Extract columns 1 and 2 and put them in ‘a’
test1 <- apply(a,1, function(x){
  str_split(x[2], '///', simplify=T)       # Extract columns 2 of ‘a’ and give it to test1
} )


test2 <- apply(a, 1, function(x){          # Link the probe number and gene number with ‘---’
  paste(x[1],str_split(x[2], '///', simplify=T), sep = "---")
})


unlist(test2)                              # Convert list data into a string vector or numeric vector
x <- tibble(unlist(test2))                 # tibble，replaces traditional data.frame, reads and automatically adds column names: unlist(test2)
colnames(x) <- "lala"                      # Change the column name of x: define unlist(test2) as lala
x2 <- separate(x,lala,c("id","entrezID"),sep = '---')     # Identify '---' in 'lala', separate the data into separate columns and attach new labels
x3 <- merge(x2,expset3,by.x = "id", by.y="ID_REF", all=FALSE)  # Merge two files into one sequentially
x4<-x3[,-c(1,3)]                           # Delete the 1st and 3rd columns, and the remaining data is still in character type, with ""
zz <- as.matrix(apply(as.matrix(x4[,1]),1,function(x) as.numeric(x)))
XX <- x4[,-1]
colnames(XX)[1:3]
XX1 <- cbind(zz,XX)
colnames(XX1) <- c("entrezID",colnames(XX))


## Use the gene ID to update the gene name of the organized chip annotation file
homo<-read.table("Data/homo.txt",header=T,sep='\t')
x5 <- merge(homo, XX1, by.x="GeneID", by.y = "entrezID", all=FALSE) 
## x5 is reduced from 25088 to 24478, and the gene number is matched with the gene name.


## Probe name matches gene name: taken out the data (that multiple probes corresponding to one gene)
## Calculate the IQR, and the probe data with the largest IQR is retained.
expset4 <- x5 %>% 
  dplyr::select(-GeneID) %>%              # Remove redundant information
  mutate(rowIQR =apply(.[,-1],1,IQR)) %>% # Calculate the IQR of each row
  arrange(desc(rowIQR)) %>%               # Sort the average expression values from largest to smallest
  distinct(Symbol,.keep_all = T) %>%      # symbol, leaves the first one
  dplyr::select(-rowIQR)   %>%            # Reverse selection to remove rowIQR column
  tibble::column_to_rownames(colnames(.)[1]) # Reverse selection to remove rowIQR column...
View(expset4[1:10,1:10])
dim(expset4)    # 24478   326


## Label  ----------------------------------------------------------------------
lable2 = read.csv("Data/GSE59491_all.csv", header = T, sep=',')
dim(lable2)    # 326   2
data = rbind(as.matrix(t(lable2[,2])), as.matrix(expset4))
rownames(data) <- c('Lable', rownames(expset4))
View(data[1:10,1:10])
dim(data)    #  24479   326
# write.table(data,"Data/GSE59491_outcome.txt",quote=F,sep="\t")


