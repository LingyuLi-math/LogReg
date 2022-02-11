
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
setwd("D:\\E\\��ʿ\\R_����\\GSE59491_15")
expSet <- read.table("Data/GSE59491_series_matrix.txt", header=T, sep='\t', fill=TRUE, strip.white = T)
expSet$ID_REF <- as.character(expSet$ID_REF)  # ��ID_REF��ȫ��ת���ɷ���,Ϊ��ͬanno2�ϲ�


## load annotation
anno <- read.table("Data/GPL18964.txt",header=T,sep='\t',fill=TRUE,strip.white = T,quote = "")  # ̽��ͻ���ID�Ķ�Ӧ�ļ�
anno2 <- anno[,c('ID','ENTREZ_GENE_ID')]     # ��ȡ�����б�ǩ
colnames(anno2) <- c('ID_REF','EntrezID')    # �������еı�ǩ�滻, 'GeneSymbol''Gene.Symbol'
anno2$ID_REF <- as.character(anno2$ID_REF)   # ��ID_REF��ȫ��ת���ɷ���,Ϊ��ͬexpSet�ϲ�


## ���������������оƬע���ļ���̽�������ж�Ӧ
expset2 <- expSet %>%                      # ��>������dplyr���Ĺܵ�������
  inner_join(anno2,by='ID_REF') %>%        # �����ǽ�ǰһ���Ľ��ֱ�Ӵ��θ���һ���ĺ�����ʡ���м�ĸ�ֵ����
  select(ID_REF,EntrezID, everything())    # %>%# ��������
# View(expset2[,1:3])                      # expset2 �� expset ��ȣ���2�ж��˻����


## ����оƬע���ļ���������һ��̽���Ӧ�������Ĳ�ֿ�
expset3 <- expset2
a <- tibble(expset3[,1:2])                 # �ѵ� 1 �͵� 2 ����ȡ���������� a ��  
test1 <- apply(a,1, function(x){
  str_split(x[2], '///', simplify=T)       # �� a �ĵ�2����ȡ������test1
} )


test2 <- apply(a, 1, function(x){          # ��̽��źͻ���ţ�����---������
  paste(x[1],str_split(x[2], '///', simplify=T), sep = "---")
})


unlist(test2)                              # �� list ���ݱ���ַ�����������������������ʽ
x <- tibble(unlist(test2))                 # tibble��ȡ����ͳdata.frame����ȡ���Զ�����������unlist(test2)
colnames(x) <- "lala"                      # �ı� x ���������� unlist(test2) ����Ϊ lala
x2 <- separate(x,lala,c("id","entrezID"),sep = '---')     # ʶ�� lala �е� ---�������ݷ��룬�������в����±�ǩ
x3 <- merge(x2,expset3,by.x = "id", by.y="ID_REF", all=FALSE)  #  �������ļ���˳��ϲ�Ϊһ����
x4<-x3[,-c(1,3)]                           # �� ��1 �͵�3 ����ɾ��, ʣ�µ����ݻ����ַ��͵ģ����š� "
zz <- as.matrix(apply(as.matrix(x4[,1]),1,function(x) as.numeric(x)))
XX <- x4[,-1]
colnames(XX)[1:3]
XX1 <- cbind(zz,XX)
colnames(XX1) <- c("entrezID",colnames(XX))


## �û���id�������õ�оƬע���ļ����л������ĸ���
homo<-read.table("Data/homo.txt",header=T,sep='\t')
x5 <- merge(homo, XX1, by.x="GeneID", by.y = "entrezID", all=FALSE) 
# �ϲ��� x5 ��25088����Ϊ24478������������������ƥ��


## ̽����ƥ���������ȡ�����̽���Ӧһ����������ݼ���IQR������IQR����̽������
expset4 <- x5 %>% 
  dplyr::select(-GeneID) %>%              # ȥ��������Ϣ
  mutate(rowIQR =apply(.[,-1],1,IQR)) %>% # ����ÿ�е�IQR
  arrange(desc(rowIQR)) %>%               # �ѱ�������ƽ��ֵ���Ӵ�С����
  distinct(Symbol,.keep_all = T) %>%      # symbol���µ�һ��
  dplyr::select(-rowIQR)   %>%            # ����ѡ��ȥ��rowIQR��һ��
  tibble::column_to_rownames(colnames(.)[1]) # �ѵ�һ�б��������ɾ��
View(expset4[1:10,1:10])
dim(expset4)    # 24478   326
  

# ��ǩ ----------------------------------------------------------------------
lable2 = read.csv("Data/GSE59491_all.csv", header = T, sep=',')
dim(lable2)    # 326   2
data = rbind(as.matrix(t(lable2[,2])), as.matrix(expset4))
rownames(data) <- c('Lable', rownames(expset4))
View(data[1:10,1:10])
dim(data)    #  24479   326
# write.table(data,"Data/GSE59491_outcome.txt",quote=F,sep="\t")

