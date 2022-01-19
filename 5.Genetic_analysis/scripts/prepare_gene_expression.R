setwd("/public10/home/sci0011/projects/tomato2/12_pheno/03_exp_meta")
dat =read.table("/public10/home/sci0010/project/00_new_tomato/02.kallisto/sample_reverse_trans",head=T)
ID = read.table("332.txt",head=F)
dat.back = dat
idx = match(ID$V1,dat[,1])
dat = dat[idx,]
na.num = apply(dat[,-1],2,function(x){
  sum(x>0.5)
})
sum(na.num > 99) ## 26902; 24199
idx = na.num > 99
dat2 = dat[,-1]
dat3 = dat2[,idx]
ID = as.character(dat[,1])
exp = data.frame(taxa=ID,dat3)
write.table(exp,"pheno_exp_19353_ori.txt",col.names = T,row.names = F,quote=F,sep="\t")
exp2 = exp
exp3 = exp2[,-1]
library(edgeR)
pc = prcomp(t(exp3))
pcs = pc$x
dim(pcs)
plot(pcs[,1:2])
which(colnames(exp3) == "Solyc09T002497.1")
which(colnames(exp3) == "Solyc04T002594.1")
which(colnames(exp3) == "Solyc01T002129.1")

quantile_normalisation <- function(df){
  df_rank <- apply(df,2,rank,ties.method="min")
  df_sorted <- data.frame(apply(df, 2, sort))
  df_mean <- apply(df_sorted, 1, mean)
  
  index_to_mean <- function(my_index, my_mean){
    return(my_mean[my_index])
  }
  
  df_final <- apply(df_rank, 2, index_to_mean, my_mean=df_mean)
  rownames(df_final) <- rownames(df)
  return(df_final)
}
exp4 = quantile_normalisation(exp3)
data = data.frame(FID=exp$taxa,IID=exp$taxa,father=rep(0,332),mother=rep(0,332),sex=rep(0,332),exp4)

## local
setwd("/Users/yaozhou/Documents/MyNutstore/yaolab/Projects/GraphGWAS/results/population/01_prepare_sample")
exp = read.table("pheno_exp_18861_ori.txt",head=T)
exp2 = exp
exp3 = exp2[,-1]
library(edgeR)
pc = prcomp(t(exp3))
pcs = pc$x
dim(pcs)
plot(pcs[,1:2])
which(colnames(exp3) == "Solyc09T002497.1")
which(colnames(exp3) == "Solyc04T002594.1")

quantile_normalisation <- function(df){
  df_rank <- apply(df,2,rank,ties.method="min")
  df_sorted <- data.frame(apply(df, 2, sort))
  df_mean <- apply(df_sorted, 1, mean)
  
  index_to_mean <- function(my_index, my_mean){
    return(my_mean[my_index])
  }
  
  df_final <- apply(df_rank, 2, index_to_mean, my_mean=df_mean)
  rownames(df_final) <- rownames(df)
  return(df_final)
}
exp4 = quantile_normalisation(exp3)
data = data.frame(FID=exp$taxa,IID=exp$taxa,father=rep(0,332),mother=rep(0,332),sex=rep(0,332),exp4)
write.table(data,"pheno_exp_19353.norm.txt",col.names = T,row.names = F,quote = F,sep="\t")
