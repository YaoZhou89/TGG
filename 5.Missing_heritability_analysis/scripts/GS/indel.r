#!/bin/bash Rscript
args <- commandArgs()
rep = args[6]
fold = args[7]
k="indels"
phe = read.table("/public10/home/sci0011/projects/tomato2/37_GS/01_phenotype/meta.txt",head=T)
ksnp = read.table("/public10/home/sci0011/projects/tomato2/37_GS/02_snps_only/K.txt",head=F)
kindel = read.table("/public10/home/sci0011/projects/tomato2/37_GS/03_indels_only/K.txt",head=F)
ksv = read.table("/public10/home/sci0011/projects/tomato2/37_GS/04_svs_only/K.txt",head=F)
# K = list(ksnp= as.matrix(ksnp),kindel = as.matrix(kindel),ksv = as.matrix(ksv))
cov =read.table("/public10/home/sci0011/projects/tomato2/17_cov/01_snps/pheno.cov",head=T)
Y = read.table(paste("/public10/home/sci0011/projects/tomato2/37_GS/01_phenotype/rep",rep,"_fold",fold,".txt",sep=""),head=T)

res = matrix(NA,5,971)
res = as.data.frame(res)
res[1,1] = k
res[1,2] = c("h2")

res[2:4,1] = k
res[2:4,2] = c("h2","fit","acc")

colnames(res) = c("type",colnames(Y)[5:974])

for (p in 5:974){
  data = data.frame(y=Y[,p],gid=Y$IID)
  idx = is.na(Y[,p])
  vc <- hiblup.vc(y = Y[,p], K = list(kindel), X = cov[,3:6],method = "HE",
                  nAIiter = 1000, blup.solution = T, verbose = TRUE)
  
  fixed = as.matrix(cov[,3:6]) %*% as.vector(vc$beta)[2:5]
  fixed = NULL
  random = vc$u
  blup = apply(cbind(fixed,random),1,sum)
  res[1,p-2] = vc$vc[1]/sum(vc$vc)
  res[2,p-2] = vc$vc[2]/sum(vc$vc)
  res[3,p-2] = sum(vc$vc[1:2])/sum(vc$vc)
  res[4,p-2] = cor(phe[!idx,5],blup[!idx],use="pairwise.complete.obs")
  res[5,p-2] = cor(phe[idx,5],blup[idx],use="pairwise.complete.obs")
}
write.table(res,paste("rep",rep,"_fold",fold,".",k,".txt",sep=""),col.names = T,row.names = F,quote=F,sep="\t")


