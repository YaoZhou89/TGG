library(rrBLUP)
args <- commandArgs()
rep = args[6]
fold = args[7]
phe = read.table("/public10/home/sci0011/projects/tomato2/37_GS/01_phenotype/meta.txt",head=T)
ksnp = read.table("/public10/home/sci0011/projects/tomato2/37_GS/07_snps_linear/K.txt",head=F)
kindel = read.table("/public10/home/sci0011/projects/tomato2/37_GS/08_indels_linear/K.txt",head=F)
ksv = read.table("/public10/home/sci0011/projects/tomato2/37_GS/09_svs_linear/K.txt",head=F)
# K = list(ksnp= as.matrix(ksnp),kindel = as.matrix(kindel),ksv = as.matrix(ksv))
cov =read.table("/public10/home/sci0011/projects/tomato2/17_cov/01_snps/pheno.cov",head=T)

Y = read.table(paste("/public10/home/sci0011/projects/tomato2/37_GS/01_phenotype/rep",rep,"_fold",fold,".txt",sep=""),head=T)
res = matrix(0,9,972)
res = as.data.frame(res)
res[1:3,1] = "snps"
res[1:3,2] = c("h2","fit","acc")

res[4:6,1] = "indels"
res[4:6,2] = c("h2","fit","acc")

res[7:9,1] = "svs"
res[7:9,2] = c("h2","fit","acc")

colnames(res) = c("marker","type",colnames(Y)[5:974])

for (p in 5:974){
  data = data.frame(y=Y[,p],gid=Y$IID)
  idx = is.na(Y[,p])
  
  colnames(ksnp) = Y$IID
  rownames(ksnp) = Y$IID
  ans <- kin.blup(data=data,geno="gid",pheno="y",K=as.matrix(ksnp))
  res[1,p-2]  = ans$Vg/(ans$Vg + ans$Ve)
  res[2,p-2]  = cor(Y[!idx,p],ans$g[!idx],use="pairwise.complete.obs")
  res[3,p-2]  = cor(phe[idx,p],ans$g[idx],use="pairwise.complete.obs")
  
  
  colnames(kindel) = Y$IID
  rownames(kindel) = Y$IID
  ans <- kin.blup(data=data,geno="gid",pheno="y",K=as.matrix(kindel))
  res[4,p-2]  = ans$Vg/(ans$Vg + ans$Ve)
  res[5,p-2]  = cor(Y[!idx,p],ans$g[!idx],use="pairwise.complete.obs")
  res[6,p-2]  = cor(phe[idx,p],ans$g[idx],use="pairwise.complete.obs")
  
  colnames(ksv) = Y$IID
  rownames(ksv) = Y$IID
  ans <- kin.blup(data=data,geno="gid",pheno="y",K=as.matrix(ksv))
  res[7,p-2]  = ans$Vg/(ans$Vg + ans$Ve)
  res[8,p-2]  = cor(Y[!idx,p],ans$g[!idx],use="pairwise.complete.obs")
  res[9,p-2]  = cor(phe[idx,p],ans$g[idx],use="pairwise.complete.obs")
  
}
write.table(res,paste("rep",rep,"_fold",fold,".rrBLUP.linear.txt",sep=""),col.names = T,row.names = F,quote=F,sep="\t")



# myGS = iGS.eps(X = as.matrix(cov[,3:6]),Y= Y[,c(1,5)],K = K,alg="fs")
# phe1 = myGS$phe[,1]

# source("/public10/home/sci0011/soft/GST/calcAIvar.R")
# source("/public10/home/sci0011/soft/GST/gBLUP.eps.R")
# source("/public10/home/sci0011/soft/GST/iGS.eps.R")
# source("/public10/home/sci0011/soft/GST/Keps.R")
