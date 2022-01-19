#!/bin/bash Rscript
require("lars")
require("stats")
lassopv<-function (x,y,normalize=TRUE,H0=c("spherical","normal"),log.p=FALSE,max.predictors=NULL,trace = FALSE,Gram,eps = .Machine$double.eps,max.steps,use.Gram=TRUE) 
{
  dx=dim(x)[2]
  ds=dim(x)[1]
  if(ds!=length(y))
    stop("Incorrect dimensions for x or y. The number of rows in x should match the length of y.")
  if(is.null(max.predictors))
    max.predictors=dx
  else if(max.predictors<=0)
    stop("Incorrect max.predictors")
  
  isdone=rep(F,dx)
  if(dx==1)
  {
    ans=cor(drop(x),y)		
    isdone[1]=T
  }
  else
  {
    #Normalize
    x=t(t(x)-colMeans(x))
    sx=sqrt(colMeans(x**2))
    y=y-mean(y)
    sy=sqrt(mean(y**2))
    if(normalize)
    {
      x=t(t(x)/sx)
      sx=sqrt(colMeans(x**2))
    }
    a=lars(x,y,type='lasso',trace=trace,normalize=F,intercept=F,Gram=Gram,eps=eps,max.steps=max.steps,use.Gram=use.Gram)
    #Crop away ascending lambda values (supposed bug in lars)
    nstep=which(a$lambda[2:length(a$lambda),drop=F]>a$lambda[1:length(a$lambda)-1])
    if(length(nstep)==0)
      nstep=length(a$lambda)
    else
      nstep=nstep[1]
    #lambda values
    lambda=a$lambda[1:nstep,drop=F]/ds
    #Residue
    ypredict=predict.lars(a,x,s=a$lambda[1:nstep,drop=F],mode='lambda')$fit
    if(is.null(dim(ypredict)))
    {
      ypredict=matrix(ypredict,ncol=1)
    }
    yres=y-ypredict
    #Residue variance
    pv=t(t(yres)-colMeans(yres))
    pv=colMeans(pv**2)
    wcl=lambda/sqrt(pv)
    #ans=lambda_i/(sigma_i*sigma_yres)
    ans=rep(0.,dx)
    ntot=0
    for (i in 1:nstep)
    {
      if(!is.finite(wcl[i]))
        next
      pa=unlist(a$actions[i])
      for(j in 1:length(pa))
      {
        if((pa[j]<=0)||isdone[pa[j]])
          next
        ans[pa[j]]=wcl[i]/sx[pa[j]]
        isdone[pa[j]]=T
        ntot=ntot+1
      }
      if(ntot>=max.predictors)
        break
    }
  }
  
  if(H0[1]=="normal")
  {
    ans=pchisq((ans**2)*ds,1,lower.tail=F,log.p=log.p)
  }
  else if(H0[1]=="spherical")
  {
    ans=pt(abs(ans)*sqrt((ds-2)/(1-ans**2)),ds-2,lower.tail=F,log.p=log.p)
    if(log.p)
      ans=ans+log(2)
    else
      ans=ans*2
  }
  else
    stop("Unknown null distribution type H0.")
  if(log.p)
  {
    ans[!isdone]=0
    ans[is.nan(ans)]=-Inf
    #Account for numerical precision limitations
    ans[ans>0]=0
  }
  else
  {
    ans[!isdone]=1
    ans[is.nan(ans)]=0
    #Account for numerical precision limitations
    ans[ans>1]=1
  }
  names(ans)=colnames(x)
  ans
}
`Blink.LDRemoveBlock`<-function(GDneo=NULL,LD=NULL,Porder=NULL,bound=FALSE,model="A",orientation=NULL){
#`Blink.LDRemove`<-function(GDneo=NULL,LD=NULL,Porder=NULL,bound=FALSE,model="A",orientation=NULL){
#Objects: Calculate LD and remove the correlated SNPs
#Authors: Yao Zhou
#Last Update:  03/03/16
	if (model=="D"){
		GDneo=1-abs(GDneo-1)
	}

	GDneo=as.matrix(GDneo)
	if(min(ncol(GDneo),nrow(GDneo))<201) bound=FALSE
	if(orientation=="col"){
		n=nrow(GDneo)
		if(bound){
			GDneo=GDneo[sample(n,200,replace=F),]
		}
	}else{
		n=ncol(GDneo)
		if(bound){
			GDneo=GDneo[,sample(n,200,replace=F)]
		}
		GDneo=t(GDneo)
	}
	# cat("ncol(GDneo) is",ncol(GDneo),"\n")
	corr=cor(GDneo)
	corr[is.na(corr)]=1
	corr[abs(corr)<=LD]=0
	corr[abs(corr)>LD]=1
	Psort=as.numeric(matrix(1,1,ncol(corr)))
	# print(ncol(corr))
	for(i in 2:ncol(corr)){
		p.a=Psort[1:(i-1)]
		p.b=as.numeric(corr[1:(i-1),i])
		index=(p.a==p.b)
		index[(p.a==0)&(p.b==0)]=FALSE
		if(sum(index)!=0) Psort[i]=0
	}
	seqQTN=Porder[Psort==1]
	return(seqQTN)
}

`Blink.LDRemove`<-function(GDneo=NULL,LD=0.7,Porder=NULL,bound=FALSE,model="A",orientation="row",block=1000,LD.num =50){
#Objects: LD remove, especially length(Porder)>10000
#Authors: Yao Zhou
#Last update: 08/15/2016
  GDneo = as.matrix(GDneo)
  SNP.index = apply(GDneo,1,sd)!=0
  GDneo = GDneo[SNP.index,]
  Porder = Porder[SNP.index]
  l = block
	seqQTN=NULL
	lp=length(Porder)
	k=ceiling(lp/l)
	GDneo=as.matrix(GDneo)
	if(min(ncol(GDneo),nrow(GDneo))<201) bound=FALSE
	if(orientation=="col"){
		n=nrow(GDneo)
		if(bound){
			GDneo=GDneo[sample(n,200,replace=F),]
		}
	}else{
		n=ncol(GDneo)
		if(bound){
			GDneo=GDneo[,sample(n,200,replace=F)]
		}
		GDneo=t(GDneo)
	}
	for(i in 1:k){
		bottom=(i-1)*l+1
		up=l*i
		if(up>lp) up = lp
		Porderb=Porder[bottom:up]

		index = seq(bottom:up)
		GDneob = GDneo[,index]
		# cat("i is ",i,"\n")
		# print(length(index))
		seqQTNs = Blink.LDRemoveBlock(GDneo=GDneob,LD=LD,Porder=Porderb,orientation="col",model=model)
		# print(seqQTN)
		seqQTN = append(seqQTN,seqQTNs)
		if(k >1){
		  index1 = which(Porder %in% seqQTN)
		  Porderb = Porder[index1]
		  GDneob = GDneo[,index1]
		  if(length(index1)>1){
		    seqQTN = Blink.LDRemoveBlock(GDneo=GDneob,LD=LD,Porder=Porderb,orientation="col",model=model)
		  }else{
		    seqQTN = Porderb
		  }

		}
		if(LD.num < length(seqQTN)) break
	}
	rm(GDneob,Porderb)
	return(seqQTN)
}
num=4
dsnps = try( read.table(paste("./",num,".dat",sep=""),head=F), silent=TRUE)
if ('try-error' %in% class(dsnps)) {
  dsnps = NULL
  msnps = NULL
}else{
  msnps= read.table(paste("./",num,".filtered.bim",sep=""),head=F)
  msnps = msnps[,c(2,1,4)]
}
dat = dsnps
GM = msnps
if(nrow(dat)>1){
a = Blink.LDRemove(GDneo=dat,LD=0.84,Porder=seq(1:nrow(dat)))
}else{
a=1
}
if (!is.null(dat)){
  y = read.table("pheno.txt",head=T)
  colnames(GM) = c("ID","chr","pos")
  # EM-lasso analysis
  p = try( lassopv(x=t(dat[a,]),y=y[,4]+y[,6]+y[,8],H0=c("spherical"),use.Gram=F), silent=TRUE)
  if (!('try-error' %in% class(p))){
    res = data.frame(GM[a,],pvalue=p)
    write.table(res,paste("./",num,".gwas.txt",sep=""),col.names=T,row.names=F,quote=F,sep="\t")
  }
}

