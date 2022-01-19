# Heritability analysis

## Phenotype data preparion

### Gene expression and metabolite contents

To quantify the expression of all genes, we apply the `Kallisto` software

```shell
kallisto index ./all.gene.cdna.fa.gz -i all.gene.cdna.fa.index
kallisto quant -i ./all.gene.cdna.fa.gz -o ./Result -t 4 -b 100 PATH/Sample_R1.fq.gz PATH/Sample_R2.fq.gz 
```

Normalization of all expression results

```shell
Rscript script/prepare_gene_expression.R
```

To remove potential batch effects and cconfounding factors impacting gene expression, the Probabilistic Estimation of Expression Residuals (`PEER`) method is performed with the top four factors as covariates.

```shell
library(peer)
expr = read.csv('examples/data/expression.csv', header=FALSE)
model = PEER()
PEER_setPhenoMean(model,as.matrix(expr))
PEER_setNk(model,20)
PEER_getNk(model)
PEER_update(model)
factors = PEER_getX(model)
weights = PEER_getW(model)
precision = PEER_getAlpha(model)
residuals = PEER_getResiduals(model)
PEER_setAdd_mean(model, TRUE)
write.table(residuals,"./peer_residuals.tsv",quote=F,row.names=F,sep="\t",col.names=F)
write.table(factors,"./peer_covariates.tsv",quote=F,row.names=F,sep="\t",col.names=F)
```



## Prepare genetic variants 

Extract variants within cis regions

```shell
less expre_sl5_graph.bed | awk '{print $1"\t"$2-50000"\t"$3+50000"\t"$4}' | sort -k2,2n >tmp
sort -k1,1n -k2,2n tmp >sl5_graph_gene_bed
bedtools intersect -a ../332_snp_bed   -b ../sl5_graph_gene_bed -wo  >tmp
bedtools intersect -a ../332_sv_bed   -b ../sl5_graph_gene_bed -wo  >tmp
for i in `less ../gene_id`; do echo " grep \"$i\"  tmp | awk '{print \$4}'|sort | uniq  > $i";done  >para.sh
echo  "~/software/ParaFly -c para.sh  -CPU 60" >1.sh
```
Extract all local variants around a leading significant variants

```shell
bedtools intersect -a ../../combined_leading.bed   -b ../332_bed   -wo >tmp
less tmp|awk '{print $5}'|sort |uniq >gene_id
sed -i 's/_gene//' gene_id
for i in `less ./gene_id`; do echo "grep \"${i}_gene\"  tmp | awk '{print \$4}'|sort | uniq  > $i ";done  >para.sh
echo "~/software/ParaFly -c para.sh  -CPU 64 ">1.sh
sbatch -p amd_io 1.sh
```

Extract variants from a given module genes

```shell
for i in `less ./gene_id`; do echo "cat $i all_varaint |sort |uniq -u  >${i}_peripheral";done  >peripheral.sh
python peripheral.py
```



## Heritability estimation

The analysis refered to the heritability estimation and simulation are followed `LDAK` manual (https://dougspeed.com/ldak/).  

There is a small demo for heritability estimation.


```shell
## prepare genotype in plink format
plink  --vcf test.vcf.gz --recode --out test --double-id
plink  --file test --make-bed --out test
plink --keep ID.kept.txt --bfile test --make-bed --maf 0.005 --out kept.test
plink --bfile kept.test --recode vcf --out kept.test
bgzip kept.test.vcf

# prune the genotype dataset
ldak --bfile kept.test --cut-weights snps --window-prune 0.98

# weighting
ldak --bfile kept.test  --calc-weights-all snps 1 > weights.out 

# calculate kinship
ldak  --calc-kins-direct LDAK-Thin --bfile  kept.test  --weights snps/weights.thin --power -.5 

# estimate heritability
ldak --pheno pheno.txt --mpheno 4  --grm LDAK-Thin  --covar pheno.cov --reml 4 --constrain YES
```



## Genome-wide association study

For the mixed linear model (`MLM`), we used the leave-one-chromosome-out (LOCO) method and the mixed model implemented in `GCTA`. The pruned SNPs were used for the estimation of a kinship matrix (GRM) using LDAK and then performed GWAS using GCTA ((https://yanglab.westlake.edu.cn/software/gcta/))

```shell
## Calculate kinship matrix for each chromosome using LDAK
ldak --bfile input --thin exp --window-prune .5 --window-kb 100

for j in {1..12}; do ldak --calc-kins-direct exp$j --bfile input --ignore-weights YES --power -0.5 --extract exp.in --chr $j; done

for j in {1..12}; do echo "exp$j" >> list.All; done
ldak --add-grm expAll --mgrm list.All
for j in {1..12}; do echo "expAll exp$j" > list.$j; ldak --sub-grm expN$j --mgrm list.$j; done

## Genome-wide association study by chromosome using GCTA
for chr in {1..12}; do gcta64 --mlma --bfile input.chr$chr --pheno pheno.txt --out 4.chr$chr --thread-num 1 --mpheno 4 --grm expN$chr --qcovar pheno.cov ; done
```



## LASSO analysis

Beagle is applied for variation imputation (https://github.com/adrianodemarino/Imputation_beagle_tutorial)

Lassopv software is applied to association study (https://github.com/lingfeiwang/lassopv)

## Prepare genotype

```shell
java -Xmx4g -jar script/beagle.18May20.d20.jar gt=10000.vcf out=10000.inputed
plink --indep-pairwise 100 1 0.7 --vcf 10000.inputed.vcf.gz --out 10000.filtered
plink --vcf 10000.inputed.vcf.gz --extract 10000.filtered.prune.in --make-bed --out 10000.filtered
plink --bfile 10000.filtered --recode vcf --out 10000.filtered
WGS --model vcf --type vcf2GD --file 10000.snps.filtered.vcf --out 10000.snps.dat
```

## Using lassopv package

```shell
## this is a demo
Rscript script/lasso.demo.r
```



## Co-expression network

We applied Weighted Correlation Network Analysis (WGCNA) on the pre-filtered expression data from 332 accessions. The total tutorials followed the page.(https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/)

```R
library(WGCNA)
## checking data for excessive missing values and identification of outlier samples
options(stringsAsFactors = FALSE)  
myData = read.table("expression.data", sep=" ", header=TRUE)  
dim(myData)  
names(myData)  
datExpr = as.data.frame(t(myData[, -c(1)]))  
names(datExpr) = myData$FID  
rownames(datExpr) = names(myData)[-c(1)]  
gsg = goodSamplesGenes(datExpr, verbose = 3)  
if (!gsg$allOK)  
{  
if (sum(!gsg$goodGenes)>0)  
    printFlush(paste("Removing genes:", paste(names(datExpr)[!gsg$goodGenes], collapse = ", ")))  
if (sum(!gsg$goodSamples)>0)  
    printFlush(paste("Removing samples:", paste(rownames(datExpr)[!gsg$goodSamples], collapse = ", ")))  
    datExpr = datExpr[gsg$goodSamples, gsg$goodGenes]  
}  
write.table(names(datExpr)[!gsg$goodGenes], file="removeGene.xls", row.names=FALSE, col.names=FALSE, quote=FALSE)  
write.table(names(datExpr)[!gsg$goodSamples], file="removeSample.xls", row.names=FALSE, col.names=FALSE, quote=FALSE)  
sampleTree = hclust(dist(datExpr), method = "average")
pdf(file = "sampleClustering.pdf", width = 12, height = 9)  
par(cex = 0.6)  
par(mar = c(0,4,2,0))  
plot(sampleTree, main = "Sample clustering", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)  
dev.off()  
save(datExpr, file = "dataInput.RData")  

## automatic network construction and module detection
enableWGCNAThreads(2)  
lnames = load(file ="dataInput.RData")  
powers = c(c(1:10), seq(from = 12, to=20, by=2))  
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)  
sizeGrWindow(9, 5)  
par(mfrow = c(1,2))  
cex1 = 0.9  
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n", main = paste("Scale independence"))  
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1,col="red")  
plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"))  
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
softPower = as.integer(9)

## automatic construction of the gene network
adjacency = adjacency(datExpr, power = softPower) 
TOM = TOMsimilarity(adjacency)
dissTOM = 1-TOM  
geneTree = hclust(as.dist(dissTOM), method = "average")
minModuleSize = as.integer(10)
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM, deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = minModuleSize)
dynamicColors = labels2colors(dynamicMods)  
MEList = moduleEigengenes(datExpr, colors = dynamicColors)  
MEs = MEList$eigengenes
MEDiss = 1-cor(MEs)  
METree = hclust(as.dist(MEDiss), method = "average")  
sizeGrWindow(7, 6)  
plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")  
MEDissThres = as.double(0.15)  
pdf(file = "clusterModuleEigengenes.pdf", width = 7, height = 6)  
plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")  
abline(h=MEDissThres, col = "red")  
dev.off()  

## Identification of modules
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)  
mergedColors = merge$colors  
mergedMEs = merge$newMEs  
pdf(file = "mergedModuleTree.pdf", width = 12, height = 9)  
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)  
dev.off()  
moduleColors = mergedColors  
colorOrder = c("grey", standardColors(50))  
moduleLabels = match(moduleColors, colorOrder) - 1  
MEs = mergedMEs  
write.table(paste(colnames(datExpr), moduleColors, sep = "\t"), file="netcolor2gene.xls", row.names=FALSE, quote=FALSE)  
save(MEs, moduleLabels, moduleColors, geneTree,file = "networkConstruction.RData") 

## visualizing the network of eigengenes
options(stringsAsFactors = FALSE)  
enableWGCNAThreads()  
lnames = load(file ="dataInput.RData")  
lnames = load(file ="networkConstruction.RData")  
nGenes = ncol(datExpr)  
nSamples = nrow(datExpr)  
softpower = as.integer(readline("which softPower: "))  
dissTOM = 1 - TOMsimilarityFromExpr(datExpr, power = softpower)  
nSelect = as.integer(readline("choose a maxGeneNumber(1500): "))  
set.seed(10)  
select = sample(nGenes, size = nSelect)  
selectTOM = dissTOM[select, select]  
selectTree = hclust(as.dist(selectTOM), method = "average")  
selectColors = moduleColors[select]  
plotDiss = selectTOM^7  
diag(plotDiss) = NA  
pdf(file = "networkHeatmap.pdf", width = 15, height = 15)  
TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot")  
dev.off()  
save(dissTOM, file = "dissTOM.RData")  
MEs = moduleEigengenes(datExpr, moduleColors)$eigengenes  
MET = orderMEs(MEs)  
pdf(file = "hubGeneHeatmap.pdf", width = 6, height = 6)  
plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2), plotDendrograms = FALSE, xLabelsAngle = 90)  
dev.off()  

## extract interested modules according color Id
options(stringsAsFactors = FALSE)  
enableWGCNAThreads()  
lnames = load(file ="dataInput.RData")  
lnames = load(file ="networkConstruction.RData")  
lnames = load(file ="dissTOM.RData")  
TOM = 1 - dissTOM  
modules = c("brown", "red");
probes = names(datExpr)
inModule = is.finite(match(moduleColors, modules));
modProbes = probes[inModule];
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
cyt = exportNetworkToCytoscape(modTOM,edgeFile = paste("AS-green-FPKM-Step-by-step-CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""), nodeFile = paste("AS-green-FPKM-Step-by-step-CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),weighted = TRUE,threshold = 0.25,nodeNames = modProbes,nodeAttr = moduleColors[inModule]);

## find hubgene in each module
HubGenes <- chooseTopHubInEachModule(datExpr,moduleColors)
write.table (HubGenes,file = "HubGenes_of_each_module.xls",quote=F,sep='\t')

## according to module id extract link genes with threshold 0.02
color_txt<-read.csv("module.id",head=FALSE)
c<-as.vector(color_txt[,1])
for (modules in c) {
probes = names(datExpr)
inModule = is.finite(match(moduleColors, modules));
modProbes = probes[inModule];
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
cyt = exportNetworkToCytoscape(modTOM,edgeFile = paste("AS-edges-", paste(modules, collapse="-"), ".txt", sep=""), nodeFile = paste("AS-nodes-", paste(modules, collapse="-"), ".txt", sep=""),weighted = TRUE,threshold = 0.02,nodeNames = modProbes,nodeAttr = moduleColors[inModule]);
}
```



# Genomic selection

The package rrBLUP was used for the genomic selection of metabolites. SNPs and InDels with weight larger than zero were used to derive the kinship matrix with the A.mat function implemented in rrBLUP. The prediction accuracy was obtained by performing a five-fold cross-validation with 20 repetitions (https://github.com/cran/rrBLUP).

```shell
## a demo
Rscript script/rrBLUP.r
```

