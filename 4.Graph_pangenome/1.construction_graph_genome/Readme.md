# Graph Pangenome



## SNPs and InDels calling from HiFi reads

The 31 tomato accessions HiFi reads are mapped to SL5.0 using `minimap`  and called variants using `deepvariant`. 

```shell
cat sample.txt | parallel -j 10 'minimap2 -ax map-pb -a -k 19 -O 5,56 -E 4,1 -B 5 -z 400,50 -r 2k --eqx --secondary=no | samblaster -e | samtools sort - -@ 32  -o {}.bam
cat sample.txt | parallel -j 10 'sh DeepVariant.sh {}.bam'
```

DeepVariant.sh

```shell
ID=$1
singularity exec \
  -B /public10/home/sci0011/projects/tomato/genotype/06_snps/:/input \
  -B /public10/home/sci0011/projects/tomato/genotype/06_snps/:/output \
 /public10/home/sci0011/soft/deepvariant-1.0.0.sif \
   --config DeepVariant \
  /opt/deepvariant/bin/run_deepvariant \
   --sample_name=$ID\
  --model_type=PACBIO \
  --ref=/input/SL5.fa \
  --reads=/input/$ID/${ID}.rm.bam \
  --output_vcf=/output/$ID/${ID}.vcf.gz \
  --output_gvcf=/output/$ID/${ID}.g.vcf.gz \
  --intermediate_results_dir=/output/deepvariant_tmp_output/$ID \
  --vcf_stats_report=false \
  --num_shards=32
```

## SNPs and InDels fillter

Retain site depth ranged from 400 to 1,500 using `WGS`

```shell
WGS --model file --type random --file graph31.vcf.gz --headLine 31 --r 0.005 --out random.vcf
WGS --model vcf --type calTotalDP --file random.vcf --out random.dep ## mean depth is 981.3
# R code
## dat = read.table("random.dep",head=F)
## hist(dat$V3,col=c("gray"),breaks=1000,xlim=c(0,2000),xlab="Total depth",main="")
WGS --model vcf --type depthFilterDP --minDepth 400 --maxDepth 1500 --file graph31.vcf.gz --out ../depth/graph31.dep.vcf 
zcat graph31.vcf.gz | head -n 30 > header
```

Retain site quality score more than 20

```shell
cd /public4/home/sc55932/projects/tomato/deepvariant/quality
cat ../depth/graph31.dep.vcf | awk '{n=split($0,a,"\t"); if (a[6] > 19) print $0}' > graph31.qual.vcf
## 16703129 out of 17642861 left (SNPs and InDels)
```

Keep indels size less than 50 bp

```shell
WGS --model vcf --type inDel_len --file ../quality/graph31.qual.vcf --maxLength 50 --out graph31.len.vcf
```

Keep biallelic sites 

```shell
 cd /public4/home/sc55932/projects/tomato/deepvariant/biallelic
 vcftools --vcf ../len/graph31.len.vcf --min-alleles 2 --max-alleles 2 --recode --recode-INFO-all --out graph31
```

Split variants vcf by chromosome

```shell
cd /public4/home/sc55932/projects/tomato/deepvariant/final
cat ../code/chr.txt | parallel -j 13 "WGS --model vcf --type splitByChr --chr {} --file ../biallelic/graph31.recode.vcf --out chr{}.vcf"
```

Transfer variants vcf to plink format

```shell
cd /public4/home/sc55932/projects/tomato/deepvariant/final
cat ../code/chr.txt | parallel -j1 'plink --vcf chr{}.vcf --make-bed --out chr{} --vcf-half-call m --memory 4600'
cat ../code/chr.txt | parallel -j1 "grep 'varaint'  chr{}.log | head -n 1 | awk '{print $1}'"
cat ../code/chr.txt | parallel -j12 "bgzip chr{}.vcf"
```

Caculate the accession heterozygousity

```shell
cd /public4/home/sc55932/projects/tomato/deepvariant/het
cat ../code/chr.txt | parallel -j1 'plink --bfile ../final/chr{} --het --out chr{}'
```



## Structural variants (SVs) detection

To detect SVs using HiFi reads from the 31 accessions. We apply `NGLMR` mapping software and a total of four callers (`Sniffles`,  `SVIM`, `CuteSV`, `PBSV`). The complete pipeline is packed in a SV_Snakfile. There is a small demo in 01_ccs_sv_pipeline. 

1. The input HiFi reads must named as {sample}.ccs.fa

   SAMPLE_INDEX = {"sampleA","sampleB"}

2. Reference store
   INDEX_REF = "01_raw_data/ref/sv_refer.fa"

```shell
source /public/agis/huangsanwen_group/chenglin/softwares/miniconda3/bin/activate base
snakemake -s SV_Snakefile -j 2 --stats snakejob.stats >&2 2>>snakejob.log
```

- workflow

  ![image-20220102154738977](/Users/zhiyangzhang/Library/Application Support/typora-user-images/image-20220102154738977.png)
  Figure 1. The pipeline of SVs indentification.


## Structure variants (SVs) filtering

We retain variants with a “pass” flag and read depth of at least three. Deletions ranging from 51 bp to 100 kb in length, and insertions ranging from 51 bp to 20 kb in length are retained by `WGS` 

```shell
cd /public4/home/sc55932/projects/tomato/SVs/03_vcf/01_sniffles
## maxLength set 20000 for SV called from reads;
cat ../06_summary/ID.txt | parallel 'WGS --model vcf --type SVfilter_reads --maxLength 20000 --file {}.sniffles.vcf --out {}.filtered.vcf' 

cd /public4/home/sc55932/projects/tomato/SVs/03_vcf/02_cuteSV
cat ../06_summary/ID.txt | parallel 'WGS --model vcf --type SVfilter_reads --maxLength 20000 --file {}.cuteSV.vcf --out {}.filtered.vcf'

cd /public4/home/sc55932/projects/tomato/SVs/03_vcf/03_svim
cat ../06_summary/ID.txt | parallel 'WGS --model vcf --type SVfilter_reads --maxLength 20000 --file {}/variants.vcf --out {}/filtered.vcf'

cd /public4/home/sc55932/projects/tomato/SVs/03_vcf/04_pbsv
cat ../06_summary/ID.txt | parallel 'WGS --model vcf --type SVfilter_reads --maxLength 20000 --file {}/pbsv.vcf --out {}/filtered.vcf'
```

## Detect SVs based on assembly methods

To detect SVs more than 100kb, we apply `Assemblytics` software

```shell
cd /public4/home/sc55932/projects/tomato/SVs/code
cat ID.txt | parallel '/public4/home/sc55932/soft/SURVIVOR/Debug/SURVIVOR bedtovcf ../assemblies_SV/{}/03_SV_result/03_assemblytics/{}/{}.Assemblytics_structural_variants.bed DEL ../assemblies_SV/{}/03_SV_result/03_assemblytics/{}/{}.DEL.vcf'

cat ID.txt | parallel '/public4/home/sc55932/soft/SURVIVOR/Debug/SURVIVOR bedtovcf ../assemblies_SV/{}/03_SV_result/03_assemblytics/{}/{}.Assemblytics_structural_variants.bed INS ../assemblies_SV/{}/03_SV_result/03_assemblytics/{}/{}.INS.vcf'

cat ID.txt | parallel 'WGS --model vcf --type bed2vcf --file ../assemblies_SV/{}/03_SV_result/03_assemblytics/{}/{}.Assemblytics_structural_variants.bed --file2 /public4/home/sc55932/data/tomato/ref/SL5.0/SL5.0_chr_number.fa --file3 /public4/home/sc55932/projects/tomato/assembly/gala/{}/08_final/{}.contigs.fasta --out ../assemblies_SV/{}/03_SV_result/03_assemblytics/{}/{}.vcf'
```
## Structural variants (SVs) stats 

```shell
##### 01_sniffles
cd /public4/home/sc55932/projects/tomato/SVs/03_vcf/01_sniffles
cat ../06_summary/ID.txt | parallel "grep -vP '^0\t' {}.sniffles.vcf | grep 'PASS' | grep -o -P '(?<=SVLEN=).*?(?=;)' > {}.len.txt"
cat cat ZY*.len.txt > ../06_summary/01_length/01_sniffle.len.txt
##### 02_cuteSV
cd /public4/home/sc55932/projects/tomato/SVs/03_vcf/02_cuteSV
cat ../06_summary/ID.txt | parallel "grep -vP '^0\t' {}.cuteSV.vcf | grep 'PASS' | grep -o -P '(?<=SVLEN=).*?(?=;)' > {}.len.txt"
cat ZY*.len.txt > ../06_summary/01_length/02_cuteSV.len.txt
##### 03_svim
cd /public4/home/sc55932/projects/tomato/SVs/03_vcf/03_svim
cat ../06_summary/ID.txt | parallel "grep -vP '^0\t' {}/variants.vcf | grep 'PASS' | grep -o -P '(?<=SVLEN=).*?(?=;)' > {}.len.txt"
cat ZY*.len.txt > ../06_summary/01_length/03_svim.len.txt
##### 04_pbsv
cd /public4/home/sc55932/projects/tomato/SVs/03_vcf/04_pbsv
cat ../06_summary/ID.txt | parallel "grep -vP '^0\t' {}/pbsv.vcf | grep 'PASS' | grep -o -P '(?<=SVLEN=).*?(?=;)' > {}.len.txt"
cat ../06_summary/ID.txt | parallel "grep -vP '^0\t' {}/pbsv.vcf | grep 'PASS' | grep -o -P '(?<=SVLEN=).*?(?=\t)' >> {}.len.txt"
cat ZY*.len.txt > ../06_summary/01_length/04_pbsv.len.txt
###### consensus
cd /public4/home/sc55932/projects/tomato/SVs/04_consensus_vcf
cat ../03_vcf/06_summary/ID.txt | parallel "grep -vP '^0\t' {}.consensus.vcf | grep 'PASS' | grep -o -P '(?<=SVLEN=).*?(?=;)' > {}.len.txt"
cat ZY*.len.txt > ../03_vcf/06_summary/01_length/05_consensus.len.txt
```
## Structural variants (SVs) merge

For the 31 accessions with SVs from the five callers, we merge all SVs shorter than 100 kb with `SURVIVOR`  using a maximum allowed distance of 1 kb, only to report calls supported by two callers and agreed-upon regarding the type of variant. SVs longer than 100 kb dectected by `Assemblytics` are kept.

```shell
cd /public4/home/sc55932/projects/tomato/SVs/code
while read line;do echo "/public4/home/sc55932/projects/tomato/SVs/assemblies_SV/$line/03_SV_result/03_assemblytics/$line/$line.vcf" > ../05_sample_list/$line.list; echo "/public4/home/sc55932/projects/tomato/SVs/03_vcf/01_sniffles/$line.filtered.vcf" >> ../05_sample_list/$line.list; echo "/public4/home/sc55932/projects/tomato/SVs/03_vcf/02_cuteSV/$line.filtered.vcf" >> ../05_sample_list/$line.list; echo "/public4/home/sc55932/projects/tomato/SVs/03_vcf/03_svim/$line/filtered.vcf" >> ../05_sample_list/$line.list; echo "/public4/home/sc55932/projects/tomato/SVs/03_vcf/04_pbsv/$line/filtered.vcf">> ../05_sample_list/$line.list; done < ID.txt

while read line;do /public4/home/sc55932/soft/SURVIVOR/Debug/SURVIVOR merge ../05_sample_list/$line.list 1000 2 1 0 0 50  ../06_merge_five_caller/$line.five_caller.vcf; done < ID.txt

## subtract long insertion
cat ID.txt | parallel 'WGS --model vcf --type SVfilter_long --file ../assemblies_SV/{}/03_SV_result/03_assemblytics/{}/{}.vcf --out ../assemblies_SV/{}/03_SV_result/03_assemblytics/{}/{}.long.vcf'
```
## Filter merge vcf and prepare input variation vcf for vg

There are some disharmony contents from merge vcf files constructed by `SURVIVOR`, we used the custome scripts (`correct_for_vg.py`) to correct these and prepare input variation vcf for vg

```shell
cat ../code/ID.txt | parallel 'python3 correct_for_vg.py  SL5.0_chr_number.fa {}.five_caller.vcf {}.five_caller.vcf_vg_out'

cat ../code/ID.txt | parallel 'WGS --model vcf --type mergeSVs --file /public4/home/sc55932/projects/tomato/SVs/06_seq_solved/{}.five_caller.vcf_vg_out --file2 /public4/home/sc55932/projects/tomato/SVs/assemblies_SV/{}/03_SV_result/03_assemblytics/{}/{}.long.vcf --out /public4/home/sc55932/projects/tomato/SVs/07_merge_long_ins/{}.all.vcf'
```





## Construction of tomato graph genome (TGG)

31 HiFi accessions tomato SVs vcf

```shell
## 214401 SVs
cd /public10/home/sci0011/projects/tomato/genotype/01_graph31
ln -s ../../vg/vcf/graph31.vcf.gz* .
WGS --model vcf --type cleanSVs --file graph31.vcf.gz --out graph31.final.vcf # 214395 left
bgzip graph31.final.vcf
tabix graph31.final.vcf.gz
```

As the publicly 100 tomato SVs (https://solgenomics.net/projects/tomato100) are identified by different version of the reference genome (SL4.0), we transformed the coordinates to SL5.0 using `liftover`

```shell
Rscript filter.R
awk '! /\#/' tomato100.filtered.vcf | awk '{print $1"\t"($2-1)"\t"($2+length($4)-1)"\t"$3}' > tomato100.filtered.bed
liftOver tomato100.filtered.bed SL4ToSL5.over.chain.gz tomato100.filtered.lifted.bed tomato100.filtered.unlifted.bed
python3 liftVcf.py -v tomato100.filtered.vcf -l tomato100.filtered.lifted.bed -o tomato100.filtered.lifted.vcf -r SL5.fa
python removeInfo.py -v tomato100.filtered.lifted.vcf -o tomato100.clean.vcf  -p tomato100
bcftools norm -m -both tomato100.clean.vcf | bcftools norm -d none -c x --fasta-ref SL5.0.fa | bcftools sort | bgzip > tomato100.norm.vcf.gz
WGS --model vcf --type cleanSVs --file tomato100.norm.vcf.gz --out tomato100.final.vcf
bgzip tomato100.final.vcf
tabix tomato100.final.vcf.gz
```

31 HiFi accessions tomato SNPs and InDels vcf

 ````shell
## after filtering, 16201825 SNPs and InDels
cd /public10/home/sci0011/projects/tomato/genotype/03_snps
ln -s ../../snps/graph31.clean.vcf.gz* .
 ````

SVs from the 31 accessions with HiFi reads and previously identified SVs from the 100 tomatoes are merged

```shell
cd /public10/home/sci0011/projects/tomato/genotype/04_merge
bcftools merge ../01_graph31/graph31.final.vcf.gz ../02_tomato100/tomato100.final.vcf.gz | bcftools norm -m -any -N | bcftools norm -d none --fasta-ref ~/data/ref/SL5.0/SL5.0_chr_number.fa  | bcftools sort | bgzip > graphSV.merged.vcf.gz ## 312852 lines
tabix -f graphSV.merged.vcf.gz
```

Remove sites on the same positon

```shell
cd /public10/home/sci0011/projects/tomato/genotype/05_rmdup
WGS --model vcf --type dupPos --file ../04_merge/graphSV.merged.vcf.gz --out graphSV.vcf ## 256783 SVs; 42388 added
WGS --model vcf --type dupPos --file ../04_merge/graphAll.merged.vcf.gz --out graphAll.vcf #  16440660 variants;
bgzip graphSV.vcf
bgzip graphAll.vcf
tabix -f graphSV.vcf.gz
tabix -f graphAll.vcf.gz
```


## Construction of the graph genome TGG1.0

Remap 20 diverse accessions (>30X) to dedup merged SVs followed the protocols (https://github.com/vgteam/giraffe-sv-paper/blob/master/scripts/sv/remap-to-dedup-merged-svs)

```shell
Rscript findDups.R 
snakemake -s dedup_Snakefile --cores 64
grep "#" called_vcf/cluster-1-*.vcf >> called_merged.vcf
grep "#" -vh called_vcf/*vcf >> called_merged.vcf
Rscript mergeRemapResults.R
```

Final SVs vcf merged with 31 HiFi accessions tomato SNPs and InDels vcf

```shell
cd /public10/home/sci0011/projects/tomato/genotype/04_merge
WGS --model vcf --type dupPosSNPs --file ../03_snps/graph31.clean.vcf.gz --file2 graphSV.merged.vcf --out SNPs.vcf # 16201825
bgzip SNPs.vcf 
tabix -f SNPs.vcf.gz
bcftools merge SNPs.vcf.gz  graphSV.merged.vcf.gz | bcftools norm -m -any -N | bcftools norm -d none --fasta-ref ~/data/ref/SL5.0/SL5.0_chr_number.fa  | bcftools sort | bgzip > graphAll.merged.vcf.gz ##  lines
```

To build a variantion graph `TGG1.0` from a reference FASTA file and variants in a VCF file.  A modified `Snakefile`  and `config.yaml`  are download from the `vgteam` (https://github.com/vgteam/vg_snakemake). The correponding graph index `gcsa`, `xg`,`snarls`,`gbwt` files will generate 

```shell
cd /public10/home/sci0011/projects/tomato/genotype/06_snps
snakemake -s vg_Snakefile --configfile config.yaml --config graph="SL5-graphAll" --resources mem_mb=230000 --cores 64 bam

```



## Construction of the graph genome TGG1.1

Short reads from 706 tomato accessions (>6X) are mapped to TGG1.0 with `vg giraffe`

```shell
cat 706.sample.txt | parrell -j 10 | 'vg giraffe -p -t 64 --sample {} -m SL5-graphAll.k39.w15.N16.min -d SL5-graphAll.dist --gbwt-name SL5-graphAll.N16.gbwt -x SL5-graphAll.xg -N {} -f {}/{}_1.fq.gz -f {}/{}_2.fq.gz -o bam > {}/{}-SL5-graphAll.giraffe39k15w16N.bam 2> logs/{}-SL5-graphAll-giraffe39k15w16N.log.txt'
```

 SNP and InDels are called using `DeepVariants` with the NGS model

```shell
#!/bin/bash
while read line
do
sbatch -N 1 -p amd_io ./deepvariant.sh $line
./detect.sh
sleep 1s
done < step3.txt
```

detect.sh

```shell
#!/bin/bash
n=$(echo $(squeue | wc -l))
while [ $n -gt 5 ];
do
sleep 60s
n=$(echo $(squeue | wc -l))
done
```

deepvariant.sh

```shell
#!/bin/bash
ID=$1
source /public10/soft/modules/cn-module.sh
module load singularity/3.5.3-kd
mkdir -p deepvariant_tmp_output/$ID
for i in {0..31}
do
singularity exec   -B /public10/home/sci0011/projects/tomato/07_vg_snps/:/input/   -B /public10/home/sci0011/projects/tomato/07_vg_snps/:/output/  /public10/home/sci0011/soft/deepvariant-1.0.0.sif  /opt/deepvariant/bin/make_examples --mode calling --ref /input/SL5.fa --reads /input/${ID}/${ID}-SL5-graphAll.giraffe39k15w16N.rm.bam --examples /output/deepvariant_tmp_output/${ID}/make_examples.tfrecord@32.gz --gvcf /output/deepvariant_tmp_output/${ID}/gvcf.tfrecord@32.gz --sample_name ${ID} --task $i &
done
wait
echo "make example done" > /public10/home/sci0011/projects/tomato/07_vg_snps/deepvariant_tmp_output/${ID}/example.done
```

### SNPs filtering for 706 samples

According to above criterias, SNPs and InDels in 706 accessions are filltered using `WGS`

```shell
cd /public10/home/sci0011/projects/tomato/10_snps/01_706/00_code

cat chr.txt | parallel -j 13 'WGS --model vcf --type depthFilterDP --minDepth 3000 --maxDepth 9000 --file  ../01_raw/chr{}.vcf.gz --out ../02_depth/chr{}.vcf'
./bgzip.sh 02_depth
cd ../00_code

cat chr.txt | parallel -j 13 'WGS --model vcf --type qualityFilter --threshold 20 --file ../02_depth/chr{}.vcf.gz --out ../03_quality/chr{}.vcf'
./bgzip.sh 03_quality
cd ../00_code

cat chr.txt | parallel -j 13 'WGS --model vcf --type inDel_len --file ../03_quality/chr{}.vcf.gz --out ../04_indel_size/chr{}.vcf'
./bgzip.sh 04_indel_size
cd ../00_code

cat chr.txt | parallel -j 13 'WGS --model vcf --type nameSNPs --file ../04_indel_size/chr{}.vcf.gz --out ../05_rename/chr{}.vcf'
./bgzip.sh 05_rename
cd ../00_code

cat chr.txt | parallel -j 13 'plink --vcf ../05_rename/chr{}.vcf.gz --maf 0.005 --biallelic-only --geno 0.6 --make-bed --vcf-half-call m --memory 2000 --out ../06_plink/chr{}'

cat chr.txt | parallel -j 13 'plink --vcf ../04_indel_size/chr{}.vcf.gz --maf 0.005 --biallelic-only --geno 0.6 --make-bed --vcf-half-call m --memory 2000 --out ../06_plink/chr{}'

```



### Genotypes of Structural variants (SVs) in 706 accessions

Genotypes of SVs for the 706 accessions were called by `Paragraph` using default parameters (https://github.com/Illumina/paragraph). 

```shell
#!/bin/bash
sample=$1
mkdir -p $sample
cd $sample
multigrmpy.py -i /public10/home/sci0011/projects/tomato2/08_paragraph/02_vcf/SV.paragraph.vcf -m /public10/home/sci0011/projects/tomato2/08_paragraph/01_bam/samples_${sample}.txt -r ~/data/ref/SL5.0/SL5.0_chr_number.fa -o . --threads 64
tabix -f genotypes.vcf.gz
```

We futher process the 706 accessions vcf with merging, normalization, renaming, maf, missing rate using `WGS`.

```shell
cd /public10/home/sci0011/projects/tomato/genotype2/02_svs/07_all
cd  01_ori
bcftools merge ()  | bgzip > sv706.vcf.gz

cd ../02_norm
WGS --model vcf --type normVariant --file ../01_ori/sv706.vcf.gz --out sv706.vcf
bgzip sv706.vcf

cd ../03_rename
WGS --model vcf --type nameSNPs --file ../02_norm/sv706.vcf.gz --out sv706.vcf
bgzip sv706.vcf

cd ../04_plink
plink --vcf ../03_rename/sv706.vcf.gz --make-bed --out sv706

## filter with maf < 0.01 and missing rate > 0.2, only 61434 left
cd ../05_filter/
plink --bfile ../04_plink/sv706 --maf 0.01 --geno 0.2 --make-bed --out sv706 ## 61434 left

cd ../06_maf/
plink --bfile ../04_plink/sv706 --freq --out sv706
plink --bfile ../05_filter/sv706 --freq --out sv706.filtered
```





# Benchmark study of graph simulation

To be as consistent as possible with empirical results, we simulated the short-reads with with different coverage (5x, 10x, 15x, 20x, and 25x) followed vg team code. (https://github.com/vgteam/vg)


```shell
cat ID.txt | /public10/home/sci0010/anaconda3/bin/parallel  -j5  'vg sim -r -n 13333333   -a -s 12345 -p 570 -v 165 -i 0.00029 -x {}.xg -g {}.gbwt --sample-name ref   -F SRR7279644_1.fastq.gz  -F SRR7279644_2.fastq.gz | vg annotate -p -x {}.xg -a - | vg view -X -a  - | /public10/home/sci0010/anaconda3/bin/pigz > 5X_{}'

cat ID.txt | /public10/home/sci0010/anaconda3/bin/parallel  -j8  'vg sim -r -n 26666666   -a -s 12345 -p 570 -v 165 -i 0.00029 -x {}.xg -g {}.gbwt --sample-name ref   -F SRR7279644_1.fastq.gz  -F SRR7279644_2.fastq.gz | vg annotate -p -x {}.xg -a - | vg view -X -a  - | /public10/home/sci0010/anaconda3/bin/pigz > 10X_{}'

cat ID.txt | /public10/home/sci0010/anaconda3/bin/parallel  -j5  'vg sim -r -n 40000000   -a -s 12345 -p 570 -v 165 -i 0.00029 -x {}.xg -g {}.gbwt --sample-name ref   -F SRR7279644_1.fastq.gz  -F SRR7279644_2.fastq.gz | vg annotate -p -x {}.xg -a - | vg view -X -a  - | /public10/home/sci0010/anaconda3/bin/pigz > 15X_{}'

cat ID.txt | /public10/home/sci0010/anaconda3/bin/parallel  -j5  'vg sim -r -n 53333333   -a -s 12345 -p 570 -v 165 -i 0.00029 -x {}.xg -g {}.gbwt --sample-name ref   -F SRR7279644_1.fastq.gz  -F SRR7279644_2.fastq.gz | vg annotate -p -x {}.xg -a - | vg view -X -a  - | /public10/home/sci0010/anaconda3/bin/pigz > 20X_{}'

cat ID.txt | /public10/home/sci0010/anaconda3/bin/parallel  -j5  'vg sim -r -n 66666666   -a -s 12345 -p 570 -v 165 -i 0.00029 -x {}.xg -g {}.gbwt --sample-name ref   -F SRR7279644_1.fastq.gz  -F SRR7279644_2.fastq.gz | vg annotate -p -x {}.xg -a - | vg view -X -a  - | /public10/home/sci0010/anaconda3/bin/pigz > 25X_{}'
```

We simulate the genome with known SNPs using custom scripts and further simulated InDels and SVs using `simuG`

```shell
```

SNPs and InDels of linear genome were mapped and called by `bwa` and `DeepVariant` respectively. SVs of linear genome was called same as previous pipeline. To compare the performance between linear and graph genome. We used `hap.py` and `truvari` software

Use `hap.py `

```shell
#!/bin/bash
while read line
do
for cov in 5 10 15 20 25
do
hap.py ../02_genome_simulated/$line.indel.vcf.gz ../05_NGS_calling_DeepVariant/cov$cov.$line/cov$cov.$line.indels.vcf.gz -o cov$cov.$line.indel -r /public10/home/sci0011/projects/tomato2/35_simulation_calling/05_NGS_calling_DeepVariant/SL5.0_chr_number.fa --unhappy --threads 6 --no-roc --no-json -T chr.bed &
hap.py ../01_vcf/$line.snp.recode.vcf.gz ../05_NGS_calling_DeepVariant/cov$cov.$line/cov$cov.$line.snps.vcf.gz -o cov$cov.$line.snps -r /public10/home/sci0011/projects/tomato2/35_simulation_calling/05_NGS_calling_DeepVariant/SL5.0_chr_number.fa --unhappy --threads 6 --no-roc --no-json -T chr.bed &
done
wait
done < sample.txt
```

 Use `truvari`

```shell
#!/bin/bash
while read line
do
for cov in 5 10 15 20 25
do
truvari bench -b ./../../02_genome_simulated/$line.sv.only.vcf1.gz -c ../04_para/vcf/cov$cov.$line.vcf.gz -f ../../05_NGS_calling_DeepVariant/SL5.0_chr_number.fa -o cov$cov.$line -r 2000 -p 0.8 -P 0.7 -O 0.5 --sizemax 200000 -t --multimatch --includebed ../06_eval_ngs/cov$cov.$line.bed &
done
wait
done < sample.txt
```

