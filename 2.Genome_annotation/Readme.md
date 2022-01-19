# Genome annotation methods

## Transposable elements (TEs) annotation

The Extensive de-novo TE Annotator (`EDTA`) is applied to annotate TEs which includes `LTRharvest_parallel`,`LTR_FINDER_parallel`,`LRR_retriever`,`TIR-Learner`,`HelitronScanner`,`RepeatModeler`,`RepeatMasker`. This pipeline is integrated  in current directory `./workflow/Snakefile` and the final TE.gff3 is as input file for following gene annotation (**Figure 1**).

## Protein-coding gene annotation

The gene annotation pipeline:

- Homology-based prediction

  ​	83 public RNA-seq dataset including various tissues from BIG, CER, PIM groups are mapped to assembly using `hisat2` and subsequently assembled into transcripts by `stringtie`. Then, `TACO` is applied to merge all stringtie gtf.

  ​	Proteins from SwissProt Viridiplantae, LA2093, DM v6.1, Ath, ITAG4.0 are integrated and remove redundant proteins by `CD-HIT` . ESTs are retrieved from the NCBI.

- Ab initio prediction

  ​	`GeneMark-ET/-ES`,`BRAKER` ,`SNAP` and `Augustus` are applied to predict gene structure.

- Integration homology-based and ab initio prediction

​			The integrated gene models with AED values < 0.5 are retained



The TEs and protein-coding gene annotation pipeline is deposited in current directory `1.Repeat_gene_annotaion`

To configure this workflow, you can modify ``config/config.yaml`` according to needs, following the explanations provided in the file:

* Add samples to `config/samples.tsv`. Only the column `sample` is mandatory, but any additional columns can be added.
* For each sample, add one or more sequencing units (runs, lanes or replicates) to the unit sheet `config/units.tsv`. For each unit, define platform, and either one (column `fq1`) or two (columns `fq1`, `fq2`) FASTQ files (these can point to anywhere in your system).

## Workflow

![image-20220101234140544](/Users/zhiyangzhang/Library/Application Support/typora-user-images/image-20220101234140544.png)
Figure 1. Overview of the genome annotation pipeline. 



## Pooled RNA-seq and assembly model (PRAM) annotation
To discover novel transcripts located in the intergenic region of SL5.0, , which based on on a single step that encompasses all transcript construction, is used. A total of 332 RNA-seq sequences are equipped.

Build SL5.0 genome index by  `STAR`

```shell
STAR --runThreadN 20 \
    --runMode genomeGenerate \
    --genomeDir /your_path/STAR/genome.index \
    --genomeSAindexNbases 13
    --genomeFastaFiles /your_path/SL5.0.genome.fa \
    --sjdbGTFfile /your_path/genome.annotation.gtf \
    --sjdbOverhang 149
```

Map the 332 RNA-seq clean data

```shell
STAR --runThreadN 20 --genomeDir /your_path/STAR/genome.index \
--readFilesCommand gunzip -c --outSAMtype BAM SortedByCoordinate \
--outFileNamePrefix sample --outSAMattrIHstart 0 \
--readFilesIn ./clean_fq/{sample}_1.fq.gz  ./clean_fq/{sample}_2.fq.gz
```

`PRAM` is used to identify the intergenic region (https://bioconductor.riken.jp/packages/3.9/bioc/vignettes/pram/inst/doc/pram.pdf)

```shell
source activate RNA
library("pram")
in_gtf=system.file('extdata/SL5.0.annotation.gtf', package='pram')
bam_data<-read.table("/public/home/zhangzhiyang/anaconda3/envs/RNA/lib/R/library/pram/extdata/399_bam/bam_list",header=F)
pred_out_gtf=tempfile(fileext='.gtf')
pram::runPRAM( in_gtf, in_bamv,pred_out_gtf,method='plst',stringtie='~/anaconda3/envs/RNA/bin/stringtie')
```

Filter the potential transposable elements, according to information from pram.fasta.rexdb.dom.tsv to remove the potential transposable elements

```shell
/public10/home/sci0009/scripts/TEsorter.sh
source /public10/home/sci0009/miniconda3/bin/activate edta
TEsorter pram.fasta   -p 64
```

`Transdecoder` is used to predict the likely coding regions in remains

```shell
TransDecoder.Predict -t target_transcripts.fasta
```



### Gene function annotation

Homolog annotation of 46 tomato genes by using `UniProtKB/SwissProt` database

```shell
~/anaconda3/bin/diamond makedb --in swissprot.fasta  -d swissprot
~/anaconda3/bin/diamond blastp --more-sensitive -d swissprot -q 46_genome_protein.fasta -k 20 -e 0.00001 -o swissprot.xml -f 5

python xml2tab.py
python UniProt2GO_annotate.py idmapping.tb.gz swissprot_sprot.tab swissprot_sprot.go
```

The motifs and domains of predicted genes were predicted using `InterProScan`

```shell
source activate interproscan
for i in `less  /public10/home/sci0010/project/17.function/all_prot/MM/MM.pep.fa.split/id`; do echo "/public10/home/sci0009/software/InterProscan/interproscan-5.48-83.0/interproscan.sh    -i /public10/home/sci0010/project/17.function/all_prot/MM/MM.pep.fa.split/${i} -o tsv/${i}  -goterms -iprloo1kup -pa -f TSV -dp -cpu 2 &";done > run.sh
sbatch -p amd_io run.sh
/public10/home/sci0009/software/InterProscan/interproscan-5.48-83.0/interproscan.sh -i  /public10/home/sci0010/project/17.function/all_prot/SL5/SL5.pep.fa.split/SL5.pep.part_001.fa -o tsv/SL5.pep.part_001.fa  -goterms -iprloo1kup -pa -f TSV -dp -cpu 2
```



