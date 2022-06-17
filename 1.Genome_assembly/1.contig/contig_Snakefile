

```shell
#configfile: "config.yaml"
########################################################
## This is for diploid genome assembly (homozygous species)
## Author: Yao Zhou
## Last updated: 12-17-2020
## Contact: Yao Zhou (zhouyao@caas.cn)
########################################################

SAMPLE_INDEX = {"testSample"}
THREADS = "64"

rule all:
    input:
        expand("08_final/{sample}.contigs.fasta",sample=SAMPLE_INDEX)

rule step1:
    input:
       "draft_names_paths.txt" 
    output:
       "01_comparison/asm1vsasm2.paf",
       "01_comparison/asm1vsasm3.paf",
       "01_comparison/asm2vsasm1.paf",
       "01_comparison/asm2vsasm3.paf",
       "01_comparison/asm3vsasm1.paf",
       "01_comparison/asm3vsasm2.paf"
shell:
    """
    comp draft_names_paths.txt 
sed "s/comparison/01_comparison/g" -i draft_comp.sh
    sed "s/^minimap2/minimap2 -t 64/g" -i draft_comp.sh
    sed "1 i#\!\/bin\/bash" -i draft_comp.sh
    sh draft_comp.sh > draft_comp.out 2> draft_comp.err
    """

rule step2:
  input:
    "01_comparison/asm1vsasm2.paf",
    "01_comparison/asm1vsasm3.paf",
    "01_comparison/asm2vsasm1.paf",
    "01_comparison/asm2vsasm3.paf",
    "01_comparison/asm3vsasm1.paf",
    "01_comparison/asm3vsasm2.paf"
  output:
    "02_gathering/gathering_asm1_cuts.txt",
    "02_gathering/gathering_asm2_cuts.txt",
    "02_gathering/gathering_asm3_cuts.txt"
  shell:
    """
      mdm 01_comparison 3 -f 02_gathering
      cd 02_gathering
      mv 02_gathering_asm1_cuts.txt gathering_asm1_cuts.txt
      mv 02_gathering_asm2_cuts.txt gathering_asm2_cuts.txt
      mv 02_gathering_asm3_cuts.txt gathering_asm3_cuts.txt
      mv 02_gathering_asm1.txt gathering_asm1.txt
      mv 02_gathering_asm2.txt gathering_asm2.txt
      mv 02_gathering_asm3.txt gathering_asm3.txt
    """
rule step3:
  input:
    "draft_names_paths.txt",
    "02_gathering/gathering_asm1_cuts.txt",
    "02_gathering/gathering_asm2_cuts.txt",
    "02_gathering/gathering_asm3_cuts.txt"
  output:
    "03_new_genome/new_asm1.fa",
    "03_new_genome/new_asm2.fa",
    "03_new_genome/new_asm3.fa",
    "03_new_genome/new_draft_names_paths.txt"
  run:
    shell("mkdir 03_new_genome")
    shell("cd 03_new_genome")
    shell("newgenome ../draft_names_paths.txt ../02_gathering")

rule step4:
  input:
    "03_new_genome/new_draft_names_paths.txt"
  output:
    "04_comparison/new_asm1vsnew_asm2.paf",
    "04_comparison/new_asm2vsnew_asm1.paf",
    "04_comparison/new_asm1vsnew_asm3.paf",
    "04_comparison/new_asm2vsnew_asm3.paf",
    "04_comparison/new_asm3vsnew_asm1.paf",
    "04_comparison/new_asm3vsnew_asm2.paf"
  shell:
    """
    cd 03_new_genome
    comp new_draft_names_paths.txt 
    sed "s/comparison/..\/04_comparison/g" -i draft_comp.sh
    sed "s/^minimap2/minimap2 -t 50/g" -i draft_comp.sh
    sed "1 i#\!\/bin\/bash" -i draft_comp.sh
    sh draft_comp.sh > draft_comp.out 2> draft_comp.err
    """
rule step5:
  input:
    "03_new_genome/new_asm3.fa"
  output:
    "03_new_genome/new_asm3.fa.sa"
  run:
    shell("cd 03_new_genome")
    shell("bwa index new_asm3.fa")

rule step6:
  input:
    "04_comparison/new_asm1vsnew_asm2.paf",
    "04_comparison/new_asm1vsnew_asm3.paf",
    "04_comparison/new_asm2vsnew_asm1.paf",
    "04_comparison/new_asm2vsnew_asm3.paf",
    "04_comparison/new_asm3vsnew_asm1.paf",
    "04_comparison/new_asm3vsnew_asm2.paf"
  output:
    "06_mapping/ccm.done"
  shell:
    """
    ccm 04_comparison 3 -f 05_scaffold
    cd 05_scaffold
    mv 05_scaffold_new_asm1.scaff scaffold_new_asm1.scaff
    mv 05_scaffold_new_asm2.scaff scaffold_new_asm2.scaff
    mv 05_scaffold_new_asm3.scaff scaffold_new_asm3.scaff
    cd ..
    mkdir -p 06_mapping
    cd 06_mapping
    echo "ccm fnished" > ccm.done
    """
rule step7:
    input:
        "03_new_genome/new_asm3.fa.sa",
    output:
        "06_mapping/new_asm3.bam",
        "06_mapping/new_asm3.bam.bai"
    shell:
      """
      mkdir -p 06_mapping
      bwa mem -t 64 03_new_genome/new_asm3.fa /public4/home/sc55932/data/tomato/hifi/{sample}.fq.gz | samtools sort -@ 64 | samtools view -Sb > 06_mapping/new_asm3.bam
      cd 06_mapping
      samtools index new_asm3.bam
      """
rule step8:
    input:
        "06_mapping/new_asm3.bam",
        "06_mapping/new_asm3.bam.bai",
        "06_mapping/ccm.done"
    output:
        "07_new_assembly/all.fa",
        "07_new_assembly/all.fa.fai"
    shell:
      """
        mkdir -p 07_new_assembly
        cd 06_mapping
        samtools view -H new_asm3.bam |grep "SQ"|cut -f 2|cut -d : -f 2 > contig_names
        mkdir -p bams
        cd bams
        cat ../contig_names | parallel -j 64 'samtools view -b ../new_asm3.bam {} > {}.bam'
        cat ../contig_names | parallel -j 64 'samtools view {}.bam | cut -f 1 > {}.bam.read_names'
        grep "^scaff_" ../../05_scaffold/scaffolds_new_asm3.scaff | awk '{s="cat";n=split($0,a,"\t"); for(i=2;i<=(n-2);i++) s=s" "a[i]".bam.read_names";s=s" > "a[1]".con.read_names"; print s}' > cat.sh
        sh cat.sh
        ls *.con.read_names | sed 's/.con.read_names//g' > groups.txt
        cat groups.txt | parallel -j 64  'WGS --model fastq --type subtract --threshold 100 --file /public4/home/sc55932/data/tomato/hifi/{sample}.fq.gz --file2 {}.con.read_names --out {}.con.read.fq'
        cat groups.txt | parallel -j 2 "hifiasm -o ../../07_new_assembly/{} -t32 {}.con.read.fq 2> {}.asm.log"
        cd ../../07_new_assembly
        ls *.p_ctg.gfa | sed 's/.p_ctg.gfa//g' > groups.txt
        cat  groups.txt | parallel -j 64 'gfatools gfa2fa {}.p_ctg.gfa > {}.p_ctg.fasta'
        
rule step9:
    input:
        "07_new_assembly/all.fa",
        "07_new_assembly/all.fa.fai"
    output:
        "07_new_assembly/final.200K.fasta",
        "07_new_assembly/final.less_200K.fasta"
    shell:
      """
        cd 07_new_assembly
        cat all.fa.fai|sort -k2n | awk '{n=split($0,a,"\t");if (a[2] > 1000000) print a[1]}' > all.larger_than_1M.txt
        WGS --model fasta --type getChrs --file all.fa --file2 all.larger_than_1M.txt --out all.GT1M.fasta
        minimap2 -x asm5 -t 64 all.GT1M.fasta all.GT1M.fasta > all.GT1M.paf 2> all.GT1M.paf.log
        WGS --model paf --type removeHS --file all.GT1M.paf --out all.GT1M.paf.txt
        WGS --model fasta --type getChrs --file all.GT1M.fasta --file2 all.GT1M.paf.txt --out final.GT1M.fasta
        cat all.fa.fai|sort -k2n | awk '{n=split($0,a,"\t");if (a[2] < 1000001) print a[1]}' > all.smaller_than_1M.txt
   			WGS --model fasta --type getChrs --file all.fa --file2 all.smaller_than_1M.txt --out all.LT1M.fasta
        minimap2 -x asm5 -t 64 all.LT1M.fasta all.LT1M.fasta > all.LT1M.paf 2> all.LT1M.paf.log
        WGS --model paf --type removeHS --file all.LT1M.paf --out all.LT1M.paf.txt
        WGS --model fasta --type getChrs --file all.LT1M.fasta --file2 all.LT1M.paf.txt --out final.LT1M.pre.fasta
        minimap2 -x asm5 -t 64 final.GT1M.fasta final.LT1M.pre.fasta > GT1M_LT1M.paf 2> GT1M_LT1M.paf.log
        WGS --model paf --type removeHS --file GT1M_LT1M.paf --out GT1M_LT1M.txt
        WGS --model fasta --type getChrs --file all.fa --file2 GT1M_LT1M.txt --out GT1M_LT1M.fasta

        samtools faidx GT1M_LT1M.fasta
        cat GT1M_LT1M.fasta.fai|sort -k2n | awk '{n=split($0,a,"\t");if (a[2] > 200000) print a[1]}' > final.200K.txt
        WGS --model fasta --type getChrs --file GT1M_LT1M.fasta --file2 final.200K.txt --out final.200K.fasta
        cat GT1M_LT1M.fasta.fai|sort -k2n | awk '{n=split($0,a,"\t");if (a[2] < 200001) print a[1]}' > final.less_200K.txt
        WGS --model fasta --type getChrs --file GT1M_LT1M.fasta --file2 final.less_200K.txt --out final.less_200K.fasta
      """ 


rule step10:
    input:
        "07_new_assembly/final.200K.fasta",
        "07_new_assembly/final.less_200K.fasta",
        "/public4/home/sc55925/database/plastid/plastid.fasta",
        "/public4/home/sc55925/database/mitochondrion/mitochondrion.fasta"
    output:
        "08_final/{sample}.contigs.fasta"
    shell:
      """
        cd 07_new_assembly
        minimap2 -x asm20 -t 64 /public4/home/sc55925/database/plastid/plastid.fasta final.less_200K.fasta > toChl.paf 2> toChl.paf.log
        WGS --model paf --type removeVsRef --file toChl.paf --out toChl.txt
        minimap2 -x asm20 -t 64 /public4/home/sc55925/database/mitochondrion/mitochondrion.fasta final.less_200K.fasta > toMt.paf 2> toMt.paf.log
        WGS --model paf --type removeVsRef --file toMt.paf --out toMt.txt
        cat toChl.txt toMt.txt | sort | uniq -d | sort -n > plastid.txt
        WGS  --model fasta --type getChrs --file final.less_200K.fasta --file2 plastid.txt --out final.remove_plastid.fasta
        cd ..
        mkdir 08_final
        cat 07_new_assembly/final.200K.fasta  07_new_assembly/final.remove_plastid.fasta > 08_final/HDA.contigs.fasta
  """
```
