include: "rule/00_common.smk"

rule all:
	input:
		###################  HiC  ##################
		expand("cleandata/03_hic/{sample}_clean_{n}.fq.gz", sample = Hsample, n = HAP),
		expand("results/02_assembly/03_hic/01_qc/02_fastqc/{sample}_clean_{n}_fastqc.html", sample = Hsample, n = HAP),
		expand("results/02_assembly/03_hic/01_qc/02_fastqc/{sample}_clean_{n}_fastqc.zip", sample = Hsample, n = HAP),
		expand("results/02_assembly/03_hic/02_ref/{sample}_{enzyme}.txt", sample=Psample,enzyme=enzyme),
		expand("results/02_assembly/03_hic/04_merge/{sample}.merged_nodups.txt",sample=Psample),
		expand("results/02_assembly/03_hic/05_3d-dna/{sample}_hifiasm.hic.p_ctg_clean.0.hic", sample=Psample),


include: "rule/01_clean.smk"
include: "rule/02_hic.smk"
