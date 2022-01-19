configfile: "config.yaml"


#############################   QC   ########################
rule fastp_HiC:
	input:
		R1 = "rawdata/03_hic/{sample}_1.fq.gz",
		R2 = "rawdata/03_hic/{sample}_2.fq.gz"
	output:
		R1 = "cleandata/03_hic/{sample}_clean_1.fq.gz",
		R2 = "cleandata/03_hic/{sample}_clean_2.fq.gz",
		json = "logs/03_hic/01_qc/{sample}_fastp.json",
		html = "logs/03_hic/01_qc/{sample}_fastp.html"
	log:
		"logs/03_hic/01_qc/{sample}_fastp.log"
	threads:
		16
	shell:
		"""
		singularity exec GenomeAssemblyContainer_v0.2 \
			fastp --detect_adapter_for_pe \
			-w {threads} \
			-i {input.R1} \
			-I {input.R2} \
			-o {output.R1} \
			-O {output.R2} \
			--json {output.json} \
			--html {output.html} >&2 2>{log}
		"""


##########################    QC fastqc     ############################
rule fastqc_Hic:
	input:
		R1 = rules.fastp_HiC.output.R1,
		R2 = rules.fastp_HiC.output.R2
	output:
		html1 = "results/02_assembly/03_hic/01_qc/02_fastqc/{sample}_clean_1_fastqc.html",
		html2 = "results/02_assembly/03_hic/01_qc/02_fastqc/{sample}_clean_2_fastqc.html",
		zip1 = "results/02_assembly/03_hic/01_qc/02_fastqc/{sample}_clean_1_fastqc.zip",
		zip2 = "results/02_assembly/03_hic/01_qc/02_fastqc/{sample}_clean_2_fastqc.zip"
	params:
		outdir = "results/02_assembly/03_hic/01_qc/02_fastqc"
	log:
		"logs/03_hic/01_qc/{sample}_fastqc.log"
	threads:
		64
	shell:
		"""
		singularity exec GenomeAssemblyContainer_v0.2 \
			fastqc -t {threads} {input.R1} {input.R2} \
			-o {params.outdir} >&2 2>{log}
		"""


#########################    data split    ###########################
rule fastq_hic_split:
	input:
		R1 = rules.fastp_HiC.output.R1,
		R2 = rules.fastp_HiC.output.R2
	output:
		split_R1 = expand("results/02_assembly/03_hic/01_qc/03_fastq_split/{{sample}}_clean_1.part_{part}.fq.gz",part=LIST),
		split_R2 = expand("results/02_assembly/03_hic/01_qc/03_fastq_split/{{sample}}_clean_2.part_{part}.fq.gz",part=LIST)
	threads:
		64
	params:
		part = {part},
		dir = "results/02_assembly/03_hic/01_qc/03_fastq_split"
	shell:
		"""
		singularity exec GenomeAssemblyContainer_v0.2 \
			seqkit split2 -f -j {threads} \
			-p {params.part} \
			-O {params.dir} \
			-1 {input.R1} \
			-2 {input.R2}
		"""



#######################     contig index and enzyme site position    #######################
rule bwa_index:
	input:
		ref = "refenrence/"+ref
	output:
		ref_touch = "results/02_assembly/03_hic/02_ref/{sample}_hifiasm_index_ok"
	params:
		ok = "{sample}_hifiasm_index_ok"
	shell:
		"""
		singularity exec GenomeAssemblyContainer_v0.2 \
			bwa index {input.ref} && touch results/02_assembly/03_hic/02_ref/{params.ok}
		"""

rule get_site:
	input:
		ref = "refenrence/"+ref
	output:
		site = expand("results/02_assembly/03_hic/02_ref/{{sample}}_{enzyme}.txt",enzyme = enzyme)
	shell:
		"""
		singularity exec GenomeAssemblyContainer_v0.2 \
			python /home/software/juicer/juicer/misc/generate_site_positions.py \
			{enzyme} \
			{wildcards.sample} \
			{input.ref}
		mv {wildcards.sample}_{enzyme}.txt results/02_assembly/03_hic/02_ref
		"""



########################    mapping    ########################
rule bwa_mem:
	input:
		ref = "refenrence/"+ref,
		ok = rules.bwa_index.output.ref_touch,
		R1 = "results/02_assembly/03_hic/01_qc/03_fastq_split/{sample}_clean_1.part_{part}.fq.gz",
		R2 = "results/02_assembly/03_hic/01_qc/03_fastq_split/{sample}_clean_2.part_{part}.fq.gz"
	output:
		sam=temp("results/02_assembly/03_hic/03_mapping/{sample}_clean_part_{part}.sam")
	threads:
		4
	shell:
		"""
		singularity exec GenomeAssemblyContainer_v0.2 \
			bwa mem -SP5M -t {threads} {input.ref} {input.R1} {input.R2} > {output.sam}
		"""


########################## chimeric reads dealing ############################
rule chimeric_blacklist:
	input:
		sam = rules.bwa_mem.output.sam
	output:
		abnorm_sam=temp("results/02_assembly/03_hic/04_merge/{sample}_part_{part}_abnorm.sam"),
		unmapped_sam=temp("results/02_assembly/03_hic/04_merge/{sample}_part_{part}_unmapped.sam"),
		norm_txt=temp("results/02_assembly/03_hic/04_merge/{sample}_part_{part}_norm.txt")
	shell:
		"""
		SINGULARITY_LC_ALL=C singularity exec GenomeAssemblyContainer_v0.2 \
			awk -v "fname1"={output.norm_txt} -v "fname2"={output.abnorm_sam} -v "fname3"={output.unmapped_sam} -f /home/software/juicer/juicer/CPU/common/chimeric_blacklist.awk {input.sam}
		"""


########################## get enzyme site region reads ############################
rule frag:
	input:
		norm_txt=rules.chimeric_blacklist.output.norm_txt,
		site=rules.get_site.output.site
	output:
		frag_txt=temp("results/02_assembly/03_hic/04_merge/{sample}_part_{part}_frag.txt")
	shell:
		"""
		singularity exec GenomeAssemblyContainer_v0.2 \
			/home/software/juicer/juicer/CPU/common/fragment.pl {input.norm_txt} {output.frag_txt} {input.site}
		"""

rule sort_frag:
	input:
		frag_txt=rules.frag.output.frag_txt
	output:
		sort_txt=temp("results/02_assembly/03_hic/04_merge/{sample}_part_{part}_sort.txt")
	threads:
		4
	shell:
		"""
		mkdir -p HIC_tmp
		SINGULARITY_LC_ALL=C singularity exec GenomeAssemblyContainer_v0.2 \
			sort --parallel={threads} -S 5G -T HIC_tmp -k2,2d -k6,6d -k4,4n -k8,8n -k1,1n -k5,5n -k3,3n {input.frag_txt} > {output.sort_txt}
		"""



########################## merge and remove duplications ###############################
rule merged_nodups:
	input:
		txt=expand("results/02_assembly/03_hic/04_merge/{{sample}}_part_{part}_sort.txt",part=LIST)
	output:
		merged_sort = temp("results/02_assembly/03_hic/04_merge/{sample}_merged_sort.txt"),
		merged_nodups = "results/02_assembly/03_hic/04_merge/{sample}.merged_nodups.txt"
	params:
		prex="results/02_assembly/03_hic/04_merge/{sample}."
	threads:
		64
	shell:
		"""
		SINGULARITY_LC_ALL=C singularity exec GenomeAssemblyContainer_v0.2 \
			sort -S 5G --parallel={threads} -T HIC_tmp -m -k2,2d -k6,6d -k4,4n -k8,8n -k1,1n -k5,5n -k3,3n {input.txt}  > {output.merged_sort}
		SINGULARITY_LC_ALL=C singularity exec GenomeAssemblyContainer_v0.2 \
			touch {output.merged_nodups} && SINGULARITY_LC_ALL=C singularity exec GenomeAssemblyContainer_v0.2 \
			awk -f /home/software/juicer/juicer/CPU/common/dups.awk -v name={params.prex} -v nowobble=1 {output.merged_sort}
		"""

########################## scaffolding ##############################
rule scaffolding:
	input:
		ref = "refenrence/"+ref,
		merged_nodups = rules.merged_nodups.output.merged_nodups
	output:
		fa ="results/02_assembly/03_hic/05_3d-dna/{sample}_hifiasm.hic.p_ctg_clean.0.hic"
	params:
		ref = "../../../../refenrence/"+ref,
		dir = "results/02_assembly/03_hic/05_3d-dna",
		merged_nodups = "../04_merge/{sample}.merged_nodups.txt",
		round_num = 0,
		q = 1
	threads:
		64
	shell:
		"""
		mkdir -p {params.dir}
		cd {params.dir}
		ln -s ../../../../GenomeAssemblyContainer_v0.2 .
		SINGULARITY_LC_ALL=C singularity exec GenomeAssemblyContainer_v0.2 \
			bash /home/software/3d-dna/run-asm-pipeline.sh -q {params.q} -r {params.round_num} --editor-repeat-coverage 50 -e {params.ref} {params.merged_nodups}
		"""

