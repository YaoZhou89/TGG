

rule seqkit_stats:
	input:
		ref_fa = "refenrence/"+ref
	output:
		ref_fa_stats = "results/03_assembly_assement/01_seqkit_stats/{sample}_hifiasm.hic.p_ctg_clean"
	shell:
		"""
		singularity exec GenomeAssemblyContainer_v0.2 \
		seqkit stats -a -j 64 -t dna -o {output.ref_fa_stats} {input.ref_fa}
		"""


rule stats_n10_n90:
	input:
		ref_fa = "refenrence/"+ref
	output:
		ref_fa_stats = "results/03_assembly_assement/01_seqkit_stats/{sample}_hifiasm.hic.p_ctg_clean_n10_n90"
	shell:
		"""
		singularity exec GenomeAssemblyContainer_v0.2 \
		/home/software/seq_n50.pl {input.ref_fa} > {output.ref_fa_stats}
		"""


rule busco:
	input:
		ref_fa = "refenrence/"+ref
	output:
		ref_short_summary = "results/03_assembly_assement/02_busco/{sample}_hifiasm.hic.p_ctg_clean/short_summary.specific.solanales_odb10.{sample}_hifiasm.hic.p_ctg_clean.txt"
	shell:
		"""
		cd results/03_assembly_assement/02_busco

		singularity exec ../../../busco busco -m genome \
					-i ../../../{input.ref_fa} \
					-o {Psample}_hifiasm.hic.p_ctg_clean \
					-l ../../../../solanales_odb10 \
					-c 64 -f \
					--offline

		cd ../../../
		"""


