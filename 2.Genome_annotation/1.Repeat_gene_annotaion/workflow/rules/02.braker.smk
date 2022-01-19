configfile: "config.yaml"

###########################   QC   ##############################
rule fastp_RNA:
    input:
        R1 = "rawdata/04.rna-seq/{sample}.R1.fastq.gz",
        R2 = "rawdata/04.rna-seq/{sample}.R2.fastq.gz"
    output:
        R1 = "cleandata/04.rna-seq/{sample}.clean.R1.fq.gz",
        R2 = "cleandata/04.rna-seq/{sample}.clean.R2.fq.gz",
        json = "results/03.annotation/04.RNA-seq/01.QC/01.fastp/{sample}.QC.fastp.json",
        html = "results/03.annotation/04.RNA-seq/01.QC/01.fastp/{sample}.QC.fastp.html"
    log:
        "logs/annotation/04.rna-seq/01.QC/{sample}.fastp.log"
    benchmark:
        "benchmarks/annotation/04.rna-seq/01.QC/{sample}.fastp.txt"
    params:
        nbase = 5,
        window = 4,
        mean_qual = 20,
        length = 75,
        quality = 15
    threads: 4
    shell:
        """
        {module}
        {QC} || echo "true"
        fastp -a auto --adapter_sequence_r2 auto --detect_adapter_for_pe -w {threads} -i {input.R1} -I {input.R2} --n_base_limit {params.nbase} --cut_window_size {params.window} --cut_mean_quality {params.mean_qual} --length_required {params.length} --qualified_quality_phred {params.quality} -o {output.R1} -O {output.R2} --json {output.json} --html {output.html} >&2 2>{log}
        """



rule md5_RNA:
    input:
        R1 = rules.fastp_RNA.output.R1,
        R2 = rules.fastp_RNA.output.R2
    output:
        md5_1 = "cleandata/04.rna-seq/MD5/{sample}.clean.R1.fq.gz.md5",
        md5_2 = "cleandata/04.rna-seq/MD5/{sample}.clean.R2.fq.gz.md5"
    benchmark:
        "benchmarks/annotation/04.rna-seq/01.QC/{sample}.md5.txt",
    shell:
        """
        md5sum {input.R1} 1> {output.md5_1}
        md5sum {input.R2} 1> {output.md5_2}
        """

############# QC fastqc #############
rule fastqc_RNA:
    input:
        R1 = rules.fastp_RNA.output.R1,
        R2 = rules.fastp_RNA.output.R2
    output:
        html1 = "results/03.annotation/04.RNA-seq/01.QC/02.Fastqc/{sample}.clean.R1_fastqc.html",
        html2 = "results/03.annotation/04.RNA-seq/01.QC/02.Fastqc/{sample}.clean.R2_fastqc.html",
        zip1 = "results/03.annotation/04.RNA-seq/01.QC/02.Fastqc/{sample}.clean.R1_fastqc.zip",
        zip2 = "results/03.annotation/04.RNA-seq/01.QC/02.Fastqc/{sample}.clean.R2_fastqc.zip"
    params:
        outdir = "results/03.annotation/04.RNA-seq/01.QC/02.Fastqc/"
    log:
        "logs/annotation/04.rna-seq/01.QC/{sample}.Fastqc.log"
    benchmark:
        "benchmarks/annotation/04.rna-seq/01.QC/{sample}.Fastqc.txt"
    threads: 4
    shell:
        """
        {module}
        {QC}
        fastqc  -t {threads} {input.R1} {input.R2} -o {params.outdir} >&2 2>{log}
        """



############################   hisat2    ##################
# rule sl:
#     input:
#         fa = "..../genome.fa",
#         fai = "..../genome.fa.fai"
#     output:
#         fa = "results/03.annotation/02.RNA-seq/ref/genome.fa",
#         fai = "results/03.annotation/02.RNA-seq/ref/genome.fa.fai"
#     shell:
#         """
#         ln -sfr {input.fa} {output.fa}
#         ln -sfr {input.fai} {output.fai}
#         """

rule hisat_build:
    input:
        ref = "results/03.annotation/01.ref/{name}.fa"
    output:
        ht = "results/03.annotation/01.ref/{name}.8.ht2"
    params:
        DIR = "results/03.annotation/01.ref"
    benchmark:
        "benchmarks/annotation/04.rna-seq/hisat/{name}.build.txt"
    threads: 12
    shell:
        """
        {mapping}
        cp ctl/.gm_key ~/
        cp -r /public/home/baozhigui/miniconda3/envs/braker2/config/ ./
        hisat2-build -p {threads} {input.ref} {params.DIR}/{wildcards.name}
        """ 

rule hisat_mapping:
    input:
        R1 = rules.fastp_RNA.output.R1,
        R2 = rules.fastp_RNA.output.R2,
        index = rules.hisat_build.output.ht
    output:
        bam = "results/03.annotation/04.RNA-seq/hisat/{name}.{sample}.sorted.bam",
        summary = "results/03.annotation/04.RNA-seq/hisat/{name}.{sample}.mapping.summary"
    params: 
        index = "results/03.annotation/01.ref/{name}"
    log:
        "logs/annotation/04.rna-seq/hisat/hisat.{name}.{sample}.log"
    benchmark:
        "benchmarks/annotation/04.rna-seq/hisat/hisat.{name}.{sample}.txt"
    threads: 12
    shell:
        """
        {mapping}
        hisat2 --dta -x {params.index} -1 {input.R1} -2 {input.R2} --rna-strandness RF --summary-file {output.summary} --new-summary -p {threads} | samtools view -@ 4 -Sb - | samtools sort -@ 6 -o {output.bam} - >&2 2>{log} 
        """ 

rule stringtie:
    input:
        bam = rules.hisat_mapping.output.bam
    output:
        gtf = "results/03.annotation/04.RNA-seq/stringtie/{name}.{sample}.gtf",
        gff3 = "results/03.annotation/04.RNA-seq/stringtie/{name}.{sample}.gff3"
    benchmark:
        "benchmarks/annotation/04.rna-seq/stringtie/{name}.{sample}.txt"
    threads: 4
    shell:
        """
        {mapping}
        stringtie -p {threads} --rf -l {wildcards.sample} -o {output.gtf} {input.bam}
        gffread -E {output.gtf} -o - | sed "s#transcript#match#g" | sed "s#exon#match_part#g" > {output.gff3}
        """ 


#############################   hardmask  genome  ####################

rule hardmask_hisat_build:
    input:
        ref = "results/03.annotation/01.ref/hardmask.{name}.fa"
    output:
        ht = "results/03.annotation/01.ref/hardmask.{name}.8.ht2"
    params:
        DIR = "results/03.annotation/01.ref"
    benchmark:
        "benchmarks/annotation/04.rna-seq/hisat/hardmask.{name}.build.txt"
    threads: 12
    shell:
        """
        {mapping}
        hisat2-build -p {threads} {input.ref} {params.DIR}/hardmask.{wildcards.name}
        """ 

rule hardmask_hisat_mapping:
    input:
        R1 = rules.fastp_RNA.output.R1,
        R2 = rules.fastp_RNA.output.R2,
        index = rules.hardmask_hisat_build.output.ht
    output:
        bam = "results/03.annotation/04.RNA-seq/hisat-hardmask/{name}.{sample}.sorted.bam",
        summary = "results/03.annotation/04.RNA-seq/hisat-hardmask/{name}.{sample}.mapping.summary"
    params: 
        index = "results/03.annotation/01.ref/hardmask.{name}"
    log:
        "logs/annotation/04.rna-seq/hisat-hardmask/hisat.{name}.{sample}.log"
    benchmark:
        "benchmarks/annotation/04.rna-seq/hisat-hardmask/hisat.{name}.{sample}.txt"
    threads: 20
    shell:
        """
        {mapping}
        hisat2 --dta -x {params.index} -1 {input.R1} -2 {input.R2} --rna-strandness RF --summary-file {output.summary} --new-summary -p {threads} | samtools view -@ 4 -Sb - | samtools sort -@ 6 -o {output.bam} - >&2 2>{log} 
        """ 



##########################   ISO-seq   #######################
if Isamples:
    rule hardmask_iso_mapping:
        input:
            ref = "results/03.annotation/01.ref/hardmask.genome.fa",
            fa = "results/03.annotation/03.iso-seq/03.cluster/{sample}/{sample}.polished.hq.fasta.gz"
        output:
            bam = "results/03.annotation/03.iso-seq/04.hardmask_mapping/{sample}.sorted.bam"
        benchmark:
            "benchmarks/annotation/03.iso-seq/hardmask_mapping/minimap2.{sample}.txt"
        threads: 8
        shell:
            """
            {mapping}
            minimap2 -ax splice:hq -uf {input.ref} {input.fa} -t {threads} | samtools view -@ 4 -Shb - | samtools sort -@ 4 -o {output.bam}
            """ 

    ####################  braker   #############
    rule braker:
        input:
            ref = "results/03.annotation/01.ref/softmask.{name}.fa",
            bam = expand("results/03.annotation/04.RNA-seq/hisat-hardmask/{name}.{sample}.sorted.bam",name=name,sample=Rsamples),
            iso_bam = expand("results/03.annotation/03.iso-seq/04.hardmask_mapping/{name}.{sample}.sorted.bam",name=name,sample=Isamples)
        output:
            gff3 = "results/03.annotation/05.braker_{name}/augustus.hints.gff3"
        params:
            DIR = "results/03.annotation/05.braker_{name}",
            bam = ",".join("{input.bam}") + "," + ",".join("{input.iso_bam}")
        benchmark:
            "benchmarks/annotation/05.braker/braker.{name}.txt"
        threads: 20
        run:
            bam = ",".join(input.bam) + "," + ",".join(input.iso_bam)
            shell("{module} && {genemark} && {braker2} && braker.pl --species=braker_{name} --genome={input.ref} --workingdir={params.DIR} --gff3 --AUGUSTUS_CONFIG_PATH={augustus_conf}  --nocleanup --softmasking --bam={bam} --cores {threads}")

else :
    ####################  braker   #############
    rule braker:
        input:
            ref = "results/03.annotation/01.ref/softmask.{name}.fa",
            bam = expand("results/03.annotation/04.RNA-seq/hisat-hardmask/{name}.{sample}.sorted.bam",name=name,sample=Rsamples)
        output:
            gff3 = "results/03.annotation/05.braker_{name}/augustus.hints.gff3"
        params:
            DIR = "results/03.annotation/05.braker_{name}"
        benchmark:
            "benchmarks/annotation/05.braker/braker.{name}.txt"
        threads: 20
        run:
            bam = ",".join(input.bam)
            shell("{module} && {genemark} && {braker2} && braker.pl --species=braker_{name} --genome={input.ref} --workingdir={params.DIR} --gff3 --AUGUSTUS_CONFIG_PATH={augustus_conf}  --nocleanup --softmasking --bam={bam} --cores {threads}")
