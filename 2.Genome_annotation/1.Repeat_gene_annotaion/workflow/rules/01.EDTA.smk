
configfile: "config.yaml"

############################   EDTA    ##################
rule sl:
    input:
        fa = "{name}.fa",
    output:
        in_fai = "{name}.fa" + ".fai",
        fa = "results/03.annotation/02.repeat/{name}.fa",
        fai = "results/03.annotation/02.repeat/{name}.fa.fai",
        ref = "results/03.annotation/01.ref/{name}.fa",
        ref_fai = "results/03.annotation/01.ref/{name}.fa.fai"
    params:
        fa = outdir + "/{name}.fa", 
        fai = outdir + "/{name}.fa.fai", 
    shell:
        """
        {module}
        {asm}
        samtools faidx {input.fa}
        ln -sf {params.fa} {output.fa}
        ln -sf {params.fa} {output.ref}
        ln -sf {params.fai} {output.fai}
        ln -sf {params.fai} {output.ref_fai}
        """


rule LTR:
    input:
        fa = rules.sl.output.fa,
        fai = rules.sl.output.fai
    output:
        fa = "results/03.annotation/02.repeat/{name}.fa.mod.EDTA.raw/{name}.fa.mod.LTR.raw.fa"
    params:
        DIR = "results/03.annotation/02.repeat"
    benchmark: 
        "benchmarks/annotation/02.repeat/{name}.LTR.txt" 
    threads: 12
    shell:
        """
        {module}
        {EDTA}
        cd {params.DIR} && EDTA_raw.pl --genome {wildcards.name}.fa --species {species} --type ltr --threads {threads}
        """

rule TIR:
    input:
        fa = rules.sl.output.fa,
        fai = rules.sl.output.fai
    output:
        fa = "results/03.annotation/02.repeat/{name}.fa.mod.EDTA.raw/{name}.fa.mod.TIR.raw.fa"
    params:
        DIR = "results/03.annotation/02.repeat"
    benchmark: 
        "benchmarks/annotation/02.repeat/{name}.TIR.txt" 
    threads: 12
    shell:
        """
        {module}
        {EDTA}
        cd {params.DIR} 
        EDTA_raw.pl --genome {wildcards.name}.fa --species {species} --type tir --threads {threads}
        """


rule Helitron:
    input:
        fa = rules.sl.output.fa,
        fai = rules.sl.output.fai
    output:
        fa = "results/03.annotation/02.repeat/{name}.fa.mod.EDTA.raw/{name}.fa.mod.Helitron.raw.fa"
    params:
        DIR = "results/03.annotation/02.repeat"
    benchmark: 
        "benchmarks/annotation/02.repeat/{name}.Helitron.txt" 
    threads: 12
    shell:
        """
        {module}
        {EDTA}
        cd {params.DIR} && EDTA_raw.pl --genome {wildcards.name}.fa --species {species} --type helitron --threads {threads}
        """


rule final:
    input:
        LTR = rules.LTR.output.fa,
        TIR = rules.TIR.output.fa,
        Helitron = rules.Helitron.output.fa
    output:
        gff = "results/03.annotation/02.repeat/{name}.fa.mod.EDTA.anno/{name}.fa.mod.EDTA.TEanno.gff3",
        lib = "results/03.annotation/02.repeat/{name}.fa.mod.EDTA.anno/{name}.fa.mod.EDTA.TElib.fa"
    params:
        DIR = "results/03.annotation/02.repeat"
        #problib = "/EDTA/20191204/database/alluniRefprexp082813"
        #blast = "/data/software/blast+/ncbi-blast-v2.7.1+/bin/",
        #masker = "/data/software/RepeatMasker/RepeatMasker/",
        #modeler = "/data/software/RepeatModeler/1.0.11/"
    benchmark: 
        "benchmarks/annotation/02.repeat/{name}.EDTA.final.txt" 
    threads: 12
    run:
        #cmd = EDTA + " && " + RepeatMasker + " && " + RepeatModeler + " && " + blast + " && "
        cmd = module + "&&" + EDTA + " && " 
        cmd = cmd + " cd " + params.DIR + "  && EDTA.pl --genome " + wildcards.name + ".fa  --species " + species
        if lib :
           cmd = cmd + " --curatedlib " + lib
        if cds :
           cmd = cmd + " --cds " + cds
        cmd = cmd + " --sensitive " + str(sensitive) + " --anno 1  --step " + step + " -t " + str(threads) 
        #cmd = cmd + " --blast " + params.blast + " --repeatmasker " + params.masker + " --repeatmodeler " + params.modeler + " -protlib " + params.problib  
        shell("{cmd}")

rule mask:
    input:
        fa = rules.sl.output.fa,
        gff = rules.final.output.gff
    output:
        hard = "results/03.annotation/01.ref/hardmask.{name}.fa",
        soft = "results/03.annotation/01.ref/softmask.{name}.fa"
    benchmark: 
        "benchmarks/annotation/02.repeat/{name}.mask.txt"
    threads: 4 
    shell:
        """
        {module}
        {EDTA}
        grep -v -P "Satellite|rich|Simple_repeat|Low_complexity|tRNA|rRNA" {input.gff} | bedtools maskfasta -fi {input.fa} -bed - -fo {output.soft} -soft
        grep -v -P "Satellite|rich|Simple_repeat|Low_complexity|tRNA|rRNA" {input.gff} | bedtools maskfasta -fi {input.fa} -bed - -fo {output.hard} 
        """
