configfile: "config.yaml"

#### maker round1
rule maker_round1:
    input:
        gff3 = expand("results/03.annotation/04.RNA-seq/stringtie/{name}.{sample}.gff3",name=name,sample=Rsamples),
        pro = "rawdata/05.Annotation/proteins.fasta",
        lib = "results/03.annotation/02.repeat/{name}.fa.mod.EDTA.anno/{name}.fa.mod.EDTA.TElib.fa"
    output:
        maker_gff = "results/03.annotation/06.maker/round1/{name}.maker.output/{name}.maker.gff"
    log:
        "{name}.maker_round1.log"
    benchmark:
        "benchmarks/annotation/06.maker/{name}.maker_round1.benchmark.txt"
    threads: 120
    params:
        fasta = outdir + "/{name}.fa",
        pro = outdir + "/" + "rawdata/05.Annotation/proteins.fasta",
        te = outdir + "/" + "results/03.annotation/02.repeat/{name}.fa.mod.EDTA.anno/{name}.fa.mod.EDTA.TElib.fa"
    shell:
        """
        {module}
        {maker}
        mkdir -p results/03.annotation/06.maker/round1 && cd results/03.annotation/06.maker/round1
        cp ../../../../ctl/.gm_key ~/
        cp ../../../../ctl/maker/round1/*.ctl ./
        num=`ls {outdir}/results/03.annotation/04.RNA-seq/stringtie/*.gff3|wc -l`
        if [ $num -gt 1 ];then
            gff=`ls {outdir}/results/03.annotation/04.RNA-seq/stringtie/*.gff3|tr '\\n' ','|sed "s/,$//g"`
        else
            gff=`ls {outdir}/results/03.annotation/04.RNA-seq/stringtie/*.gff3`
        fi
        echo $gff

        sed -i 's#^genome=#genome={params.fasta}#g' maker_opts.ctl
        sed -i -e "s|^est_gff=|est_gff=$gff|g" maker_opts.ctl
        sed -i 's#^protein=#protein={params.pro}#g' maker_opts.ctl
        sed -i 's#^rmlib=#rmlib={params.te}#g' maker_opts.ctl
        mpiexec -n {threads} maker -base {wildcards.name} maker_bopts.ctl maker_exe.ctl maker_opts.ctl  >&2 2>{log}
        cd {wildcards.name}.maker.output
        gff3_merge -d ./{wildcards.name}_master_datastore_index.log -o {wildcards.name}.maker.gff 
        fasta_merge -d ./{wildcards.name}_master_datastore_index.log -o {wildcards.name}
        """

rule SNAP_Training:
    input:
        rules.maker_round1.output.maker_gff
    output:
        hmm = "results/03.annotation/06.maker/round2/Training/SNAP/{name}.hmm"
    params:
        aed = "0.1",
        hmm = "{name}.hmm"
    shell:
        """
        mkdir -p results/03.annotation/06.maker/round2/Training/SNAP
        cd results/03.annotation/06.maker/round2/Training/SNAP
        {module}
        {maker} 
        maker2zff -x {params.aed} ../../../round1/{wildcards.name}.maker.output/{wildcards.name}.maker.gff
        # gather some stats and validate
        fathom genome.ann genome.dna -gene-stats > gene-stats.log 2>&1
        fathom genome.ann genome.dna -validate > validate.log 2>&1
        # collect the training sequences and annotations, plus 1000 surrounding bp for training
        fathom genome.ann genome.dna -categorize 1000 > categorize.log 2>&1
        fathom uni.ann uni.dna -export 1000 -plus > uni-plus.log 2>&1
        mkdir params
        cd params
        forge ../export.ann ../export.dna > ../forge.log 2>&1
        cd ..
        hmm-assembler.pl genome params > {params.hmm}
        """



rule maker_round2:
    input:
        gff = rules.maker_round1.output.maker_gff,
        hmm = rules.SNAP_Training.output.hmm
    output:
        maker_gff = "results/03.annotation/06.maker/round2/{name}.maker.output/{name}.maker.gff"
    #log:
    #    "logs/annotation/06.maker/{name}.maker_round2.log"
    benchmark:
        "benchmarks/annotation/06.maker/{name}.maker_round2.benchmark.txt"
    threads: 20
    params:
        fasta = outdir + "/{name}.fa",
        config = outdir + "/config",
        snap = outdir + "/" + rules.SNAP_Training.output.hmm,
        gm = outdir + "/" + "results/03.annotation/05.braker_{name}/GeneMark-ET/output/gmhmm.mod", 
        gff = outdir + "/" + rules.maker_round1.output.maker_gff
    shell:
        """
        {module}
        {maker}
        cp -r results/03.annotation/05.braker_{wildcards.name}/species/* config/species/
        export AUGUSTUS_CONFIG_PATH={params.config}
        mkdir -p results/03.annotation/06.maker/round2 && cd results/03.annotation/06.maker/round2
        cp ../../../../ctl/maker/round2/*.ctl ./
        sed -i 's#^genome=#genome={params.fasta}#g' maker_opts.ctl
        sed -i 's#^maker_gff=#maker_gff={params.gff}#g' maker_opts.ctl
        sed -i "s#^snaphmm=#snaphmm={params.snap}#g" maker_opts.ctl
        sed -i "s#^gmhmm=#gmhmm={params.gm}#g" maker_opts.ctl
        sed -i "s#^augustus_species=#augustus_species=braker_{wildcards.name}#g" maker_opts.ctl
        mpiexec -n {threads} maker -base {wildcards.name} maker_bopts.ctl maker_exe.ctl maker_opts.ctl  >&2 2> maker_round2.log
        gff3_merge -d ./{wildcards.name}.maker.output/{wildcards.name}_master_datastore_index.log -o {wildcards.name}.maker.all.gff
        awk '$2 == "maker"' {wildcards.name}.maker.all.gff > {wildcards.name}.maker.only.gff
        fasta_merge -d ./{wildcards.name}.maker.output/{wildcards.name}_master_datastore_index.log -o {wildcards.name}
        """
