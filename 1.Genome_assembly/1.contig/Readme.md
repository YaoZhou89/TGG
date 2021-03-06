
# Diploid genome assembly
This directory includes the workflow and files refer to homozygous species assemble. An analysis pipeline in `contig_Snakefile` to streamline the processing of assembly. 

```shell
snakemake -s contig_Snakefile -j 64 --stats snakejob.stats >&2 2>>snakejob.log
```

The following software are applied:

`Flye` , `Hicanu` and `Hifiasm` are used to assemble primary genome.
`GALA` is used to fillter the potential miss-assembly regions
`WGS` is manual software intergrated in pipeline (https://github.com/YaoZhou89/WGSc) 

## Workflow

![1_2genome_assembly_scaffold](https://github.com/YaoZhou89/TGG/blob/main/figs/1_2genome_assembly_scaffold.png)

Figure1.  Pipeline of genome assembly.





