# Scaffold to chromosome-level assembly



## Tomato Heinz 1706 with Hi-C reads 

This directory includes the workflow and files refer to scaffold. An analysis pipeline in `snakemake` to streamline the processing for contig assembly scaffolding.

The pipeline contain following steps:

- Quality control for Hi-C reads by `fastp` 
- Identify the enzyme site position, merge and remove duplication by `juicer`
- `HIC-Pro` and `JABT`  extract valid pairs 
- Scaffolding the contigs by `3d-DNA`

```shell
#!/bin/bash
#SBATCH --exclusive
#SBATCH -p amd_256
#SBATCH -N 1
#SBATCH -n 64
snakemake -s scaffold_Snakefile --stats snake.stats --latency-wait 120 -k -j 32
```

## Workflow

![image-20220103183236565](https://github.com/YaoZhou89/TGG/blob/main/figs/1_2genome_assembly_scaffold.png)


Figure1.  Pipeline of genome scaffold.



## Tomato accessions without Hi-C reads

The accessions belonging to BIG or CER groups are directly guided using the Heinz 1706 assembly, while the remaining PIM groups are guided using LA2093 assembly using `Ragtag`.

```shell
conda install -c bioconda ragtag
# scaffold a Heinz_1706 assembly
ragtag.py scaffold Heinz_1706.fasta query.fasta
# scaffold a Heinz_1706 assembly
ragtag.py scaffold LA2093.fasta query.fasta
```





