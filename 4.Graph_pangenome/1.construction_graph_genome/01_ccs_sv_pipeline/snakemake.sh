
source /public/agis/huangsanwen_group/chenglin/softwares/miniconda3/bin/activate base


snakemake -s Snakefile -j 2 --stats snakejob.stats >&2 2>>snakejob.log
