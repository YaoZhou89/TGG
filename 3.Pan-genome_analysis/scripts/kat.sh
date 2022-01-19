#/bin/bash
out=Heinz_1706
threads=96
options="-H 10000000000 -I 10000000000 -m 31 -h"
kat comp -o $out -t $threads $options <(gunzip -c Heinz_1706.fq.gz) Heinz_1706.fasta