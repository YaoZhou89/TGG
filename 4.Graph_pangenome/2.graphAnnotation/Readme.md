#  Tomato Graph Annotation (TGA1.0)
Mapping gene to graph and get unmapped gene from the gmap result to get non-redundant gene list for tomato Pangenome. The corresponding scripts are deposited in directory.



Get each accesion gene

```shell
for i in `cat genome.list`;awk '$3 == "mRNA"' ${i}.longest.gff3|cut -f1 -d ";"|awk '{print $1"\t"$4-1"\t"$5"\t"$9"\t.\t"$7}'|sed "s#ID=##g" > ${i}_longest_mRNA.bed
for i in `cat genome.list`;do bedtools getfasta -fi ${i}.fasta -fo ${i}_longest_mRNA.cdna.fasta -bed ${i}_longest_mRNA.bed -nameOnly -s;done
```

Using `minigraph` to  map
```shell
for i in `cat genome.list`;do ~/software/minigraph/minigraph -x lr -t 64 ./G46.gfa ./00.mRNA/mRNA_fasta/${i}_longest_mRNA.fasta > 01.gaf/${i}_mRNA_xlr.gaf
for i in `cat ../00.mRNA/genome.list`;do grep -P ">|<" ../01.gaf/${i}_mRNA_xlr.gaf >> graph_gene.gaf;done
```

Remap gene to SL5.0 ref to filter false mapping
```shell
cut -f1 graph_gene.gaf |sort|uniq > graph_path_gene.list
seqkit grep -f graph_path_gene.list all.cds.fa > graph_path_gene.cds.fa
gmap -B 5 -A -O -n 1 --no-chimeras -D ./ -d SL5 -t 12 -f gff3_gene --gff3-swap-phase 1 graph_path_gene.cds.fa  > graph_path_gene.gff3
perl scripts/filter_gmap.pl graph_path_gene.gff3 90 90 > filter.gff3
awk '$3 == "mRNA"' filter.gff3|cut -f9|cut -f1 -d ";"|sort|uniq|sed "s#.mrna1##g" > graph_path_gene.filter.list
```

Use gmap to SL5.0 ref
```shell
for i in `cat ../00.mRNA/genome.list`;do perl scripts/filter_gmap.pl gff3/${i}_SL5.gff3 90 90 > ${i}_filter.gff3;done
for i in `cat ../00.mRNA/genome.list`;do awk '$3 == "mRNA"' filter/${i}_filter.gff3|cut -f1 -d ";"|awk '{print $1"\t"$4-1"\t"$5"\t"$9"\t.\t"$7}'|sed "s#ID=##g" > ${i}_filter.bed;done
for i in `cat ../00.mRNA/genome.list`;do bedtools intersect -b SL5_longest_mRNA.bed -a ./03.bed/${i}_filter.bed -wao> ${i}_SL5_intersect.bed;done
for i in `cat ../00.mRNA/genome.list`;do awk '$8 == -1' 04.intersect/${i}_SL5_intersect.bed > 04.intersect/remain/$i.bed;done
cat *.bed|sort -k1,1n -k2,2n > cluster/sort.bed

perl scripts/cluster.pl sort.bed |cut -f1-6 >  sort_filter.bed
bedtools merge -d 1000 -i sort_filter.bed > cluster.bed
bedtools intersect -a cluster.bed -b sort_filter.bed -wao > cluster_sort_filter.bed
perl scripts/select_gene.pl cluster_sort_filter.bed order.txt > cluster_select.bed
cut -f7 cluster_select.bed|sed "s#.mrna1##g" > cluster_select.gene.list
```
Combine and filter the result
```shell
cat graph_path_gene.filter.list cluster_select.gene.list > graph_gene.list
seqkit grep -f graph_gene.list all.pep.fa > graph_gene.pep.fa
cd-hit -i ./graph_gene.pep.fa -o graph_gene_c90.pep.fa -c 0.9 -n 5 -M 16000 -d 0 -T 64
perl scripts/cdhit_select.pl graph_gene_c90.pep.fasta.clstr order.txt > graph_gene_c90.select.fasta
```