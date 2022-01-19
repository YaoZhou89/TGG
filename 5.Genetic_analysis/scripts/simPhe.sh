#!/bin/bash
for h2 in 0.15 0.3 0.45 0.6 0.75 0.9
do
for qtn in 10 100 500 1000
do
ldak --make-phenos h2_${h2}-qtn_${qtn} --bfile exp --weights /public10/home/sci0011/projects/tomato/pop/01_snps/01_exp/10_kinship_indel/snps/weights.all --power -.5 --her $h2 --num-phenos 100 --num-causals $qtn &
done
done