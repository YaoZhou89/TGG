#!/bin/bash
for rep in {1..20}
do
for fold in {1..5}
do
Rscript indel.r $rep $fold &
done
wait
done
wait
