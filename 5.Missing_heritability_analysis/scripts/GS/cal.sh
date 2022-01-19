#!/bin/bash
for fold in {1..5}
do
for rep in {1..20}
do
  Rscript rrblup.r $rep $fold &
done
wait
done
wait
