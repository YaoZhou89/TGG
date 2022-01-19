## orignal code: https://github.com/vgteam/giraffe-sv-paper/blob/master/scripts/sv/merge-svs/filter.R

library(VariantAnnotation)
vcf = readVcf('/public10/home/sci0011/projects/tomato/100tomato/merged.ont.v1.0.vcf.gz')
## PASS filter
table(rowRanges(vcf)$FILTER)
vcf = vcf[which(rowRanges(vcf)$FILTER == 'PASS')]
seqlevels(vcf) = paste0('chr', seqlevels(vcf))
writeVcf(vcf, '100tomato.filtered.vcf')
