library(sveval)
library(VariantAnnotation)
library(GenomicRanges)
library(dplyr)
library(ggplot2)

## full catalog with variants ids
type.df = lapply( c(1:12), function(chrn){
  message(chrn)
  svs = readSVvcf(paste0('../02_redup/graph-fordedup-', chrn, '.vcf'),out.fmt=c("vcf"), keep.ids=TRUE)
  tibble(id=names(svs), type=info(svs)$SVTYPE)
})
type.df = do.call(rbind, type.df)

## cluster info
cl.df = read.table('../02_redup/neardups-clusters.tsv', as.is=TRUE, header=TRUE)
ids.l = strsplit(cl.df$ids, ',')
cl.df = lapply(1:length(ids.l), function(ii) tibble(cmp=cl.df$cmp[ii], id=ids.l[[ii]]))
cl.df = do.call(rbind, cl.df)

## vg calls
vcf = readVcf('called_vcf/called_merged.vcf')
calls.df = tibble(id=rownames(geno(vcf)$GT), gt=geno(vcf)$GT[,1], ad=sapply(geno(vcf)$AD, '[', 2), dp=info(vcf)$DP)

calls.cl = merge(calls.df, cl.df) %>% merge(type.df) %>% arrange(cmp)

write.table(calls.cl, file='called_vcf/graph-fordedup-support.tsv', sep='\t', quote=FALSE, row.names=FALSE)

## explore results
calls.cl.s = calls.cl %>% group_by(cmp, type) %>%
  mutate(ad.ratio=ad/dp, ad.ratio.cmp=ad/max(ad)) %>%
  summarize(nvar=n(), 
            top.ad.ratio=sort(ad.ratio, decreasing=TRUE)[1],
            top2.ad.ratio=sort(ad.ratio, decreasing=TRUE)[2],
            top.ad.ratio.cmp=sort(ad.ratio.cmp, decreasing=TRUE)[1],
            top2.ad.ratio.cmp=sort(ad.ratio.cmp, decreasing=TRUE)[2])
                                
pdf('graph-fordedup-support.pdf', 9, 7)
ggplot(calls.cl.s, aes(x=top2.ad.ratio, y=top.ad.ratio, size=nvar, colour=type)) +
  geom_point(alpha=.3) + theme_bw() + xlim(0,1) + ylim(0,1)
ggplot(calls.cl.s, aes(x=top.ad.ratio-top2.ad.ratio, fill=factor(nvar))) +
  geom_histogram() + theme_bw()
ggplot(calls.cl.s, aes(x=top.ad.ratio-top2.ad.ratio)) +
  geom_histogram() + theme_bw() + facet_grid(nvar~., scales='free')
ggplot(calls.cl.s, aes(x=top2.ad.ratio.cmp, fill=factor(nvar))) +
  geom_histogram() + theme_bw()
ggplot(calls.cl.s, aes(x=top2.ad.ratio.cmp)) +
  geom_histogram() + theme_bw() + facet_grid(nvar~., scales='free')
dev.off()

write.table(calls.cl.s, file='graph-fordedup-support-summary.tsv', row.names=FALSE, sep='\t')
