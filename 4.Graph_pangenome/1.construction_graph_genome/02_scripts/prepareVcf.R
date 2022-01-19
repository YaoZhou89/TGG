library(sveval)
library(VariantAnnotation)
library(GenomicRanges)
library(dplyr)
library(igraph)

## HSVLR catalog
svs = lapply(c(0:12), function(chrn){
  message(chrn)
  readSVvcf(paste0('../02_redup/graph-fordedup-', chrn, '.vcf'),keep.ins.seq=T,keep.ref.seq=T,out.fmt=c("vcf") , keep.ids=TRUE)
})
svs = do.call(rbind, svs)

## Remove near-duplicates that were less supported by short-reads
calls.cl = read.table('../03_dedup/called_vcf/graph-fordedup-support.tsv',as.is=TRUE, header=TRUE)
v.torm = calls.cl %>% group_by(cmp) %>%
  mutate(ad.ratio=ad/dp) %>% arrange(desc(ad.ratio)) %>%
  do({.[-1,]})

length(svs) ## 179115
svs = svs[which(!(names(svs) %in% v.torm$id))]
length(svs) ## 123785

## rename to include size information
names(svs) = paste0('sv', 1:length(svs), '_', info(svs)$SVLEN)

## Make GRanges object
svs.gr = rowRanges(svs)
svs.gr$size = info(svs)$SVLEN
svs.gr$type = info(svs)$SVTYPE
end(svs.gr) = info(svs)$END


##
#### Separate small insertions for augmentation
##

## Less stringent cluster insertions
ins.gr = svs.gr[which(svs.gr$type=='INS' & svs.gr$size<1000)]
ol.ins = findOverlaps(ins.gr, ins.gr, maxgap=100) %>%
  as.data.frame %>% filter(queryHits<subjectHits) %>% 
  mutate(qs=ins.gr$size[queryHits], ss=ins.gr$size[subjectHits],
         qid=names(ins.gr)[queryHits], sid=names(ins.gr)[subjectHits],
         rol=ifelse(qs>ss, ss/qs, qs/ss)) %>%
  select(-queryHits, -subjectHits) %>% 
  filter(rol > .3)

## Overlap deletions together
del.gr = svs.gr[which(svs.gr$type=='DEL' & svs.gr$size<1000)]
## Cluster deletions
ol.del = findOverlaps(del.gr, del.gr) %>%
  as.data.frame %>% filter(queryHits<subjectHits) %>% 
  mutate(qs=del.gr$size[queryHits], ss=del.gr$size[subjectHits],
         qss=width(pintersect(del.gr[queryHits], del.gr[subjectHits])),
         qid=names(del.gr)[queryHits], sid=names(del.gr)[subjectHits],
         rol=ifelse(qs>ss, qss/qs, qss/ss)) %>%
  select(-queryHits, -subjectHits, -qss) %>% 
  filter(rol > .3)

ol.df = rbind(ol.ins, ol.del)


## Find connected components
ol.g = ol.df %>% select(qid, sid) %>% as.matrix %>% graph_from_edgelist(directed=FALSE)
ol.c = components(ol.g)
ol.c.df = tibble(id=names(ol.c$membership), cmp=as.numeric(ol.c$membership))

## Compute the average reciprocal overlap between each variant with the other variants in the component
ol.s = rbind(ol.df %>% mutate(id=qid) %>% select(-qid, -sid),
             ol.df %>% mutate(id=sid) %>% select(-qid, -sid)) %>%
  group_by(id) %>% summarize(nol=n(), mean.rol=mean(rol))

## Find near-duplicate variants to remove
## All variants except the variant with highest average reciprocal overlap in each component
ol.cs = merge(ol.s, ol.c.df)
v.toaug = ol.cs %>% group_by(cmp) %>%
  arrange(desc(nol), desc(mean.rol)) %>%
  do({.[-1,]})


##
#### Output files
##

svs.for.cons = svs[which(!(names(svs) %in% v.toaug$id))]
length(svs.for.cons) ## 82736
svs.for.aug = svs[which(names(svs) %in% v.toaug$id)]
length(svs.for.aug) ## 41049

## Prepare chunk regions: BED file with chunk boundaries and name of the augmented graph
flank.l = 5000
cl.g = findOverlaps(svs, svs, maxgap=2 * flank.l + 1000) %>% as.data.frame %>%
  mutate(qid=names(svs)[queryHits], sid=names(svs)[subjectHits]) %>% 
  select(qid, sid) %>% as.matrix %>% graph_from_edgelist(directed=FALSE)
cl.c = components(cl.g)
cl.c.df = tibble(id=names(cl.c$membership), cmp=as.numeric(cl.c$membership))

svreg.df = tibble(id=names(svs), chr=as.character(seqnames(svs)),
                  start=start(svs), end=end(svs)) %>% merge(cl.c.df) %>%
  group_by(cmp, chr) %>% summarize(start=min(start)-flank.l, end=max(end)+flank.l)

## chr lengths
cyto.df = read.table('/public10/home/sci0011/data/ref/SL5.0/SL5.cyto.bed', as.is=TRUE, sep='\t')
cyto.df = cyto.df[, 1:3]
colnames(cyto.df) = c('chr', 'start', 'end')

aug.dir = 'chunks'
chunks.df = lapply(c(0:12), function(chrn){
  svreg.df = subset(svreg.df, chr==chrn)
  chunk.bks = sort(c(as.integer(svreg.df$start), as.integer(svreg.df$end),
                     0, max(subset(cyto.df, chr==chrn)$end)))
  tibble(chr=chrn, start=chunk.bks[-length(chunk.bks)], end=chunk.bks[-1]) %>%
    mutate(aug.vg=paste0(aug.dir, '_', chr, '/chunk_', chr, '_', start, '_', end, '.vg'))
}) %>% bind_rows()

write.table(chunks.df, file='graph_srdedup20_augment_chunks.bed', row.names=FALSE, col.names=FALSE, sep='\t', quote=FALSE)

## Write VCF to convert to fasta for augmentation
info(svs.for.aug)$GT = info(svs.for.aug)$END = info(svs.for.aug)$QUAL = info(svs.for.aug)$SVLEN = NULL
names(geno(svs.for.aug)) = character()
writeVcf(svs.for.aug, filename='graph_srdedup20_foraug.vcf')

## Write VCF for constructinfo(svs.for.cons)$GT = info(svs.for.cons)$END = info(svs.for.cons)$QUAL = info(svs.for.cons)$SVLEN = NULL
info(svs.for.cons)$GT = info(svs.for.cons)$END = info(svs.for.cons)$QUAL = info(svs.for.cons)$SVLEN = NULL
writeVcf(svs.for.cons, filename='graph_srdedup20_forcons.vcf')
