```{r}
library(ggplot2)
library(dplyr)

TEpal <- c("#5fc0a7","#abc457","#403dcc","#4771f5",
           "#4da3dc","#f5c456","#ef8b3b","#ec5a29","grey")
SVpal <- c("#5fc0a7","#4771f5","#f5c456","#ec5a29","grey80")
mypal <- c(TEpal[6],TEpal[5],TEpal[1])
data1<- read.table("/Users/zhiyangzhang/server/draw/46_genome_order.txt",header=T)
my_order <- (unique(data1$Name))

library(tidyr)
library(dplyr)
data <- read.table("/Users/zhiyangzhang/server/draw/Tomato_46genomes.gene.stats2.tsv",header=F,sep="\t")
colnames(data) <- c("Name","Group","Stats","Number")
head(data)

newdata <- spread(data,Stats,Number)
write.table(newdata,file="/Users/zhiyangzhang/server/draw/Tomato_46genomes.gene.stats_long2.xls",sep="\t",quote=F)
plotdata <- filter(data,Stats=="Mean exon size" | Stats == "Mean gene locus size (first to last exon)" | Stats == "Mean number of distinct exons per gene" | Stats == "Mean transcript size (UTR CDS)" | Stats == "Number of genes")
plotdata$Stats <- gsub("Mean gene locus size (first to last exon)","Mean gene locus size",plotdata$Stats)
plotdata$Stats <- gsub("Mean transcript size (UTR CDS)","Mean transcript size",plotdata$Stats)
plotdata$Name <- factor(plotdata$Name,levels=my_order)
plotdata$Number <- as.numeric(plotdata$Number)

ggplot(plotdata)+
  geom_bar(stat="identity",aes(x=Name,y=Number,fill=Group))+
  facet_wrap(~Stats,ncol = 1,scales = "free_y")+
  theme_bw()+
  theme(
    axis.text.x  = element_text(size=6,color="black",angle=90),
    axis.text.y  = element_text(size=6,color="black")
  )+
  scale_fill_manual(values=mypal)+
  scale_y_continuous(breaks = scales::pretty_breaks(3), limits = c(0, NA))

ggsave("/Users/zhiyangzhang/server/draw/Tomato-46genomes-gene.stats.pdf",width=7,
       height=7*0.618,units="in")


```