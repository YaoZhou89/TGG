```{r}
library(ggplot2)
library(dplyr)

TEpal <- c("#5fc0a7","#abc457","#403dcc","#4771f5",
           "#4da3dc","#f5c456","#ef8b3b","#ec5a29","grey")
SVpal <- c("#5fc0a7","#4771f5","#f5c456","#ec5a29","grey80")
mypal <- c(TEpal[6],TEpal[5],TEpal[1])
data1<- read.table("/Users/zhiyangzhang/server/draw/46_genome_order.txt",header=T)
order <- (unique(data1$Name))
data <- read.table("/Users/zhiyangzhang/server/draw/46genome_RNA_stats.tsv",header = T)
data$Name <-factor(data$Name,levels=(order))
ggplot(data,aes(x=Name,y=Mapping_rate,fill=Group,order(data[,4])))+
  geom_boxplot()+
  theme_bw()+
  scale_fill_manual(values=mypal)+
  scale_y_continuous(labels = scales::percent_format(scale = 100))+
  labs(x="",y="83 SRA reads mapping rate")+
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    axis.text.y = element_text(size=6,color="black"),
    axis.text.x=element_text(size=6,color="black",angle=90)
    
  )

ggsave("/Users/zhiyangzhang/server/draw/46genome_RNA_stats2.pdf",height = 7*0.618+2,width=7,units="in")
```