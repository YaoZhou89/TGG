library(ggplot2)
library(dplyr)

TEpal <- c("#5fc0a7","#abc457","#403dcc","#4771f5",
           "#4da3dc","#f5c456","#ef8b3b","#ec5a29","grey")
SVpal <- c("#5fc0a7","#4771f5","#f5c456","#ec5a29","grey80")
mypal <- c(TEpal[6],TEpal[5],TEpal[1])
data1<- read.table("/Users/zhiyangzhang/server/draw/46_genome_order.txt",header=T)
order <- (unique(data1$Name))

data <- read.table("/Users/zhiyangzhang/server/draw/sum_LAI_rename2.xls",header=T)
data$Name <- factor(data$Name,levels=order)
ggplot(data,aes(x=Name,y=LAI))+
  geom_boxplot(aes(fill=Group),outlier.size = 0.5)+
  scale_fill_manual(values=mypal)+
  theme_bw()+
  theme(
    axis.text.x=element_text(size=8,color="black",angle=90),
    axis.text.y = element_text(color="black",size=8)
  )+
  ylim(0,30)+
  labs(x="",y="LTR Assembly Index")
ggsave("/Users/zhiyangzhang/server/draw/LAI-Tomato.pdf",height = 7*0.618,width=7,units="in")

