```{r}
library(ggplot2)
library(tidytree)
library(treeio)
library(patchwork)
pal<-c("#C2749E","#57AADB")
data<-read.csv("C:/Users/Younger/Desktop/TOMATO/SNP_evaluation/indel_distance",header = FALSE,sep = "\t")
p1<-ggplot(data)+geom_histogram(aes(x=V1,fill=V2),binwidth = 10,alpha=0.8)+theme_classic()+scale_x_continuous(expand=c(0,0),limits = c(-1000,1000))+scale_y_continuous(expand=c(0,0))+labs(x="Distance",y="Indel counts")+theme_classic()+theme(legend.position = "none",
   axis.title.y = element_text(size=6,color="black"),                                                                      axis.title.x = element_text(size=6,color="black"),
    axis.text= element_text(size=6,color="black"),
    axis.line = element_line(size=0.3),
    axis.ticks = element_line(size=0.3)
  )+scale_fill_manual(values=pal)

data<-read.csv("C:/Users/Younger/Desktop/TOMATO/SNP_evaluation/snp_distance",header = FALSE,sep = "\t")
p2<-ggplot(data)+geom_histogram(aes(x=V1,fill=V2),binwidth =10,position="identity",alpha=0.8)+theme_classic()+scale_x_continuous(expand=c(0,0),limits = c(-1000,1000))+scale_y_continuous(expand=c(0,0))+labs(x="Distance",y="SNP counts")+theme(axis.title.x = element_text(size=6,color="black"),
axis.title.y = element_text(size=6,color="black"),
axis.text= element_text(size=6,color="black"),
axis.line = element_line(size=0.3),
axis.ticks = element_line(size=0.3)
)+scale_fill_manual(values=pal)
p1|p2
p3<-p1|p2
ggsave("C:/Users/Younger/Desktop/indel_snp_evluation2.pdf",p3,width=16,height=6,units="cm")
```