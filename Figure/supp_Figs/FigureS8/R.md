
library(ggplot2)
TEpal <- c("#5fc0a7","#abc457","#403dcc","#4771f5",
           "#4da3dc","#f5c456","#ef8b3b","#ec5a29","grey")

data <- read.table("/Users/zhiyangzhang/server/tmp/sum_repeat_rename2.xls",header=T)
data1<- read.table("/Users/zhiyangzhang/server/draw/46_genome_order.txt",header=T)
order <- rev(unique(data1$Name))
data$Name <- factor(data$Name,levels=order)
data$Repeat <- factor(data$Repeat,levels=c("DNA/TIR","DNA/Helitron","LTR/Copia","LTR/Gypsy","LTR/unknown","LINE","Low_complexity","Simple_repeat","Unknown"))

ggplot(data,aes(x=Name,y=Percentage))+
  geom_bar(stat = "identity",aes(fill=Repeat))+
  scale_fill_manual(values=TEpal)+
  theme_bw()+
  theme(
    axis.text=element_text(size=8,color="black")
  )+
  coord_flip()+
  labs(x="",y="Repeat Content (%)")+
  scale_y_continuous(expand=c(0,0),labels = scales::percent,limits=c(0,1))
ggsave("/Users/zhiyangzhang/server/tmp/Repeat-Content-Tomato.pdf",height = 7*0.618+2,width=7,units="in")

   data$Base <- data$Base/1e6

ggplot(data,aes(x=Name,y=Base))+
  geom_bar(stat = "identity",aes(fill=Repeat))+
  scale_fill_manual(values=TEpal)+
  theme_bw()+
  theme(
    axis.text=element_text(size=8,color="black")
  )+
  coord_flip()+
  labs(x="",y="Repeat Content (Mb)")+
  scale_y_continuous(expand=c(0,0),limits=c(0,800))
ggsave("/Users/zhiyangzhang/server/tmp/Repeat-Content-Bases-Tomato.pdf",height = 7*0.618+2,width=7,units="in")
