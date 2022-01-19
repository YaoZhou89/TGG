### Figure1B_1C

1B- pangenome曲线

1C core / softcore统计

```{r}

library(dplyr)
data <- read.table("Finale_Figure/Figure1/1B_1C/Orthogroups_simulation_100.stats",header=T)
data$Type <- factor(data$Type,levels=c("Pan","Core"))


summ_data <- data %>%
   group_by(Number,Type) %>%
   summarise(qs = quantile(Gene, c(0.05, 0.95)),prob = c("p05", "p95"),mean=mean(Gene))

summ_data_wide <- spread(summ_data,key = prob,value=qs)


p1b <- ggplot()+
  geom_line(data=summ_data,aes(x=Number,y=mean,color=Type))+
  geom_ribbon(data=summ_data_wide,aes(x=Number,ymin=p05,ymax=p95,fill=Type),alpha=0.5)+
  theme_classic()+
  theme(
    axis.text = element_text(color="black",size=6),
    legend.title = element_blank(),
    legend.position = c(0.8,0.5),
    legend.background = element_blank(),
    legend.key.size = unit(1,"line"),
    legend.text = element_text(size=6),
    axis.line = element_line(size=0.3,color="black"),
    axis.title=element_text(color="black",size=6),
    axis.ticks = element_line(color="black",size=0.25)
  )+
  labs(x="# Sample",y="# Gene Family")+
  scale_y_continuous(labels = scales::comma,breaks = scales::pretty_breaks(n=6))+
  scale_fill_manual(values=c("#3288bd","#d53e4f"))+
  scale_color_manual(values=c("#3288bd","#d53e4f"))
#ggsave("./Finale_Figure/Figure1/1B_1C/Figure1B.pdf",width=8.9,height=4,units="cm")

library(patchwork)
data <- read.table("./Finale_Figure/Figure1/1B_1C/Orthogroups_freq.tsv",header=T)
head(data)
data$Type <- factor(data$Type,levels=c("Core","Softcore","Dispensable","Private"))

PANpal <- c(TEpal[7],TEpal[6],TEpal[5],TEpal[1])
p1 <- ggplot(data=data,aes(x=Freq,y=Number,fill=Type))+
  geom_bar(stat="identity")+
  geom_text(data=data,aes(label=Number), position=position_dodge(width =0.9),hjust=1,angle=90,size=1.5) +
  theme_classic()+
  theme(
    axis.text = element_text(color="black",size=6),
    legend.title = element_blank(),
    legend.position =  c(0.5,0.8),
    legend.key.size = unit(0.5,"line"),
    legend.text = element_text(size=6),
    legend.background = element_blank(),
    axis.line = element_line(size=0.3,color="black"),
    axis.title=element_text(color="black",size=6),
    axis.ticks = element_line(color="black",size=0.25)
  )+
  labs(x="Frequency",y="# Gene Family")+
  scale_x_continuous()+
  scale_y_continuous(expand=c(0,0),trans="log10",limits=c(1,50000),breaks=c(100,500,1000,5000,10000),labels=scales::comma)+
#breaks = scales::trans_breaks("log10", function(x) 10^x),labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  #annotation_logticks(sides="l",outside = TRUE,size=0.1)+
  scale_fill_manual(values=PANpal)+
  guides(fill = guide_legend(nrow = 1))

pie_data <- read.table("./Finale_Figure/Figure1/1B_1C/Orthogroups_stats.tsv",header=T)
gene_sum <- sum(pie_data$Number)
pie_data$Perc <- round(pie_data$Number/gene_sum*100,2)
pie_data$Type <- factor(pie_data$Type,levels=rev(c("Core","Softcore","Dispensable","Private")))

pie_data$ypos <- (cumsum(pie_data$Perc)-0.5*pie_data$Perc)

p2 <- ggplot(pie_data, aes(x="", y=Perc, fill=Type)) +
  geom_bar(stat="identity", width=1,position = "stack")+
  geom_text(aes(label=paste(Perc,"%",sep="")),position = position_stack(vjust = 0.5),size=2)+
  #coord_polar("y", start=0) +
  theme_void() + 
  coord_flip()+
  theme(legend.position="none")+
  #geom_text(aes(label = paste(Type," : ",Perc,"%",sep="")),position = position_stack(vjust = 0.5),color = "black", size=1.5) +
  scale_fill_manual(values=PANpal)

p1c <- p2 + p1 + plot_layout(height = c(1, 10))
p1c
p1bc <- p1b / p1c + plot_layout(height = c(1, 1))
p1abc <- p1a + p1bc
ggsave("./Finale_Figure/Figure1/Figure1abc-v2.pdf",p1abc,width=18.3,height=8.9,units="cm")
```