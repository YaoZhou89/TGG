```R
library(ggtree)
library(treeio)
library(tidytree)
mypal <- c("#f5c456","#4da3dc","#5fc0a7","#403dcc")
tree<-read.newick("C:/Users/Younger/Desktop/RAxML_bestTree.raxml.ml") 
ggtree(tree,branch.length = 'none',layout='circular') 
group_file <- read.table("C:/Users/Younger/Desktop/TOMATO/phylogeny/tree_group",header = T,row.names = 1)
groupInfo <- split(row.names(group_file), group_file$Group)
tree <- groupOTU(tree, groupInfo)
data2<-read.csv("C:/Users/Younger/Desktop/TOMATO/phylogeny/ID.txt",header = FALSE,sep = "\t")
c<-data2[data2$V2 =="select",]$V1
p2<-ggtree(tree,branch.length='none',layout='circular',aes(color=group)) + geom_point2(aes(subset=(label %in% c),fill=group),shape=21, size=1.5,hjust=1,legend.position = "none") +geom_tiplab(hjust = -0.6,size=0.5)+scale_fill_manual(values =c("#f5c456","#4da3dc","#5fc0a7","#403dcc"))
p2
ggsave("C:/Users/Younger/Desktop/main_figure/supp/phylogeny.pdf",p2,width=28,height=40,units="cm")
```

