library(ggtree)
library(treeio)
library(phytools)
require(tidyverse)
setwd("/Users/fmarletaz/Dropbox (UCL)/Genomes/Hagfish/Gene_families/test5-edited/test5-edited-iq/")

files <- list.files(path=".", pattern="*.contree",  recursive=FALSE)
print(files)

print(files)
MRCA(clgc, 'A', 'E')
ggtree(clgc) + geom_text2(aes(subset=!isTip, label=node), hjust=-.3) + geom_tiplab()
og<-MRCA(clgc, .node1='Strpur_NA',.node2='Braflo_NA')
og
root(clgc, node=og, edgelabel = TRUE)
clg<-read.tree("CLGA_C20.contree")
og<-MRCA(clg, .node1='Strpur_NA',.node2='Braflo_NA')
clgr<-root(clg, node=og, edgelabel = TRUE)
ggplot(clgr) + geom_tiplab(size=3) + theme(legend.position="right") +
  geom_nodelab(label=clgr$node.label,hjust=0, size=2.5)+  hexpand(.1, direction = 1)+#xlim(NA,clgr$edge.length)/5)+
  geom_treescale()+
  geom_tree() + theme_tree()+ggtitle(f)

#clgc%>% as.treedata %>% as_tibble
pdf(file="paralogons_t5e_c20.pdf")
for (f in files){
  t<-read.tree(f)
  #t$tip.label=='Braflo_NA'
  #rn<-which(t$tip.label=="Braflo_NA")
  #len <- t$edge.length[t$edge[, 2] == rn]
  #rt<-root(t,rn,len/10,edgelabel = TRUE)
  og<-MRCA(t, .node1='Strpur_NA',.node2='Braflo_NA')
  rt<-root(t, node=og, edgelabel = TRUE)
  prt<-ggplot(rt) + geom_tiplab(size=3) + theme(legend.position="right") +
    geom_nodelab(label=rt$node.label,hjust=0, size=2.5)+ hexpand(.1, direction = 1)+#xlim(0,sum(rt$edge.length)/5)+
    geom_tree() + theme_tree()+ggtitle(f)#+theme(plot.margin = unit(c(25,25,25,25), "mm"))
  print(prt)
}
dev.off()
rt$node.label
?pdf
sum(rt$edge.length)

clg<-read.tree("CLGB_ln_sample.chronogram")
ggplot(clg)
spln=read.tree("reCycl01_cir_0_sample.chronogram")
spln=read.tree("reCycl01_0_sample.chronogram")


#reCycl01_cir_0_sample.chronogram
clg=spln
chrm=read.tree(u)
chrm=rotateNodes(chrm,c(53,54,47,39,37,33,34))
chrm=rotateConstr(chrm,chrm$tip.label[28])
chrm$tip.label<-listTax

clg=rotateNodes(clg,c(63,95,87))

data.frame(clg$node.label) %>% separate(clg.node.label,c('lower','upper'),sep='_') %>%
  transform(upper=as.numeric(upper),lower=as.numeric(lower)) -> CI

nodeheight(chrm,16)
chrm$edge 
nodeHeights(chrm)
chrm$ed
plotTree.errorbars(clg,bar.width=7, CI,at=seq(0,600,by=50),bar.col="dodgerblue3",cex=0.6,fsize=0.8)
H=round(unique(max(nodeHeights(clg))-nodeHeights(clg)[,1]),1)
nodelabels(H,node=1:clg$Nnode+Ntip(clg),adj=c(1.1,-0.4),frame="none",cex=0.7,col="dodgerblue3")
?nodelabels
