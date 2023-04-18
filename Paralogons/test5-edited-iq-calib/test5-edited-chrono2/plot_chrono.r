setwd("/Users/fmarletaz/Dropbox/Genomes/Hagfish/Gene_families/test5-edited-iq-calib/test5-edited-chrono2")
library(phytools)
library(tidyverse)
library(wesanderson)
library(RColorBrewer)


files <- list.files(path=".", pattern="*.chronogram",  recursive=FALSE)
pdf('clg_dating.pdf')
dev.off()
x=28

rotateNodes(clg,c(63,95,87))

clg<-read.tree("CLGB_ln_sample.chronogram")

#clg=rotateNodes(clg,c(63,95,87))

data.frame(clg$node.label) %>% separate(clg.node.label,c('lower','upper'),sep='_') %>%
  transform(upper=as.numeric(upper),lower=as.numeric(lower)) -> CI
plotTree.errorbars(clg,bar.width=7, CI,at=seq(0,600,by=50),bar.col="dodgerblue3",cex=0.6,fsize=0.8)
text(600,31,label=fsp[1],adj=0)


ages=data.frame()
for (f in files){
  print(f)  
  chrm=read.tree(f)
  #chrm=rotateNodes(chrm,c(53,54,47,39,37,33,34)) # causes problems with subsequent functions
  #chrm=rotateConstr(chrm,chrm$tip.label[28])
  #tree=rotateConstr(chrm,chrm$tip.label[x])
  #chrm$tip.label<-listTax
  data.frame(chrm$node.label) %>% separate(chrm.node.label,c('lower','upper'),sep='_') %>%
    transform(upper=as.numeric(upper),lower=as.numeric(lower)) -> CI
  plotTree.errorbars(chrm,bar.width=7, CI,at=seq(0,600,by=50),bar.col="dodgerblue3",fsize=0.7)
  H=round(unique(max(nodeHeights(chrm))-nodeHeights(chrm)[,1]),1)
  nodelabels(H,node=1:chrm$Nnode+Ntip(chrm),adj=c(1.1,-0.4),frame="none",cex=0.7,col="dodgerblue3")
  #nodelabels(1:chrm$Nnode,node=1:chrm$Nnode+Ntip(chrm),adj=c(1.1,-0.4),frame="none",cex=0.9,col="dodgerblue3")
  #sel_nodes=c(1,4,8,10,18,11,9,5)
  f0=gsub("_sample.chronogram", "",f)
  fsp=strsplit(f0,"_")[[1]]
  #f1=gsub("pteropods_May18_", "",f0) 
  text(600,0.25,label=fsp[1],adj=0)
  
  #node_age=data.frame(node=name_nodes,age=H[sel_nodes],lower=CI[sel_nodes,]$lower,
  #                    upper=CI[sel_nodes,]$upper,scheme=fsp[3],model=fsp[5])
  #print(node_age)
  #ages <- rbind(ages, node_age)
}
dev.off()
write.table(ages,"ages_clades_priors.txt",quote=F,sep='\t')
ages$node=factor(  ages$node,levels=name_nodes)
ggplot(ages, aes(node, age,colour=scheme,shape=model))+ geom_boxplot(color='grey',fill='ghostwhite') +
  geom_jitter(width=0.2,size=2.5)+theme_light()+#facet_grid(rows =vars(model))+
  scale_colour_brewer(palette = "Set2")+geom_hline(c()


ggplot(ages, aes(node, age,colour=scheme,shape=model))+ #geom_boxplot(color='grey',fill='ghostwhite') +
  theme_bw()+facet_grid(rows =vars(model))+
  scale_colour_brewer(palette = "Set2")+geom_pointrange(aes(ymin=lower, ymax=upper),position=position_jitter(width=0.2))


#scale_color_manual(values = c(wes_palette("Darjeeling2"),wes_palette("Zissou1")))

# REVISED 20/05/2019
packageVersion("phytools")
?plotTree
c='pteropods_May18_s7_cir_sb_bd_sample.chronogram'
l='pteropods_May18_s7_ln_sb_bd_sample.chronogram'
u='pteropods_May18_s7_ugam_sb_bd_sample.chronogram'

chrm=read.tree(u)
chrm=rotateNodes(chrm,c(53,54,47,39,37,33,34))
chrm=rotateConstr(chrm,chrm$tip.label[28])
chrm$tip.label<-listTax

data.frame(chrm$node.label) %>% separate(chrm.node.label,c('lower','upper'),sep='_') %>%
  transform(upper=as.numeric(upper),lower=as.numeric(lower)) -> CI

nodeheight(chrm,16)
chrm$edge 
nodeHeights(chrm)
chrm$ed
plotTree.errorbars(chrm,bar.width=7, CI,at=seq(0,250,by=50),bar.col="dodgerblue3",cex=0.6,fsize=0.8)

?plotTree.errorbars
H=round(unique(max(nodeHeights(chrm))-nodeHeights(chrm)[,1]),1)
H
max(nodeHeights(chrm))-nodeHeights(chrm)[,1]
nodeHeights(chrm)
chrm$tip.label
1:chrm$Nnode+Ntip(chrm)
1:chrm$Nnode+Ntip(chrm)
nodeHeights(chrm)
nodelabels(H,adj=c(1.1,-0.4),frame="none",cex=0.9,col="dodgerblue3")
nodelabels(1:chrm$Nnode,node=1:chrm$Nnode+Ntip(chrm),adj=c(1.1,-0.4),frame="none",cex=0.9,col="dodgerblue3")

chrm$node.label
sel_nodes=c(1,4,8,10,18,11,9,5)
H[sel_nodes]
CI[sel_nodes,]
node_age=data.frame(node=name_nodes,age=H[sel_nodes],lower=CI[sel_nodes,]$lower,upper=CI[sel_nodes,]$upper)
write.table(node_age,'ugam.txt',quote=F,sep='\t')
name_nodes=c('Root','Pteropoda','Thecosomata','Euthecosomata',
  'Cavolinoidea','Limancinoidea','Pseudothecosomata','Gymnosomata')
#f0=gsub("_sample.chronogram", "",f)
#f1=gsub("pteropods_May18_", "",f0)
text(100,28,label=f1,adj=0)
?text
### PREVIOUS STUFF ### 

add_column(int_node,!!(f1):=H[int_node$node-Ntip(tree)])-> int_node

chrm$tip.label

for (f in files){
   print(f)  
   chrm=read.tree(f)
   tree=rotateConstr(chrm,chrm$tip.label[x])
   data.frame(tree$node.label) %>% separate(tree.node.label,c('lower','upper'),sep='_') %>%
   transform(upper=as.numeric(upper),lower=as.numeric(lower)) -> CI
   plotTree.errorbars(tree,bar.width=7,cex=1, CI,at=seq(0,250,by=50),bar.col="dodgerblue3",cex=0.6)
   H=round(unique(max(nodeHeights(tree))-nodeHeights(tree)[,1]),1)
   nodelabels(H,node=1:tree$Nnode+Ntip(tree),adj=c(1.1,-0.4),frame="none",cex=0.8)
   f0=gsub("_sample.chronogram", "",f)
   f1=gsub("pteropods_May18_", "",f0)
   text(240,25,label=f1,adj=0)
   add_column(int_node,!!(f1):=H[int_node$node-Ntip(tree)])-> int_node
}
H[]
nodelabels(node=1:tree$Nnode+Ntip(tree),adj=c(1.1,-0.4),frame="none",cex=0.8)

H[int_node$node-Ntip(tree)]
int_node=data.frame(node=c(29,32,36,38,46,39,37,33),
                    name=c('Root','Pteropoda','Thecosomata','Euthecosomata',
                           'Cavolinoidea','Limancinoidea','Pseudothecosomata','Gymnosomata'))
int_node
write.table(int_node,file='ages_test.tsv',quote=F,sep='\t')
int_node %>% add_column(int_node,!!(f1):=H[int_node$node-Ntip(tree)])
dev.off()
H[,1]
nodelabels(round(unique(204.2615-nodeHeights(tree)[,1]),1),node=1:tree$Nnode+Ntip(tree),
          adj=c(1.1,-0.4),frame="none",cex=0.8)
H[29,1]
?