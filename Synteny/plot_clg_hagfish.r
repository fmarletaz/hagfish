require(RColorBrewer)
require(pheatmap)
require(tidyverse)

plotSynt <- function(synt,labx,laby){
  synt  %>% select(s2gi,s2chr) %>%  group_by(s2chr) %>%  summarise(max=max(s2gi))%>%
    arrange(max) %>% mutate(mids=lag(max,default=0)+0.5*(max-lag(max,default=0)))-> s2chLim
  synt  %>% select(s1gi,s1chr) %>%  group_by(s1chr) %>%  summarise(max=max(s1gi))%>%
    arrange(max) %>% mutate(mids=lag(max,default=0)+0.5*(max-lag(max,default=0)))-> s1chLim
  #palette(colorRampPalette(brewer.pal(6,'BuPu'))(nrow(s2chLim)))
  #palette(colorRampPalette(brewer.pal(12,'Paired'))(nrow(s1chLim)))
  plab=paste("(",nrow(synt)," orthologues)",sep="")
  labxp=paste(labx,plab,sep=" ")
  #plot(synt$s1gi,synt$s2gi,pch=20,cex=0.4,col=synt$s2chr,xlab=labxp,ylab=laby)
  plot(synt$s1gi,synt$s2gi,pch=20,cex=0.4,col=synt$s1chr,xlab=labx,ylab=laby,xlim=c(0,max(s1chLim$max)))
  plot(synt$s1gi,synt$s2gi,pch=20,cex=0.4,col=ifelse(synt$scol=='S',as.factor(synt$s1chr),'lightgrey'),xlab=labx,ylab=laby,xlim=c(0,max(s1chLim$max)))
  #plot(synt$s1gi,synt$s2gi,pch=20,cex=0.4,col=ifelse(synt$scol=='S','#810f7c','lightgrey'),xlab=labx,ylab=laby,xlim=c(0,max(s1chLim$max)))
  
  abline(v=0,lty=3,lwd=0.5)
  abline(h=0,lty=3,lwd=0.5)
  
  abline(v=s1chLim$max,lty=3,lwd=0.5)
  abline(h=s2chLim$max,lty=3,lwd=0.5)
  
  axis(3,at=s1chLim$mids,labels=s1chLim$s1chr,las=2,cex.axis=0.6)
  axis(4,at=s2chLim$mids,labels=s2chLim$s2chr,las=1,cex.axis=0.6)
  
}
plotSyntCLGr <- function(synt,labx,laby){
  synt  %>% select(s2gi,s2chr) %>%  group_by(s2chr) %>%  summarise(max=max(s2gi))%>%
    arrange(max) %>% mutate(mids=lag(max,default=0)+0.5*(max-lag(max,default=0)))-> s2chLim
  synt  %>% select(s1gi,s1chr) %>%  group_by(s1chr) %>%  summarise(max=max(s1gi))%>%
    arrange(max) %>% mutate(mids=lag(max,default=0)+0.5*(max-lag(max,default=0)))-> s1chLim
  palette(c(brewer.pal(10,'Paired'),brewer.pal(7,'Dark2')))
  #palette(colorRampPalette(brewer.pal(12,'Paired'))(length(levels(synt$clg))))
  #palette(colorRampPalette(brewer.pal(12,'Paired'))(nrow(s1chLim)))
  
  #plot(synt$s1gi,synt$s2gi,pch=20,cex=0.4,col=synt$clg,xlab=labx,ylab=laby)
  #plot(synt$s2gi,synt$s1gi,pch=20,cex=0.4,col=synt$clg,xlab=labx,ylab=laby)
  plot(synt$s2gi,synt$s1gi,pch=20,cex=0.4,col=ifelse(synt$scol=='S',synt$clg,'lightgrey'),xlab=labx,ylab=laby)
  
  abline(v=0,lty=3,lwd=0.5)
  abline(h=0,lty=3,lwd=0.5)
  
  abline(h=s1chLim$max,lty=3,lwd=0.5)
  abline(v=s2chLim$max,lty=3,lwd=0.5)
  
  axis(4,at=s1chLim$mids,labels=s1chLim$s1chr,las=1,cex.axis=0.6)
  axis(3,at=s2chLim$mids,labels=s2chLim$s2chr,las=2,cex.axis=0.6)
  
}
plotSyntCLGoc <- function(synt,labx,laby,c){
  synt  %>% select(s2gi,s2chr) %>%  group_by(s2chr) %>%  summarise(max=max(s2gi))%>%
    arrange(max) %>% mutate(mids=lag(max,default=0)+0.5*(max-lag(max,default=0)))-> s2chLim
  synt  %>% select(s1gi,s1chr) %>%  group_by(s1chr) %>%  summarise(max=max(s1gi))%>%
    arrange(max) %>% mutate(mids=lag(max,default=0)+0.5*(max-lag(max,default=0)))-> s1chLim
  #palette(c(brewer.pal(10,'Paired'),brewer.pal(7,'Dark2')))
  #palette(colorRampPalette(brewer.pal(12,'Paired'))(length(levels(synt$clg))))
  #palette(colorRampPalette(brewer.pal(12,'Paired'))(nrow(s1chLim)))
  
  plot(synt$s1gi,synt$s2gi,pch=20,cex=0.4,col=ifelse(synt$scol=='S',synt$clg,'lightgrey'),xlab=labx,ylab=laby)
  #plot(synt$s1gi,synt$s2gi,pch=20,cex=0.4,col=synt$s1chr,xlab=labx,ylab=laby)
  
  abline(v=0,lty=3,lwd=0.5)
  abline(h=0,lty=3,lwd=0.5)
  
  abline(v=s1chLim$max,lty=3,lwd=0.5)
  abline(h=s2chLim$max,lty=3,lwd=0.5)
  
  axis(3,at=s1chLim$mids,labels=s1chLim$s1chr,las=2,cex.axis=0.6)
  axis(4,at=s2chLim$mids,labels=s2chLim$s2chr,las=1,cex.axis=0.6)
  
}


ordSynt <- function(synt){
  synt %>% dplyr::count(s1chr,s2chr) %>% spread(s2chr,n,fill=0) %>% column_to_rownames(var="s1chr")-> chr_mat
  hmp<-pheatmap(as.matrix(chr_mat))
  y_chrom<-rownames(as.matrix(chr_mat))[hmp$tree_row$order]
  x_chrom<-colnames(as.matrix(chr_mat))[hmp$tree_col$order]
  #print(x_chrom)
  synt %>% mutate(s1chrO=match(s1chr,y_chrom)) %>% 
    arrange(s1chrO,s1gp) %>% mutate(s1gi=row_number(s1chrO)) %>%
    mutate(s2chrO=match(s2chr,x_chrom)) %>% 
    arrange(s2chrO,s2gp) %>%  mutate(s2gi=row_number(s2chrO)) -> synt.ord
  return(synt.ord)
  #return(c(x_chrom,y_chrom))
}

ordSyntRw <- function(synt){
  synt %>% count(s1chr,s2chr) %>% spread(s2chr,n,fill=0) %>% column_to_rownames(var="s1chr")-> chr_mat
  hmp<-pheatmap(as.matrix(chr_mat))
  #y_chrom<-rownames(as.matrix(chr_mat))[hmp$tree_row$order]
  x_chrom<-colnames(as.matrix(chr_mat))[hmp$tree_col$order]
  #print(x_chrom)
  synt %>% 
    #mutate(s1chrO=match(s1chr,y_chrom)) %>% 
    #arrange(s1chrO,s1gp) %>% mutate(s1gi=row_number(s1chrO)) %>%
    mutate(s2chrO=match(s2chr,x_chrom)) %>% 
    arrange(s2chrO,s2gp) %>%  mutate(s2gi=row_number(s2chrO)) -> synt.ord
  return(synt.ord)
  #return(c(x_chrom,y_chrom))
}



ordSyntJO <- function(synt){
  synt %>% count(s1chr,s2chr) %>% spread(s2chr,n,fill=0) %>% column_to_rownames(var="s1chr")-> chr_mat
  hmp<-pheatmap(as.matrix(chr_mat))
  y_chrom<-rownames(as.matrix(chr_mat))[hmp$tree_row$order]
  x_chrom<-colnames(as.matrix(chr_mat))[hmp$tree_col$order]
  #print(x_chrom)
  synt %>% mutate(s1chrO=match(s1chr,y_chrom)) %>% 
    arrange(s1chrO,s1gp) %>% mutate(s1gi=row_number(s1chrO)) %>%
    mutate(s2chrO=match(s2chr,x_chrom)) %>% 
    arrange(s2chrO,s2gp) %>%  mutate(s2gi=row_number(s2chrO)) -> synt.ord
  #return(synt.ord)
  return(list(x_chrom,y_chrom))
}

testEnrichCLG <- function(synt){
  synt %>% 
    group_by(clg) %>% mutate(clgtot=n()) %>%
    group_by(s2chr) %>% mutate(chrtot=n())%>%
    group_by(s2chr,clgtot,chrtot,clg) %>% tally() %>% arrange(clg) -> synt.chrn
  synt.no=sum(synt.chrn$n)
  synt.chrn %>% rowwise() %>% 
    mutate(fishp=fisher.test(matrix(c(n,clgtot-n,chrtot-n,synt.no),ncol=2),alternative='greater')$p.value) %>%
    ungroup() %>% mutate(padj=p.adjust(fishp,method="bonferroni")) %>% mutate(scol=ifelse(padj<.05,'S','NS')) -> synt.test
  return(synt.test)
}


testEnrich <- function(synt){
  synt %>% 
    group_by(s1chr) %>% mutate(s1tot=n()) %>%
    group_by(s2chr) %>% mutate(s2tot=n())%>%
    group_by(s1chr,s2chr,s1tot,s2tot) %>% tally() -> synt.chrn
  synt.no=sum(synt.chrn$n)
  synt.chrn %>% rowwise() %>% 
    mutate(fishp=fisher.test(matrix(c(n,s1tot-n,s2tot-n,synt.no),ncol=2),alternative='greater')$p.value) %>%
    ungroup() %>% mutate(padj=p.adjust(fishp,method="bonferroni")) %>% mutate(scol=ifelse(padj<.05,'S','NS')) -> synt.test
  return(synt.test)
} 
gchkExp <- function(synt.exp){
  synt.exp %>% arrange(s1chr,s1gi) %>% #filter(s1chr=="Scaffold_2100") %>% 
    group_by(s1chr) %>% mutate(chunk=cut_width(s1gi,width=20)) %>% ungroup() %>% filter(scol=="S") %>% 
    select(s2chr,chunk,s1chr,s1gp) %>% 
    group_by(s1chr,chunk) %>% mutate(xmin=min(s1gp),xmax=max(s1gp)) %>%
    group_by(s2chr,s1chr,chunk,xmin,xmax) %>% tally() %>% 
    group_by(chunk) %>% mutate(ymin=(cumsum(n)-n)/sum(n),ymax=cumsum(n)/(sum(n))) -> synt.gchk
  return(synt.gchk)
}

gchkExpCLG <- function(synt.exp){
  synt.exp %>% arrange(s2chr,s2gi) %>% #filter(s1chr=="Scaffold_2100") %>% 
    group_by(s2chr) %>% mutate(chunk=cut_width(s2gi,width=20)) %>% ungroup() %>% filter(scol=="S") %>% 
    select(clg,chunk,s2chr,s2gp) %>% 
    group_by(s2chr,chunk) %>% mutate(xmin=min(s2gp),xmax=max(s2gp)) %>%
    group_by(clg,s2chr,chunk,xmin,xmax) %>% tally() %>% 
    group_by(chunk) %>% mutate(ymin=(cumsum(n)-n)/sum(n),ymax=cumsum(n)/(sum(n))) -> synt.gchk
    synt.gchk$s2chr=fct_relevel(synt.gchk$s2chr,str_sort(levels(synt.gchk$s2chr), numeric = TRUE))
  return(synt.gchk)
}


par(mar=c(5,5,5,5))


setwd("~/Dropbox/Genomes/Synteny/synt-forge/Hagfish-redux/")


bflo_col=c("#1F78B4","#1B9E77","#33A02C","#B2DF8A","#FB9A99","#E31A1C","#FFFF99","#B15928","#D95F02","#7570B3","#E7298A","#66A61E","#CAB2D6","#6A3D9A","#E6AB02","#A6761D","#8DD3C7","#666666","#FF7F00","#FDBF6F","#BEBADA","#ffd92f")
palette(bflo_col)
bflo_col_red=c("#1F78B4","#33A02C","#E31A1C","#FFFF99","#B15928","#D95F02","#7570B3","#E7298A","#66A61E","#6A3D9A","#E6AB02","#A6761D","#8DD3C7","#666666","#FF7F00","#BEBADA","#ffd92f")

bflo_col_red=c("#1F78B4","#33A02C","#E31A1C","#FFFF99","#B15928","#D95F02","#7570B3","#E7298A","#66A61E","#6A3D9A","#E6AB02","#A6761D","#8DD3C7","#666666","#FF7F00","#BEBADA","#ffd92f")
palette(bflo_col_red)

bflo_pata_synt=read.table("Bflo_Pata_synt3.txt",h=T,stringsAsFactors = T)
bflo_pata_synt %>% mutate(clg=gsub('[0-9]+', '',clg)) -> bflo_pata_synt
#bflo_pata_synt %>% row_update(s2chr %in% c('chr5','chr16') & clg=='A')
bflo_pata_synt.test=testEnrichCLG(bflo_pata_synt)
bflo_pata_synt %>% inner_join(bflo_pata_synt.test) -> bflo_pata_synt.exp
plotSyntCLGoc(ordSynt(bflo_pata_synt.exp),'Branchiostoma','Paramyxine')
gchkExpCLG(bflo_pata_synt.exp) -> bflo_pata_synt.gchk
bflo_pata_synt.gchk$s2chr=fct_relevel(bflo_pata_synt.gchk$s2chr,str_sort(levels(bflo_pata_synt.gchk$s2chr), numeric = TRUE))#-> bflo_pata_synt.gchk
bflo_pata_synt.gchk %>% group_by(s2chr) %>% summarise(max=max(xmax)) -> bflo_pata_synt.sizes
ggplot(bflo_pata_synt.gchk, aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax, fill = clg))+
  geom_rect(data = bflo_pata_synt.sizes, aes(xmax = max),xmin=0,ymin = 0, ymax = 1,fill='grey',alpha=0.5) +
  geom_rect()+scale_fill_manual(values = bflo_col)+ facet_grid(s2chr ~ .) + theme_minimal() +
  theme(strip.text.y = element_text(angle = 360)) 

bflo_pmar_synt=read.table("Bflo_Pmar_synt2.txt",h=T,stringsAsFactors = T)
bflo_pmar_synt %>% mutate(clg=gsub('[0-9]+', '',clg)) -> bflo_pmar_synt
bflo_pmar_synt.test=testEnrichCLG(bflo_pmar_synt)
bflo_pmar_synt %>% inner_join(bflo_pmar_synt.test) -> bflo_pmar_synt.exp
plotSyntCLGoc(ordSynt(bflo_pmar_synt.exp),'Branchiostoma','Petromyzon')

bflo_pmar_synt %>% filter(clg=='A')

gchkExpCLG(bflo_pmar_synt.exp) -> bflo_pmar_synt.gchk
bflo_pmar_synt.gchk %>% group_by(s2chr) %>% summarise(max=max(xmax)) -> bflo_pmar_synt.sizes
bflo_pmar_synt.gchk$s2chr=fct_relevel(bflo_pmar_synt.gchk$s2chr,bflo_pmar_synt.chord[[1]])#-> bflo_pata_synt.gchk
ggplot(bflo_pmar_synt.gchk, aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax, fill = clg))+
  geom_rect(data = bflo_pmar_synt.sizes, aes(xmax = max),xmin=0,ymin = 0, ymax = 1,fill='grey',alpha=0.5) +
  geom_rect()+scale_fill_manual(values = bflo_col2)+ facet_grid(s2chr ~ .) + theme_minimal() +
  theme(strip.text.y = element_text(angle = 360)) 
bflo_pmar_synt.chord=ordSyntJO(bflo_pmar_synt.exp)

##### RETENTION RATES #####
#skChrAsgn=read.table("skate_clg_assign.txt",stringsAsFactors = T)
#names(skChrAsgn)<-c('Chr','Cat','Lab')

head(bflo_pmar_synt)
bflo_pmar_synt %>% mutate(clg=gsub(pattern = "[1-3]", replacement = "", clg)) %>% 
  group_by(clg) %>% tally %>%  column_to_rownames(var='clg')-> clg_ct


bflo_pmar_synt.exp %>% mutate(clg=gsub(pattern = "[1-3]", replacement = "", clg)) %>% 
  filter(scol=='S') %>% group_by(clg,s2chr) %>% tally -> sig_clg_ct
sig_clg_ct %>% mutate(ret=n/clg_ct[as.vector(clg),])-> sig_clg_ret

ggplot(sig_clg_ret,aes(x=ret,fill=clg))+geom_histogram(alpha = 0.8)+
  geom_density(alpha = 0.1)+scale_fill_manual(values=bflo_col_red)+
  theme_light()+facet_wrap(~clg,scales='free_y')


#skChrAsgn %>% column_to_rownames(var='Chr') -> skChrAsgn
sk_clg_ct %>% mutate(cat=skChrAsgn[as.vector(s2chr),]$Cat) %>%
  separate(cat, c("l","clgf","Evt"), sep = "_") %>% select(-l) %>%
  mutate(prz=grepl(clg,clgf)) %>% filter(prz==T) -> sk_clg_ct


sk_clg_ct %>% mutate(dup=str_sub(Evt,-1,-1)) -> sk_clg_ct

sk_clg_ct %>% group_by(s2chr,clgf,Evt) %>% filter(ret==max(ret)) -> sk_clg_flt
ggplot(sk_clg_ret,aes(x=ret,fill=dup,color=dup))+geom_density(alpha = 0.5)+theme_light()
ggplot(sk_clg_ret,aes(x=ret,fill=dup,color=dup))+geom_histogram(alpha = 0.5)+facet_grid(rows=vars(dup))

ggplot(data = sk_clg_ct,aes(axis1 = clg, axis2 = s2chr,y = n)) +
  scale_x_discrete(limits = c("CLG", "Chr"), expand = c(.2, .05)) +
  geom_alluvium(aes(fill = clg)) +
  geom_stratum(width = 1/4, fill = "ivory2", color = "black",decreasing = TRUE)+
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  theme_minimal()



ggplot(pata_pmar_clg.sel,aes(y=n,axis1=clg2,axis2=s2chr,axis3=s1chr))+
  geom_alluvium(aes(fill = clg),width = 1/4,decreasing = TRUE)+ 
  geom_stratum(width = 1/4, fill = "ivory2", color = "black",decreasing = TRUE)+
  scale_x_discrete(limits = c("CLG", "Lamprey", "Hagfish"), expand = c(.2, .05)) +
  stat_stratum(geom = "text", size=3,aes(label = after_stat(stratum)), decreasing = TRUE)+
  #geom_text(stat = "stratum", size=3,aes(label = after_stat(stratum))) +
  scale_fill_manual(values =cols)+theme_bw()



bflo_locu_synt=read.table("Bflo_Locu_synt2.txt",h=T,stringsAsFactors = T)
bflo_locu_synt %>% mutate(clg=gsub('[0-9]+', '',clg)) -> bflo_locu_synt

bflo_locu_synt.test=testEnrichCLG(bflo_locu_synt)
bflo_locu_synt %>% inner_join(bflo_locu_synt.test) ->bflo_locu_synt.exp
plotSyntCLGoc(ordSynt(bflo_locu_synt.exp),'Branchiostoma','Lepidosus')
bflo_locu_synt.exp[(bflo_locu_synt.exp$s2chr=='LG4' & bflo_locu_synt.exp$s1chr %in% c('BFL_10','BFL_16','BFL_18')),]$scol <- 'S'
bflo_locu_synt.exp[(bflo_locu_synt.exp$s2chr=='LG1' & bflo_locu_synt.exp$s1chr %in% c('BFL_5')),]$scol <- 'S'
bflo_locu_synt.exp[(bflo_locu_synt.exp$s2chr=='LG2' & bflo_locu_synt.exp$s1chr %in% c('BFL_1','BFL_12','BFL_14')),]$scol <- 'S'
bflo_locu_synt.exp[(bflo_locu_synt.exp$s2chr=='LG3' & bflo_locu_synt.exp$s1chr %in% c('BFL_5')),]$scol <- 'S'
bflo_locu_synt.exp[(bflo_locu_synt.exp$s2chr=='LG14' & bflo_locu_synt.exp$s1chr %in% c('BFL_7')),]$scol <- 'S'
gchkExpCLG(bflo_locu_synt.exp) -> bflo_locu_synt.gchk
bflo_locu_synt.gchk %>% group_by(s2chr) %>% summarise(max=max(xmax)) -> bflo_locu_synt.sizes

arrange(bflo_locu_synt.sizes,-max)%>% select(s2chr) -> Loc_oS
bflo_pmar_synt.chord[[1]]
bflo_locu_synt.gchk$s2chr=fct_relevel(bflo_locu_synt.gchk$s2chr,as.vector(Loc_oS[[1]]))-> bflo_pata_synt.gchk

ggplot(bflo_locu_synt.gchk, aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax, fill = clg))+
  geom_rect(data = bflo_locu_synt.sizes, aes(xmax = max),xmin=0,ymin = 0, ymax = 1,fill='grey',alpha=0.5) +
  geom_rect()+scale_fill_manual(values = bflo_col_red)+ facet_grid(s2chr ~ .) + theme_minimal() +
  theme(strip.text.y = element_text(angle = 360)) 

pata_pmax_synt=read.table("Pata_Pmax_synt.txt",h=T,stringsAsFactors = T)
pata_pmax_synt %>% rename(s1g=s2g,s1gp=s2gp,s1chr=s2chr,s1gi=s2gi,s2g=s1g,s2gp=s1gp,s2chr=s1chr,s2gi=s1gi) -> pata_pmax_synt
pata_pmax_synt.test=testEnrichCLG(pata_pmax_synt)
pata_pmax_synt %>% inner_join(pata_pmax_synt.test) -> pata_pmax_synt.exp
plotSyntCLGoc(ordSynt(pata_pmax_synt.exp),'Branchiostoma','Lepidosus')

#### Hagfish vs Lamprey ####
pata_pmar_synt=read.table("Pata_Pmar2_synt.txt",h=T,stringsAsFactors = T)
pata_pmar_synt %>% rename(s1g=s2g,s1gp=s2gp,s1chr=s2chr,s1gi=s2gi,s2g=s1g,s2gp=s1gp,s2chr=s1chr,s2gi=s1gi) -> pata_pmar_synt
pata_pmar_synt.test=testEnrich(pata_pmar_synt)
pata_pmar_synt %>% inner_join(pata_pmar_synt.test) -> pata_pmar_synt.exp
write.table(pata_pmar_synt.exp,file="Pata_Pmar_syntExp.tsv",quote=F,sep='\t')
head(pata_pmar_synt.exp)

pata_pmar_chOrd=ordSyntJO(pata_pmar_synt.exp)
bflo_pmar_synt.exp %>% filter(scol=='S') %>% select(s2g,s2chr,clg) %>% mutate(s1g=s2g,sxchr=s2chr) %>% select(-s2g,-s2chr) -> bflo_pmar_synt.expf
bflo_pmar_synt.exp %>% select(s2g,s2chr,clg) %>% mutate(s1g=s2g,sxchr=s2chr) %>% select(-s2g,-s2chr) -> bflo_pmar_synt.expf

left_join(pata_pmar_synt.exp,bflo_pmar_synt.expf) %>% relocate(clg,.after=sxchr) -> pata_pmar_synt.expr
?left_join
write.table(pata_pmar_synt.expr,file='Pata_Pmar_syntExp.txt',sep='\t',quote=F)


