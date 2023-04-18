require(topGO)
require(tidyverse)
require(wesanderson)
require(RColorBrewer)
require(pheatmap)
library(viridis)

go_enrich <- function(selGenes,gene2GO) {
  geneList <- factor(as.integer(names(gene2GO) %in% selGenes))
  names(geneList) <- names(gene2GO)
  print(table(geneList))
  ontoList=list()
  ontologies=c("BP","CC","MF")
  for (i in 1:1) {
    print(ontologies[i])
    GOdata <- new("topGOdata", ontology = ontologies[i], allGenes = geneList,
                  annot = annFUN.gene2GO, gene2GO = gene2GO)
    resFischer <- runTest(GOdata, algorithm = "classic", statistic = "Fisher" )
    resGenes <- GenTable(GOdata,Fisher.classic=resFischer,orderBy ="Fisher.classic",topNodes=25,numChar=1000)
    resGenes %>% mutate(Fisher.classic = gsub("<","",Fisher.classic)) -> resGenes
    #resGenes %>% mutate(Fisher.classic = replace(Fisher.classic, Fisher.classic == "<1e-30", "1e-30")) -> resGenes
    resGenes$FDR=p.adjust(resGenes$Fisher.classic,method='fdr')
    resGenes$BH=p.adjust(resGenes$Fisher.classic,method='hochberg')
    resGenes$Ontology=ontologies[i]
    ontoList[[i]]<-resGenes
  }
  allRes=bind_rows(ontoList)
  return(allRes)
}


go_enrich_bp <- function(selGenes,gene2GO) {
  geneList <- factor(as.integer(names(gene2GO) %in% selGenes))
  names(geneList) <- names(gene2GO)
  print(table(geneList))
  #ontoList=list()
  #ontologies=c("BP","CC","MF")

  GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList,
                annot = annFUN.gene2GO, gene2GO = gene2GO)
  res <- runTest(GOdata, algorithm = "elim", statistic = "Fisher" )
  resGenes <- GenTable(GOdata,Fisher = res, orderBy = "Fisher",topNodes=1000,numChar=1000)
  resGenes %>% mutate(Fisher = gsub("< ","",Fisher)) -> resGenes
  genes <- sapply(resGenes$GO.ID, function(x)
  {
    genes<-genesInTerm(GOdata, x)
    genes[[1]][genes[[1]] %in% as.vector(selGenes)] # myGenes is the queried gene list
  })
  resGenes$genes<-unlist(lapply(genes,paste, collapse=','))
  #allGO = genesInTerm(GOdata)
  #selGO = lapply(allGO,function(x) x[x %in% selGenes] )  
  #resGenes %>% mutate(genes=paste(selGO[sep=','))
  #resGenes %>% mutate(Fisher.classic = replace(Fisher.classic, Fisher.classic == "<1e-30", "1e-30")) -> resGenes
  #resGenes$FDR=p.adjust(resGenes$Fisher.classic,method='fdr')
  #resGenes$BH=p.adjust(resGenes$Fisher.classic,method='hochberg')
  #resGenes$Ontology=ontologies[i]
  #ontoList[[i]]<-resGenes
  #allRes=bind_rows(ontoList)
  return(resGenes)
}

go_genes <- function(selGenes,gene2GO,term) {
  geneList <- factor(as.integer(names(gene2GO) %in% selGenes))
  names(geneList) <- names(gene2GO)
  print(table(geneList))
  #ontoList=list()
  #ontologies=c("BP","CC","MF")
  
  GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList,
                annot = annFUN.gene2GO, gene2GO = gene2GO)
  ggt<-as.character(unlist(genesInTerm(object = GOdata, whichGO = term)))
  cgt<-intersect(ggt,selGenes)
  #allGO = genesInTerm(GOdata)
  #selGO = lapply(allGO,function(x) x[x %in% selGenes] )
  
  return(cgt)
}


plot_go <- function(go_res,n_slice=20) {
  go_res %>% filter(BH<0.01 & Ontology=="BP") %>%
    mutate(MLP=-log10(BH)) %>%
    group_by(Ontology) %>% slice_head(n=n_slice) -> go.filt
  go.filt$Ontology=factor(go.filt$Ontology,levels=c('BP','MF'))
  go.filt$Term=factor(go.filt$Term,levels=go.filt[order(go.filt$Ontology),]$Term)
  ggplot(go.filt,aes(x=sort(MLP),y=Term)) + geom_bar(stat='identity')+
    scale_fill_manual(values = wes_palette("Royal1"))+theme_light()+ggtitle('NewDup')
  return(go.filt)
}

setwd("~/Dropbox/Genomes/Hagfish/Gene_families/PFAM_GO/")

locuGO <- readMappings(file = "Lepocu.emap.go.slim.txt")
pmarGO <- readMappings(file = "Petmar.emap.go.slim.txt")
pataGO <- readMappings(file = "Parata.emap.go.slim.txt")
xtroGO <- readMappings(file = "Xentro.emap.go.slim.txt")

fams=read.table("../Vert_Evt_OGr.txt",sep='\t')
names(fams)=c('fid','clg','clg_in','evt','sp','gid','gname','chrm','pos')
head(fams)

##### Gar #####

fams %>% filter(sp=='Lepocu') %>% group_by(fid) %>%
  mutate(has12=ifelse(any(grepl('1',evt)) & any(grepl('2',evt)),1,0)) %>%
  mutate(hasAB=ifelse(any(grepl('alpha',evt)) & any(grepl('beta',evt)),1,0)) %>%
  mutate(npar=n()) %>% ungroup %>% filter(npar<=10) -> gar

#gar %>% filter(npar>1 & npar<=8) %>% group_by(has12,hasAB) %>% head #tally

ggplot(gar,aes(x=npar))

#Locu.dp12=go_enrich(filter(gar,has12==1 | hasAB==1)$gid,locuGO)
#Locu_dup12=go_enrich(filter(gar,has12==1)$gid,locuGO)
Locu.dpAB=go_enrich_bp(filter(gar,has12==0 & hasAB==1)$gid,locuGO)
Locu.dp12=go_enrich_bp(filter(gar,has12==1 & hasAB==0)$gid,locuGO)
Locu.dp2B=go_enrich_bp(filter(gar,has12==1 & hasAB==1)$gid,locuGO)
Locu.dp1A=go_enrich_bp(filter(gar,has12==1)$gid,locuGO)

Locu.dpAB$lab<-'ABonly_Gar'
Locu.dp12$lab<-'12only_Gar'
Locu.dp2B$lab<-'BothDup_Gar'
Locu.dp1A$lab<-'Has12_Gar'

gar.go<-rbind(Locu.dpAB,Locu.dp12,Locu.dp2B,Locu.dp1A)
gar.go %>% mutate(Fisher=as.numeric(Fisher)) -> gar.go
gar.go %>% mutate(Fisher=as.numeric(Fisher)) %>%
  filter(Annotated>25)  %>% group_by(lab) %>% slice_min(n=40,order_by=Fisher) %>%
  ungroup() -> gar.min#%>% dplyr::select(Term) -> selTerms
#gar.go %>% mutate(Fisher=gsub(" ","",Fisher)) -> gar.go
gar.go %>% filter(Term %in% selTerms$Term) -> gar.slt

gar.slt %>% dplyr::select(Term,Fisher,lab) %>%
  pivot_wider(id_cols=Term,names_from=lab,values_from = Fisher, values_fill=1) %>%
  column_to_rownames(var="Term") -> test

pheatmap(-log10(test),color=magma(100))

##### Frog ####

fams %>% filter(sp=='Xentro') %>% group_by(fid) %>%
  mutate(has12=ifelse(any(grepl('1',evt)) & any(grepl('2',evt)),1,0)) %>%
  mutate(hasAB=ifelse(any(grepl('alpha',evt)) & any(grepl('beta',evt)),1,0)) %>%
  mutate(npar=n()) %>% filter(npar<=10) %>% ungroup -> frog

frog %>% filter(npar>1) %>% group_by(has12,hasAB) %>% tally

Xtro.dpAB=go_enrich_bp(filter(frog,has12==0 & hasAB==1)$gid,xtroGO)
Xtro.dp12=go_enrich_bp(filter(frog,has12==1 & hasAB==0)$gid,xtroGO)
Xtro.dp2B=go_enrich_bp(filter(frog,has12==1 & hasAB==1)$gid,xtroGO)
Xtro.dp1A=go_enrich_bp(filter(frog,has12==1)$gid,xtroGO)

Xtro.dpAB$lab<-'ABonly_Frog'
Xtro.dp12$lab<-'12only_Frog'
Xtro.dp2B$lab<-'BothDup_Frog'
Xtro.dp1A$lab<-'Has12_Frog'

frog.go<-rbind(Xtro.dpAB,Xtro.dp12,Xtro.dp2B,Xtro.dp1A)
frog.go %>% mutate(Fisher=as.numeric(Fisher)) -> frog.go
frog.go %>% mutate(Fisher=as.numeric(Fisher)) %>%
  filter(Annotated>25)  %>% group_by(lab) %>% slice_min(n=40,order_by=Fisher) %>%
  ungroup() %>% dplyr::select(Term) -> frog.sel
#gar.go %>% mutate(Fisher=gsub(" ","",Fisher)) -> gar.go
frog.go %>% filter(Term %in% frog.sel$Term) -> frog.slt

frog.slt %>% dplyr::select(Term,Fisher,lab) %>%
  pivot_wider(id_cols=Term,names_from=lab,values_from = Fisher, values_fill=1) %>%
  column_to_rownames(var="Term") -> frog.test

rbind(frog.slt,gar.slt) -> slt
slt %>% dplyr::select(Term,Fisher,lab) %>%
  pivot_wider(id_cols=Term,names_from=lab,values_from = Fisher, values_fill=1) %>%
  column_to_rownames(var="Term") -> test
pheatmap(-log10(test),color=magma(100))

##### Lamprey ######

fams %>% filter(sp=='Petmar') %>% group_by(fid) %>%
  mutate(has12=ifelse(any(grepl('1',evt)) & any(grepl('2',evt)),1,0)) %>%
  mutate(npar=n()) %>% filter(npar<=10) %>% ungroup() -> lamprey
filter(lamprey,has12==1 & npar==2) %>% tally
filter(lamprey,npar>2 & has12==0) %>% tally

Pmar.dp12=go_enrich_bp(filter(lamprey,has12==1)$gid,pmarGO)
Pmar.dp1o=go_enrich_bp(filter(lamprey,has12==1 & npar==2)$gid,pmarGO)
Pmar.dpCd=go_enrich_bp(filter(lamprey,npar>2 & has12==0)$gid,pmarGO)

Pmar.dp12$lab<-'12_Lamprey'
Pmar.dp1o$lab<-'12only_Lamprey '
Pmar.dpCd$lab<-'No12_Lamprey '

lamprey.go=rbind(Pmar.dp12,Pmar.dp1o,Pmar.dpCd)


lamprey %>% filter(npar>1) %>% group_by(has12) %>% ggplot(aes(x=npar)) + geom_histogram()
# & grep('2',any(evt))) %>%
gar %>% filter(has12==1 | hasAB==1) %>% nrow

###### Hagfish #####

fams %>% filter(sp=='Parata') %>% group_by(fid) %>%
  mutate(has12=ifelse(any(grepl('1',evt)) & any(grepl('2',evt)),1,0)) %>%
  mutate(npar=n()) %>% filter(npar<=10) %>% ungroup -> hagfish

Pata.dp12=go_enrich_bp(filter(hagfish,has12==1)$gid,pataGO)
Pata.dp1o=go_enrich_bp(filter(hagfish,has12==1 & npar==2)$gid,pataGO)
Pata.dpCd=go_enrich_bp(filter(hagfish,npar>2 & has12==0 &  npar<=8)$gid,pataGO)


Pata.dp12$lab<-'12_Hagfish'
Pata.dp1o$lab<-'12only_Hagfish'
Pata.dpCd$lab<-'No12_Hagfish'

hagfish.go=rbind(Pata.dp12,Pata.dp1o,Pata.dpCd)

##



all.go <- rbind(gar.go,frog.go,lamprey.go,hagfish.go)

all.go %>% mutate(Fisher=as.numeric(Fisher)) -> all.go

write.table(all.go,file="Vert2R_Go_enrich_long_wg.txt",sep='\t',quote=F,row.names=F)

all.go %>% mutate(ns=ifelse(Fisher < 0.01,1,0)) %>% 
  group_by(Term) %>% summarise(n_sig=sum(ns)) %>% filter(n_sig>1) -> multi_sig

all.go %>%  filter(Annotated>25)  %>% group_by(lab) %>%
  slice_min(n=40,order_by=Fisher) %>%
  ungroup() %>% dplyr::select(Term) -> selTerms
table(selTerms$Term %in% multi_sig$Term)
#gar.go %>% mutate(Fisher=gsub(" ","",Fisher)) -> gar.go
all.go %>% filter(Term %in% selTerms$Term) -> all.slt

all.slt %>% dplyr::select(Term,Fisher,lab) %>%
  pivot_wider(id_cols=Term,names_from=lab,values_from = Fisher, values_fill=1) %>%
  column_to_rownames(var="Term") -> all.mat

all.mat<-read.table("Vert2R_GO_enrich.txt",sep='\t',h=T)
all.mat %>% filter(include==1) %>% select(-include) %>% column_to_rownames(var='X')-> all.mat2
ghmp<-pheatmap(-log10(all.mat2),color=magma(100),breaks=seq(0,15,length=100))
rownames(all.mat[ghmp$tree_row[["order"]],])
write.table(all.mat[ghmp$tree_row[["order"]],],file="Vert2R_GO_enrich.txt",sep='\t',quote=F)


##### Get corresponding genes: 

genNfo<-read_tsv('~/Dropbox/Genomes/Hagfish/Hagfish_genome/gene_family_1/Eptata_genes_filt.txt')
gid2gnm<-genNfo %>% select(gid,gname) %>% column_to_rownames(var='gid')
gid2gnm<-as.vector(gid2gnm)
res<-go_enrich_bp(filter(hagfish,has12==1)$gid,pataGO)
gog=go_genes(filter(lamprey,has12==1)$gid,pmarGO)

selGenes=filter(hagfish,has12==1)$gid
res %>% mutate(gids=paste(sel[[GO.ID]],collapse=','))

lapply(results.tab$GO.ID, function(x) as.character(unlist(genesInTerm(object = topGOobject, whichGO = x))))
ggt<-as.character(unlist(genesInTerm(object = gog, whichGO = 'GO:0001755')))
intersect(ggt, filter(lamprey,has12==1)$gid)
sel=lapply(gog,function(x) x[x %in% selGenes] )
res %>% mutate(gids=paste(sel[[GO.ID]],collapse=','))
paste(as.vector(sel[res$GO.ID[1:6]]),collapse=',')


AnnotatedGenes = lapply(res$GO.ID, function(x) as.character(unlist(genesInTerm(object = gog, whichGO = x)))) # list containg genes annotated to significant GO terms
SignificantGenes = lapply(AnnotatedGenes, function(x) intersect(x, selGenes)) # where INT.GENES$V1 is your list of interesting genes
mrgSig<-lapply(SignificantGenes,function(x) paste(x,collapse=','))

##### Gene loss

#gfaNfo<-read_tsv("all_families_gain_loss.txt")
#gfaNfo %>% filter(grepl("Myxines",Losses)) %>% select(Gar_ids) %>% drop_na() %>% as.vector -> hag_los_gar #%>% str_split(Gar_ids, ",")
LossGarList<-read.table("gar_list.txt",h=F)
LossFrogList<-read.table("frog_list.txt",h=F)
LossPmarList<-read.table("Pmar_list.txt",h=F)

loss.gar=go_enrich_bp(as.vector(LossGarList$V1),locuGO)
loss.xtr=go_enrich_bp(LossFrogList$V1,xtroGO)
loss.pma=go_enrich_bp(LossPmarList$V1,pmarGO)
go_genes(LossPmarList$V1,pmarGO,'GO:0003341')
go_genes(LossGarList$V1,locuGO,'GO:0003341')

GO:0003341

write.table(loss.gar,file="Loss_Gar_GO_gn.tsv",quote=F,sep='\t',row.names=F)
write.table(loss.xtr,file="Loss_Frog_GO_gn.tsv",quote=F,sep='\t',row.names=F)
write.table(loss.pma,file="Loss_Lemprey_GO_gn.tsv",quote=F,sep='\t',row.names=F)


## Tweaking stuff... 

#geneList <- factor(as.integer(names(gene2GO) %in% selGenes))
#names(geneList) <- names(gene2GO)
#print(table(geneList))
#ontoList=list()
#ontologies=c("BP","CC","MF")

geneList <- factor(as.integer(names(pmarGO) %in% LossPmarList$V1))
names(geneList) <- names(pmarGO)
print(table(geneList))
#ontoList=list()
#ontologies=c("BP","CC","MF")

GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList,
              annot = annFUN.gene2GO, gene2GO = pmarGO)


loss.pma$genes
genes2 <- sapply(loss.pma$GO.ID, function(x)
{
  genes<-genesInTerm(GOdata, x)
  genes[[1]][genes[[1]] %in% as.vector(LossPmarList$V1)] # myGenes is the queried gene list
  })
loss.pma$genes<-unlist(lapply(genes2,paste, collapse=','))

sapply(genes2, paste, collapse = ",",USE.NAMES = TRUE) 
sapplt(paste(genes2[[3]],collapse=',')
goresults$genes <- sapply(goresults$GO.ID, function(x)
{
  genes<-genesInTerm(GOdata, x) 
  genes[[1]][genes[[1]] %in% SignificantGenes]
  
})


ggt<-as.character(unlist(genesInTerm(object = GOdata)))
cgt<-intersect(ggt,selGenes)

geneList <- factor(as.integer(names(gene2GO) %in% selGenes))
names(geneList) <- names(gene2GO)
print(table(geneList))
goa[loss.pma$GO.ID]
tgo<-

allRes$genes[which(allRes$topGO<0.05)] # print those only with p-value < 0.05  
  
  
GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList,
              annot = annFUN.gene2GO, gene2GO = gene2GO)
res <- runTest(GOdata, algorithm = "elim", statistic = "Fisher" )
resGenes <- GenTable(GOdata,Fisher = res, orderBy = "Fisher",topNodes=1000,numChar=1000)
resGenes %>% mutate(Fisher = gsub("< ","",Fisher)) -> resGenes


