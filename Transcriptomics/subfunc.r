library(ggplot2)
library(tidyverse)
library(ggvenn)
require(stringi)
require(stringdist)
require(ggpubr)
#ohno.Dr=read.table("Ohnologues/DANRE_ohnofam.txt")
#ohno.Mm=read.table("Ohnologues/MOUSE_ohnofam.txt")
#head(ohno.Pa)
#head(ohno.Bl)
setwd("/Users/fmarletaz/Dropbox/Genomes/Hagfish/Hagfish_analysis/tau/")

ohno.Bl=read.table("Bralan_ohnofams.txt",h=T)
ohno.Pa=read.table("Parata_ohnofams.txt",h=T)
ohno.Lo=read.table("Lepocu_ohnofams.txt",h=T)
ohno.Lj=read.table("Letjap_ohnofams.txt",h=T)



accAge=c('Vertebrata','Cyclostomata')
#levels(as.factor(xtroOH$age))
ohno.Pa %>% filter(dups>=1 & dups<=6 & age %in% accAge) -> fams.Pa #-> ohno_dups
ohno.Lj %>% filter(dups>=1 & dups<=6 & age %in% accAge) -> fams.Lj #-> ohno_dups

accAgeV=c('Vertebrata','Gnathostomata','Osteichthyes')
ohno.Lo %>% filter(dups>=1 & dups<=4 & age %in% accAgeV) -> fams.Lo #%>% summary()

accAgeB=c("Cephalochordata","Bralan","Ceph2" )
ohno.Bl %>% filter(dups>=1 & dups<=6 & age %in% accAgeB) -> fams.Bl
levels(as.factor(ohno.Bl$age))
head(fams.Bl)

fpkm.Pa=read.table('~/Dropbox/Genomes/Hagfish/Hagfish_genome/gene_family_1/Pata_organs_fpkms.tsv',h=T,row.names=1)
fpkm.Lo=read.table('~/Dropbox/Genomes/Hagfish/Hagfish_genome/gene_family_1/Locu_organs_fpkms.tsv',h=T,row.names=1)
fpkm.Lj=read.table('~/Dropbox/Genomes/Hagfish/Hagfish_genome/Data/Ljap_organs_fpkms.tsv',h=T,row.names=1)

rpkm.Bl=read.table('Bla_rpkm_sel.txt',h=T,row.names=1)
nrow(fpkm.Pa)
nrow(fpkm.Lo)
nrow(fpkm.Lj)
nrow(rpkm.Bl)

fpkm.Pa<-fpkm.Pa[rowSums(fpkm.Pa)!=0,]
fpkm.Lo<-fpkm.Lo[rowSums(fpkm.Lo)!=0,]
fpkm.Lj<-fpkm.Lj[rowSums(fpkm.Lj)!=0,]
rpkm.Bl<-rpkm.Bl[rowSums(rpkm.Bl)!=0,]


nrow(fpkm.Pa[rowSums(fpkm.Pa)!=0,])
nrow(fpkm.Lo[rowSums(fpkm.Lo)!=0,])
nrow(fpkm.Lj[rowSums(fpkm.Lj)!=0,])
nrow(rpkm.Bl[rowSums(rpkm.Bl)!=0,])

head(fpkm.Pa)
head(fpkm.Lo)
colnames(fpkm.Pa)
colnames(fpkm.Lo)

bin.Lo=ifelse(fpkm.Lo>5,1,0)
bin.Pa=ifelse(fpkm.Pa>5,1,0)
bin.Lj=ifelse(fpkm.Lj>5,1,0)
bin.Bl=data.frame(ifelse(rpkm.Bl>5,1,0))
nrow(bin.Lj)
nrow(bin.Lj[rowSums(bin.Lj)!=0,])
nrow(fpkm.Lj[rowSums(fpkm.Lj)!=0,])
head(bin.Lo)
head(bin.Pa)
colnames(bin.Lj)
#data.frame(bin.Lo) %>% select(brain,gills,heart,intestine,kidney,liver, muscle,ovary) -> bin.r.Lo
#data.frame(bin.Pa) %>% select(Brain,Gills,Heart,Intestine,Kidney,Liver,SkeletalMuscle,Ovaries) -> bin.r.Pa

data.frame(bin.Lo) %>% select(brain,gills,intestine,liver, muscle,ovary) -> bin.r.Lo
data.frame(bin.Pa) %>% select(Brain,Gills,Intestine,Liver,SkeletalMuscle,Ovaries) -> bin.r.Pa
data.frame(bin.Lj) %>% select(Lj_brain,Lj_gill,LJ_Intestine,Lj_liver,JL_muscle,JL_ovary) -> bin.r.Lj
bin.Bl %>% select(NeuralTube,Gills_b,Gut_b,Hepatic,Muscle_b,FemGonads_b) -> bin.r.Bl



#fams.Pa %>% mutate(cluster=wgcna.Pa[match(gid,wgcna.Pa$gid),1]) 

fams.Pa %>% filter(n_tissue!=0) %>% group_by(fid,Patt) %>% mutate(scat=ifelse(min(n)<=2 & max(n)>=7,) %>% tally

bin.r.Lo %>% mutate(n_tissue=rowSums(bin.r.Lo)) %>% unite('Patt',brain:ovary,remove=T,sep='') %>% rownames_to_column(var='gid')-> exp.Lo
bin.r.Pa %>% mutate(n_tissue=rowSums(bin.r.Pa)) %>% unite('Patt',Brain:Ovaries,remove=T,sep='') %>% rownames_to_column(var='gid')-> exp.Pa
bin.r.Lj %>% mutate(n_tissue=rowSums(bin.r.Lj)) %>% unite('Patt',Lj_brain:JL_ovary,remove=T,sep='') %>% rownames_to_column(var='gid')-> exp.Lj
bin.r.Bl  %>% mutate(n_tissue=rowSums(bin.r.Bl)) %>% unite('Patt',NeuralTube:FemGonads_b,remove=T,sep='') %>% rownames_to_column(var='gid')-> exp.Bl                                                          

left_join(fams.Pa,exp.Pa) -> fams.Pa
left_join(fams.Lo,exp.Lo) -> fams.Lo
left_join(fams.Lj,exp.Lj) -> fams.Lj
left_join(fams.Bl,exp.Bl) -> fams.Bl

#fams.Pa %>% mutate()                                                              
#case_when(min(n)<=2 & max(n)>=7 ~ 'special', min(n)>=6len(distinct(Patt)) > 1 ~ 'some_sub' ~clg,is.na(clg) ~ clgb,is.na(clgb) ~ clg,TRUE ~clg))

fams.Lo %>% select(fid,Patt) %>% rename(Pat_Lo=Patt) -> pats.Lo
fams.Bl %>% select(fid,Patt) %>% rename(Pat_Bl=Patt) -> pats.Bl
left_join(pats.Lo,pats.Bl) -> pats.BlLo
pats.BlLo  %>% group_by(fid) %>% expand(Pat_Lo,Pat_Bl,.name_repair="minimal") -> pats.exp.BlLo
pats.exp.BlLo %>% drop_na() %>% mutate(pdst=stringdist(Pat_Lo,Pat_Bl)*stri_compare(Pat_Lo,Pat_Bl)) %>%
  filter(pdst > -7)-> pats.cnt.BlLo

BlLo<-ggplot(pats.cnt.BlLo, aes(x=pdst))+
  geom_histogram(color="darkblue", fill="lightblue",binwidth=1)+ggtitle('Lo <  > Bl')+
  theme_minimal()+theme(plot.title = element_text(hjust = 0.5))
BlLo
fams.Pa %>% select(fid,Patt) %>% rename(Pat_Pa=Patt) -> pats.Pa
#fams.Bl %>% select(fid,Patt) %>% rename(Pat_Bl=Patt) -> pats.Bl
left_join(pats.Pa,pats.Bl) -> pats.BlPa
pats.BlPa  %>% group_by(fid) %>% expand(Pat_Pa,Pat_Bl,.name_repair="minimal") -> pats.exp.BlPa
pats.exp.BlPa %>% drop_na() %>% mutate(pdst=stringdist(Pat_Pa,Pat_Bl)*stri_compare(Pat_Pa,Pat_Bl)) %>%
  filter(pdst > -7)-> pats.cnt.BlPa

BlPa<-ggplot(pats.cnt.BlPa, aes(x=pdst))+
  geom_histogram(color="darkblue", fill="lightblue",binwidth=1)+ggtitle('Pa <  > Bl')+
  theme_minimal()+theme(plot.title = element_text(hjust = 0.5))

fams.Lj %>% select(fid,Patt) %>% rename(Pat_Lj=Patt) -> pats.Lj
#fams.Bl %>% select(fid,Patt) %>% rename(Pat_Bl=Patt) -> pats.Bl
left_join(pats.Lj,pats.Bl) -> pats.BlLj
pats.BlLj  %>% group_by(fid) %>% expand(Pat_Lj,Pat_Bl,.name_repair="minimal") -> pats.exp.BlLj
pats.exp.BlLj %>% drop_na() %>% mutate(pdst=stringdist(Pat_Lj,Pat_Bl)*stri_compare(Pat_Lj,Pat_Bl)) %>%
  filter(pdst > -7)-> pats.cnt.BlLj

BlLj<-ggplot(pats.cnt.BlLj, aes(x=pdst))+
  geom_histogram(color="darkblue", fill="lightblue",binwidth=1)+ggtitle('Lj < > Bl')+
  theme_minimal()+theme(plot.title = element_text(hjust = 0.5))

#fams.Pa %>% select(fid,Patt) %>% rename(Pat_Pa=Patt) -> pats.Pa
#fams.Bl %>% select(fid,Patt) %>% rename(Pat_Bl=Patt) -> pats.Bl
left_join(pats.Pa,pats.Lo) -> pats.LoPa
pats.LoPa  %>% group_by(fid) %>% expand(Pat_Pa,Pat_Lo,.name_repair="minimal") -> pats.exp.LoPa
pats.exp.LoPa %>% drop_na() %>% mutate(pdst=stringdist(Pat_Pa,Pat_Lo)*stri_compare(Pat_Pa,Pat_Lo)) %>% filter(pdst > -7)-> pats.cnt.LoPa
LoPa<-ggplot(pats.cnt.LoPa, aes(x=pdst))+
  geom_histogram(color="darkblue", fill="lightblue",binwidth=1)+ggtitle('Pa < > Lo')+
  theme_minimal()+theme(plot.title = element_text(hjust = 0.5))



left_join(pats.Lj,pats.Lo) -> pats.LoLj
pats.LoLj  %>% group_by(fid) %>% expand(Pat_Lj,Pat_Lo,.name_repair="minimal") -> pats.exp.LoLj
pats.exp.LoLj %>% drop_na() %>% mutate(pdst=stringdist(Pat_Lj,Pat_Lo)*stri_compare(Pat_Lj,Pat_Lo)) %>% 
  filter(pdst > -7)-> pats.cnt.LoLj

LoLj<-ggplot(pats.cnt.LoLj, aes(x=pdst))+
  geom_histogram(color="darkblue", fill="lightblue",binwidth=1)+ggtitle('Lj <  > Lo')+
  theme_minimal()+theme(plot.title = element_text(hjust = 0.5))

left_join(pats.Lj,pats.Pa) -> pats.LjPa
pats.LjPa  %>% group_by(fid) %>% expand(Pat_Lj,Pat_Pa,.name_repair="minimal") -> pats.exp.LjPa
pats.exp.LjPa %>% drop_na() %>% mutate(pdst=stringdist(Pat_Lj,Pat_Pa)*stri_compare(Pat_Lj,Pat_Pa)) %>% 
  filter(pdst > -7)-> pats.cnt.LjPa

LjPa<-ggplot(pats.cnt.LjPa, aes(x=pdst))+
  geom_histogram(color="darkblue", fill="lightblue",binwidth=1)+ggtitle('Lj <  > Pa')+
  theme_minimal()+theme(plot.title = element_text(hjust = 0.5))

ggarrange(BlLo,BlLj,BlPa,LoPa,LoLj,LjPa)

head(pats.cnt.BlPa)

pats.cnt.BlPo %>% filter(pdst< -10)

fids=list(`Hagfish`=fams.Pa$fid,`Gar`=fams.Lo$fid)
bin.Pa[rownames(bin.Pa) %in% filter(fams.Pa,fid=='OG_1016')$gid,]
ggvenn(fids, c("Hagfish", "Gar"))
table(fams.Pa$fid %in% fams.Pa$fid)

wgcna.Pa=read.table("Pata_wgnca_genes.txt")
wgcna.Lo=read.table("Locu_wgcna_genes.txt")
names(wgcna.Lo)=c('clust','gid')
names(wgcna.Pa)=c('clust','gid')



fams.Pa %>% mutate(cluster=wgcna.Pa[match(gid,wgcna.Pa$gid),1])  %>%
  group_by(fid) %>% mutate(n_clust=n_distinct(cluster))-> fams.Pa

fams.Pa %>% ungroup() %>% group_by(n_clust) %>% tally()%>% mutate(per=100*n/sum(n))

fams.Lo %>% mutate(cluster=wgcna.Lo[match(gid,wgcna.Lo$gid),1])  %>%
  group_by(fid) %>% mutate(n_clust=n_distinct(cluster)) -> fams.Lo
fams.Lo %>% ungroup() %>% group_by(n_clust)  %>% tally()%>% mutate(per=100*n/sum(n))

