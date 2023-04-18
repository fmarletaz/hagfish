library(ggplot2)
library(tidyverse)
library(ggpubr)


#################

calcTau <- function(dataset,sp) {
  normalised=t(apply(log10(dataset+1),1, function(x) x/max(x)))
  tau=apply(normalised,1,function(x) sum(1-x)/length(x))
  dtau=data.frame(unlist(tau))
  dtau$species=sp
  colnames(dtau)=c('tau','species')
  dtau=dtau[!is.na(dtau$tau),]
  return(dtau)
}

median.quartile <- function(x){
  out <- quantile(x, probs = c(0.25,0.5,0.75))
  names(out) <- c("ymin","y","ymax")
  return(out) 
}

#################


setwd('~/Dropbox/Genomes/Hagfish/Hagfish_analysis/tau/')

# Load RPKMs
rpkm.Bl=read.table('Bla_forWGCNA-orth-v2-11_09_16-RUN8.tab',h=T,row.names=1)
rpkm.Bl=select(rpkm.Bl,-NAME)
head(rpkm.Bl)
colnames(rpkm.Bl)
#rpkm.Dr=read.table('manu_cRPKMs/Dre_forWGCNA-orth-v2-11_09_16-RUN8.tsv',h=T,row.names=1,sep='\t')
#rpkm.Dr=select(rpkm.Dr,-NAME)
#colnames(rpkm.Dr)
#rpkm.Mm=read.table('manu_cRPKMs/Expr-Mmu72-ALL_GENES_AND_COUNTs.tsv',h=T,row.names=1,sep='\t')
#rpkm.Mm=select(rpkm.Mm,contains(".R"))
#rpkm.Xt=read.table('manu_cRPKMs/Expr-Xtr47-ALL_GENES_AND_COUNTs.tsv',h=T,row.names=1,sep='\t')
#rpkm.Xt=select(rpkm.Xt,-NAME)
#ct.Pa=read.table("~/Dropbox/Genomes/Hagfish/Transcriptomics/Hagfish_organs.counts.tsv",h=T,row.names=1)


fpkm.Pa=read.table('~/Dropbox/Genomes/Hagfish/Hagfish_genome/gene_family_1/Pata_organs_fpkms.tsv',h=T,row.names=1)
fpkm.Lo=read.table('~/Dropbox/Genomes/Hagfish/Hagfish_genome/gene_family_1/Locu_organs_fpkms.tsv',h=T,row.names=1)
fpkm.Lj=read.table('~/Dropbox/Genomes/Hagfish/Hagfish_genome/Data/Ljap_organs_fpkms.tsv',h=T,row.names=1)

# generate a constistent set of organs and average

rpkm.Bl %>% filter(rowSums(.)!=0)

#colnames(rpkm.Bl)
#fpkm.Pa %>% select()
#rpkm.Bl$Hepatic=rowMeans(select(rpkm.Bl, starts_with('Hepatic')), na.rm = TRUE)
#rpkm.Bl$NeuralTube=rowMeans(select(rpkm.Bl, starts_with('NeuralTube_')), na.rm = TRUE)
#rpkm.Bl=select(rpkm.Bl, -starts_with('Hepatic_'))
#rpkm.Bl=select(rpkm.Bl, -starts_with('NeuralTube_'))

rpkm.Bl %>% select(Epidermis_b,FemGonads_b,Gills_b,Gut_b,Muscle_b,Hepatic_b,NeuralTube_b,Cirri_b) %>% filter(rowSums(.)!=0) -> rpkm.Bl.f
fpkm.Pa %>% select(Brain,Gills,Heart,Intestine,Kidney,Liver,Ovaries,Skin,SkeletalMuscle) %>% filter(rowSums(.)!=0)-> fpkm.Pa.f
fpkm.Lo %>% select(brain,gills,heart,intestine,kidney,liver,muscle,ovary,bones) %>% filter(rowSums(.)!=0) -> fpkm.Lo.f
fpkm.Lj %>% select(Lj_brain,Lj_gill,Lj_heart,LJ_Intestine,Lj_kidney,Lj_liver,JL_muscle,JL_ovary) %>% filter(rowSums(.)!=0) -> fpkm.Lj.f
nrow(fpkm.Lj.f)
nrow(fpkm.Lj[rowSums(fpkm.Lj)!=0,])
head(fpkm.Lj)
#rpkm.Dr$Muscle=rowMeans(select(rpkm.Dr, starts_with('Muscle_')), na.rm = TRUE)
#rpkm.Dr$Intestine=rowMeans(select(rpkm.Dr, starts_with('Intestine_')), na.rm = TRUE)
#rpkm.Dr$Liver=rowMeans(select(rpkm.Dr, starts_with('Liver_')), na.rm = TRUE)
#rpkm.Dr$Pancreas=rowMeans(select(rpkm.Dr, starts_with('Pancreas_')), na.rm = TRUE)
#rpkm.Dr$Kidney=rowMeans(select(rpkm.Dr, starts_with('Kidney_')), na.rm = TRUE)
#rpkm.Dr$Gill=rowMeans(select(rpkm.Dr, starts_with('Gill_')), na.rm = TRUE)
#rpkm.Dr$Brain=rowMeans(select(rpkm.Dr, starts_with('Brain_')), na.rm = TRUE)
#rpkm.Dr$Testis=rowMeans(select(rpkm.Dr, starts_with('Testis_')), na.rm = TRUE)
#rpkm.Dr=select(rpkm.Dr, -starts_with('Muscle_'),-starts_with('Intestine_'),-starts_with('Liver_'),-starts_with('Pancreas_'),-starts_with('Kidney_'),-starts_with('Gill_'),-starts_with('Brain_'),-starts_with('Testis_'))
#colnames(rpkm.Dr)

head(fpkm.Pa[order(fpkm.Pa$SlimeGland,decreasing =T),])

#colnames(rpkm.Mm)
#rpkm.Mm=select(rpkm.Mm,Embr_Brain_E14.5.R,Embr_Limb_E14.5.R,Embr_Liver_E14.5.R,Oocyte.R,Embr_8C.R,Heart_a.R,Heart_b.R,Muscle_a.R,Muscle_b.R,Intestine.R,Stomach.R,Colon.R,Liver_a.R,Liver_b.R,Kidney_a.R,Kidney_b.R,Lung.R,Retina_Eye.R,Whole_Brain_a.R,Whole_Brain_b.R,Testis_a.R,Testis_b.R,Ovary.R,Spleen.R,Thymus.R,Placenta_a.R,Placenta_b.R)
#rpkm.Mm$Muscle=rowMeans(select(rpkm.Mm, starts_with('Muscle_')))
#rpkm.Mm$Heart=rowMeans(select(rpkm.Mm, starts_with('Heart_')))
#rpkm.Mm$Liver=rowMeans(select(rpkm.Mm, starts_with('Liver_')))
#rpkm.Mm$Kidney=rowMeans(select(rpkm.Mm, starts_with('Kidney_')))
rpkm.Mm$Brain=rowMeans(select(rpkm.Mm, starts_with('Whole_Brain_')))
rpkm.Mm$Testis=rowMeans(select(rpkm.Mm, starts_with('Testis_')))
rpkm.Mm$Embr14.5=rowMeans(select(rpkm.Mm, Embr_Brain_E14.5.R,Embr_Limb_E14.5.R,Embr_Liver_E14.5.R))
rpkm.Mm$Placenta=rowMeans(select(rpkm.Mm, starts_with('Placenta_')))

rpkm.Mm %>% select(-starts_with('Muscle_')) %>% 
  select(-starts_with('Heart_')) %>%
  select(-starts_with('Liver_')) %>%
  select(-starts_with('Whole_Brain_')) %>%
  select(-starts_with('Gill_')) %>%
  select(-starts_with('Testis_')) %>%
  select(-starts_with('Placenta_')) %>%
  select(-starts_with('Kidney_')) %>%
  select(-Embr_Brain_E14.5.R,-Embr_Limb_E14.5.R,-Embr_Liver_E14.5.R) -> rpkm.Mm
colnames(rpkm.Mm)

#head(rpkm.Bl %>% tbl_df() %>% rownames_to_column('gene') %>%  gather(organ,rpkm,-gene) %>% colsplit(organ,))
#colnames(select(rpkm.Mm, -starts_with('Muscle_')))
#group_by(gene,organ) %>% summarize(rpkm = mean(rpkm, na.rm = TRUE)))
#head(rpkm.Bl %>% tbl_df() %>% add_rownames('gene') %>%  gather(organ,rpkm,-gene)  )
#factor_key=TRUE)separate(col, into, sep = " ", remove = TRUE, convert = FALSE)

ohno.Bl=read.table("Bralan_ohnofams.txt",h=T)
ohno.Pa=read.table("Parata_ohnofams.txt",h=T)
ohno.Lo=read.table("Lepocu_ohnofams.txt",h=T)
ohno.Lj=read.table("Letjap_ohnofams.txt",h=T)
#ohno.Dr=read.table("Ohnologues/DANRE_ohnofam.txt")
#ohno.Mm=read.table("Ohnologues/MOUSE_ohnofam.txt")
#head(ohno.Pa)
#head(ohno.Bl)

ohno.Lo %>% group_by(fid) %>% select(gid) %>% tally
ohno.Lo %>% group_by(fid) %>% select(gid) %>% tally
ggplot(ohno.Lo,aes(x=dups))+
  geom_histogram(color="darkblue", fill="lightblue") 
accAge=c('Vertebrata','Cyclostomata')
#levels(as.factor(xtroOH$age))
ohno.Pa %>% filter(dups>=1 & dups<=8 & age %in% accAge) -> fams.Pa #-> ohno_dups
ohno.Lj %>% filter(dups>=1 & dups<=8 & age %in% accAge) -> fams.Lj #-> ohno_dups

accAgeV=c('Vertebrata','Gnathostomata','Osteichthyes')
ohno.Lo %>% filter(dups>=1 & dups<=8 & age %in% accAgeV) -> fams.Lo #%>% summary()

fams.Pa$sp<-'Pa'
fams.Lo$sp<-'Lo'
fams.Lj$sp<-'Lj'
fams<-rbind(fams.Pa,fams.Lo,fams.Lj)
fams %>%group_by(fid,sp) %>% tally -> famct
famct %>% group_by(sp,n) %>% tally
fams$sp=factor(fams$sp,levels=c('Lo','Lj','Pa'))
ggplot(fams,aes(x=dups,fill=sp))+geom_bar(stats='identity',position='dodge') +theme_light()
hist(fams.Lo$dups)
ggplot(fams.Pa,aes(x=dups))+geom_histogram() 
ggplot(fams.Lj,aes(x=dups))+geom_histogram() 


accAgeB=c("Cephalochordata","Bralan","Ceph2" )
ohno.Bl %>% filter(dups>=1 & dups<=4 & age %in% accAgeB) -> fams.Bl
levels(as.factor(ohno.Bl$age))
head(fams.Pa)


tau.Bl=calcTau(rpkm.Bl.f,'Bla')
tau.Pa=calcTau(fpkm.Pa.f,'Eat')
tau.Lo=calcTau(fpkm.Lo.f,'Loc')
tau.Lj=calcTau(fpkm.Lj.f,'Lja')
#tau.Dr=calcTau(rpkm.Dr,'Dre')
#tau.Mm=calcTau(rpkm.Mm,'Mmu')
#tau.Xt=calcTau(rpkm.Xt,'Xtr')
#head(tau.Bl)
#table(is.na(tau.Bl$tau))

ohno.Bl %>% filter(fid %in% as.vector(ohno_dups$fid)) -> ohno.Bl.d

tau.ohno.Bl=tau.Bl[as.vector(fams.Bl$gid),]
tau.ohno.Bl$fid=fams.Bl$fid
tau.ohno.Bl %>% na.omit() -> tau.ohno.Bl

tau.ohno.Pa=tau.Pa[as.vector(fams.Pa$gid),]
tau.ohno.Pa$fid=fams.Pa$fid
tau.ohno.Pa %>% na.omit() -> tau.ohno.Pa

tau.ohno.Lo=tau.Lo[as.vector(fams.Lo$gid),]
tau.ohno.Lo$fid=fams.Lo$fid
tau.ohno.Lo %>% na.omit() -> tau.ohno.Lo

tau.ohno.Lj=tau.Lj[as.vector(fams.Lj$gid),]
tau.ohno.Lj$fid=fams.Lj$fid
tau.ohno.Lj %>% na.omit() -> tau.ohno.Lj

tau.ohno=rbind(tau.ohno.Bl,tau.ohno.Pa,tau.ohno.Lo,tau.ohno.Lj)

write.table(tau.ohno,file="Tau_ohno_sp.tsv",sep='\t',quote=F)

head(tau.ohno.Pa)
tau.all=rbind(tau.Bl,tau.Pa)

head(rpkm.Bl)
head(rpkm.Dr)
head(rpkm.Mm)
head(rpkm.Xt)
pdf('tau_re.pdf')

ggplot(tau.all, aes(x=tau, fill=species)) +   geom_density(alpha=0.3)+ ggtitle("All genes")
ggplot(tau.all, aes(tau,colour=species)) + stat_ecdf()+ ggtitle("All genes")
ggplot(tau.ohno, aes(tau,colour=species)) + stat_ecdf()+ ggtitle("Ohnologues")

ggplot(tau.all, aes(x=species, y=tau)) +  geom_violin(aes(fill=species),scale = "width")+stat_summary(fun.y="median.quartile", geom="line")+stat_summary(fun.y="mean", geom="point")+theme_bw()+ ggtitle("All genes")

tau.ohno.Bl=tau.Bl[as.vector(ohno.Bl$V1),]
tau.ohno.Dr=tau.Dr[as.vector(ohno.Dr$V1),]
tau.ohno.Mm=tau.Mm[as.vector(ohno.Mm$V1),]
tau.ohno.Xt=tau.Xt[as.vector(ohno.Xt$V1),]

tau.ohno.Bl$hog=ohno.Bl$V2
tau.ohno.Dr$hog=ohno.Dr$V2
tau.ohno.Mm$hog=ohno.Mm$V2
tau.ohno.Xt$hog=ohno.Xt$V2
tau.ohno=rbind(tau.ohno.Bl,tau.ohno.Dr,tau.ohno.Mm,tau.ohno.Xt)
tau.ohno=tau.ohno[!is.na(tau.ohno$species),]

ggplot(tau.all, aes(x=species, y=tau)) +  geom_violin(aes(fill=species),scale = "width")+stat_summary(fun.y="median.quartile", geom="line")+stat_summary(fun.y="mean", geom="point")+theme_bw()+ ggtitle("All genes")

head(tau.ohno)
ggplot(tau.ohno, aes(x=species, y=tau)) +  geom_violin(aes(fill=species))+stat_summary(fun.y="median.quartile", geom="line")+stat_summary(fun.y="mean", geom="point")+theme_bw()+ ggtitle("Ohnologues")

tau.ohno.max=aggregate(tau.ohno$tau,by=list(fid=tau.ohno$fid,species=tau.ohno$species),FUN=max,na.rm=TRUE)
tau.ohno.med=aggregate(tau.ohno$tau,by=list(fid=tau.ohno$fid,species=tau.ohno$species),FUN=median,na.rm=TRUE)

tau.ohno.min=aggregate(tau.ohno$tau,by=list(fid=tau.ohno$fid,species=tau.ohno$species),FUN=min,na.rm=TRUE)

tau.ohno.max.pw=spread(tau.ohno.max,species,x)
tau.ohno.max.pw$stat<-'max'
tau.ohno.min.pw=spread(tau.ohno.min,species,x)
tau.ohno.min.pw$stat<-'min'
tau.ohno.med.pw=spread(tau.ohno.med,species,x)
tau.ohno.med.pw$stat<-'med'



tau.ohno.max.pw<-tau.ohno.max.pw[complete.cases(tau.ohno.max.pw),]
pdf('tau_max.pdf')

ggplot(tau.ohno.med.pw,aes(x=Bla,y=Lja))+theme_bw()+
  stat_density2d(aes(fill = ..density..), geom = "tile", contour = FALSE, n = 200)+
  scale_fill_continuous(low = "white", high = "deepskyblue3")+geom_point(col='dimgrey',size=.5)+
  geom_abline(intercept = 0, slope = 1, col='firebrick4')

dev.off()

bl<-ggplot(tau.ohno.max.pw,aes(x=Bla,y=Loc))+theme_bw()+
  stat_density2d(aes(fill = ..density..), geom = "tile", contour = FALSE, n = 200)+
  scale_fill_continuous(low = "white", high = "deepskyblue3")+geom_point(col='dimgrey',size=.5)+
  geom_abline(intercept = 0, slope = 1, col='firebrick4')

bp<-ggplot(tau.ohno.max.pw,aes(x=Bla,y=Eat))+theme_bw()+
  stat_density2d(aes(fill = ..density..), geom = "tile", contour = FALSE, n = 200)+
  scale_fill_continuous(low = "white", high = "deepskyblue3")+geom_point(col='dimgrey',size=.5)+
  geom_abline(intercept = 0, slope = 1, col='firebrick4')

bj<-ggplot(tau.ohno.max.pw,aes(x=Bla,y=Lja))+theme_bw()+
  stat_density2d(aes(fill = ..density..), geom = "tile", contour = FALSE, n = 200)+
  scale_fill_continuous(low = "white", high = "deepskyblue3")+geom_point(col='dimgrey',size=.5)+
  geom_abline(intercept = 0, slope = 1, col='firebrick4')

lp<-ggplot(tau.ohno.max.pw,aes(x=Loc,y=Eat))+theme_bw()+
  stat_density2d(aes(fill = ..density..), geom = "tile", contour = FALSE, n = 200)+
  scale_fill_continuous(low = "white", high = "deepskyblue3")+geom_point(col='dimgrey',size=.5)+
  geom_abline(intercept = 0, slope = 1, col='firebrick4')

lj<-ggplot(tau.ohno.max.pw,aes(x=Loc,y=Lja))+theme_bw()+
  stat_density2d(aes(fill = ..density..), geom = "tile", contour = FALSE, n = 200)+
  scale_fill_continuous(low = "white", high = "deepskyblue3")+geom_point(col='dimgrey',size=.5)+
  geom_abline(intercept = 0, slope = 1, col='firebrick4')

le<-ggplot(tau.ohno.max.pw,aes(x=Lja,y=Eat))+theme_bw()+
  stat_density2d(aes(fill = ..density..), geom = "tile", contour = FALSE, n = 200)+
  scale_fill_continuous(low = "white", high = "deepskyblue3")+geom_point(col='dimgrey',size=.5)+
  geom_abline(intercept = 0, slope = 1, col='firebrick4')

ggarrange(bl,bp,bj,lp,lj,le)

dev.off()
ggplot(tau.ohno.max.pw,aes(x=Bla,y=Xtr,col=stat))+geom_point(color='coral2')+theme_bw()+geom_abline(intercept = 0, slope = 1)
ggplot(tau.ohno.max.pw,aes(x=Bla,y=Mmu,col=stat))+geom_point(color='coral2')+theme_bw()+geom_abline(intercept = 0, slope = 1)
ggplot(tau.ohno.max.pw,aes(x=Dre,y=Mmu,col=stat))+geom_point(color='coral2')+theme_bw()+geom_abline(intercept = 0, slope = 1)
ggplot(tau.ohno.max.pw,aes(x=Xtr,y=Mmu,col=stat))+geom_point(color='coral2')+theme_bw()+geom_abline(intercept = 0, slope = 1)
ggplot(tau.ohno.max.pw,aes(x=Xtr,y=Dre,col=stat))+geom_point(color='coral2')+theme_bw()+geom_abline(intercept = 0, slope = 1)
dev.off()



library(plyr)
# Categorise
tau.all$cat=ifelse(tau.all$tau<0.4,'low',ifelse(tau.all$tau>0.6,'high','medium'))
tau.all$cat=factor(tau.all$cat,levels=c('high','medium','low'))
tau.all=filter(tau.all,!is.na(tau))
table(is.na(tau.all$tau))

tau.ohno$cat=ifelse(tau.ohno$tau<0.4,'low',ifelse(tau.ohno$tau>0.6,'high','medium'))
tau.ohno$cat=factor(tau.ohno$cat,levels=c('high','medium','low'))


tau.all %>%
  group_by(species,cat) %>%
  tally %>%
  group_by(species) %>%
  mutate(per = n / sum(n)*100) -> ct.all

tau.ohno %>%
  group_by(species,cat) %>%
  tally %>%
  na.omit() %>%
  group_by(species) %>%
  mutate(per = n / sum(n)*100) -> ct.ohno
pdf('categorised.pdf')
ggplot(ct.ohno,aes(x=species,y=per,fill=cat))+geom_bar(stat="identity")+theme_bw()+ ggtitle("Ohnologues")
ggplot(ct.all,aes(x=species,y=per,fill=cat))+geom_bar(stat="identity")+theme_bw()+ ggtitle("All genes")

ggplot(tau.ohno[complete.cases(tau.ohno),],aes(species))+ geom_bar(aes(fill = cat))+theme_bw()+ ggtitle("Ohnologues")
ggplot(tau.all[complete.cases(tau.all),],aes(species))+ geom_bar(aes(fill = cat))+theme_bw()+ ggtitle("All genes")
dev.off()




