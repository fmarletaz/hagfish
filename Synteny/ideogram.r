require(RIdeogram)
require(tidyverse)
library("wesanderson")


`%notin%` <- Negate(`%in%`)

setwd("~/Dropbox/Genomes/Synteny/synt-forge/Hagfish-redux/")
bf_colors=data.frame(row.names = c('A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q'),
                  col=c("1F78B4","33A02C","E31A1C","FFFF99","B15928","D95F02","7570B3","E7298A","66A61E","6A3D9A","E6AB02",
                        "A6761D","8DD3C7","666666","FF7F00","BEBADA","ffd92f"))
#bf_colors=data.frame(row.names=c('A','B1','B2','B3','C1','C2','D','E','F','G','H','I','J1','J2','K','L','M','N','O1','O2','P','Q'),
#                     col=c("1F78B4","1B9E77","33A02C","B2DF8A","FB9A99","E31A1C","FFFF99","B15928","D95F02","7570B3","E7298A",
#                           "66A61E","CAB2D6","6A3D9A","E6AB02","A6761D","8DD3C7","666666","FF7F00","FDBF6F","BEBADA","ffd92f"))

# EXAMPLE

data(karyotype_dual_comparison, package="RIdeogram")
head(karyotype_dual_comparison)
data(synteny_dual_comparison, package="RIdeogram")
head(synteny_dual_comparison)

##### Hagfish VS Lamprey #####

lg=read.table("Pata_Pmar_syntExp.txt",h=T)

lg %>% filter(s1chr %notin% c('chrm72','chrm71','S68_unloc.1','chrm68','chrm66')) %>%
  mutate(clg=gsub(pattern = "[1-3]", replacement = "", clg)) -> lg 
#
lg %>% select(s2chr,s1chr,clg) %>% filter(!is.na(clg)) %>% 
  group_by(s2chr,s1chr,clg) %>% tally %>% filter(n>5) %>% 
  arrange(clg,s1chr) -> lgct

lg %>% select(s2chr,s1chr,clg) %>% 
  group_by(s2chr,s1chr,clg) %>% tally %>% filter(s1chr=='chrm62')

?write.table()
write.table(lgct,file='Pata_Pmar_clg.txt',sep='\t',quote=F, row.names = F)




head(lg)
lg %>%  group_by(s1chr) %>%   summarise(End=max(s1gp)) %>% arrange(-End) %>% rename(Chr=s1chr) %>% 
  mutate(Start=1,species='Petromyzon',fill='Petromyzon',size=12,color=252525) %>% relocate(Start,.after=Chr)-> chr1

lg  %>% group_by(s2chr) %>%   summarise(End=max(s2gp)) %>% arrange(-End) %>% rename(Chr=s2chr) %>% 
  mutate(Start=1,species='Eptatretus',fill='Eptatretus',size=12,color=252525) %>% relocate(Start,.after=Chr)-> chr2

#unlist(pata_pmar_chOrd[2])

chr1_ord=  c('chrm16','chrm48','chrm41','chrm28','chrm46','chrm36','chrm42','chrm70','chrm56','chrm18','chrm29','chrm4','chrm34','chrm65','chrm40','chrm32','chrm49','chrm5','chrm61','chrm51','chrm11','chrm3','chrm74','chrm38','chrm39','chrm10','chrm2','chrm57','chrm26','chrm54','chrm6','chrm25','chrm45','chrm63','chrm64','chrm52','chrm17','chrm14','chrm53','chrm12','chrm9','chrm1','chrm35','chrm7','chrm19','chrm37','chrm58','chrm44','chrm31','chrm60','chrm69','chrm62','chrm67','chrm30','chrm33','chrm13','chrm20','chrm27','chrm8','chrm22','chrm24','chrm23','chrm47','chrm43','chrm21','chrm55','chrm15','chrm50','chrm59','chrm66','chrm68','S68_unloc.1','chrm71','chrm72')
chr2_ord= c('chr3','chr13','chr5','chr16','chr2','chr4','chr8','chr9','chr10','chr12','chr7','chr14','chr1','chr11','chr6','chr17','chr15')

lg %>% filter(s1chr=='chrm72')

chr1=chr1[order(match(chr1$Chr,chr1_ord)),]
chr2=chr2[order(match(chr2$Chr,chr2_ord)),]

rbind(chr1,chr2) -> chroms


lg$Species_1=match(lg$s1chr,chr1$Chr)
lg$Species_2=match(lg$s2chr,chr2$Chr)

lg %>% filter(scol=='S') %>% select(Species_1,s1gp,Species_2,s2gp,clg) %>% rename(Start_1=s1gp,Start_2=s2gp) %>%
  mutate(End_1=Start_1+1000,End_2=Start_2+1000,fill=bf_colors[match(clg,rownames(bf_colors),),]) %>%
  relocate(End_1,.after=Start_1) %>% relocate(End_2,.after=Start_2) %>% select(-clg) -> lgt


ideogram(karyotype = chroms, synteny = lgt)

##### Leri vs Locu #####

lg=read.table("Leri-Locu_syntr.txt",h=T)
head(lg)
lg %>% filter(s1chr %notin% c('LER50','LER45','LER44')) %>%
  filter(s2chr %notin% c('LOC24','LOC28')) -> lg 

lg %>% filter(s2chr=='LOC28')
lg %>%  group_by(s1chr) %>%   summarise(End=max(s1gp)) %>% arrange(-End) %>% rename(Chr=s1chr) %>% 
  mutate(Start=1,species='Leucoraja',fill='Leucoraja',size=12,color=252525) %>% relocate(Start,.after=Chr)-> k_leri

lg  %>% group_by(s2chr) %>%   summarise(End=max(s2gp)) %>% arrange(-End) %>% rename(Chr=s2chr) %>% 
  mutate(Start=1,species='Lepisosteus',fill='Lepisosteus',size=12,color=252525) %>% relocate(Start,.after=Chr)-> k_locu

k_leri_ord=c('LER40','LER1','LER3','LER38','LER4','LER2','LER35','LER5','LER8','LER9','LER12','LER14','LER13','LER6',
             'LER24','LER33','LER36','LER34','LER15','LER16','LER19','LER22','LER10','LER23','LER26','LER39','LER11','LER7','LER37','LER20',
             'LER17','LER25','LER28','LER21','LER18','LER31','LER30','LER29','LER27','LER32')
k_locu_ord=c('LOC4','LOC2','LOC11','LOC9','LOC1','LOC16','LOC7','LOC14','LOC17','LOC3','LOC5','LOC8','LOC10','LOC6','LOC12',
            'LOC13','LOC23','LOC20','LOC22','LOC18','LOC27','LOC21','LOC25','LOC19','LOC15','LOC26')

k_leri=k_leri[order(match(k_leri$Chr,k_leri_ord)),]
k_locu=k_locu[order(match(k_locu$Chr,k_locu_ord)),]

rbind(k_leri,k_locu) -> chroms

lg$Species_1=match(lg$s1chr,k_leri$Chr)
lg$Species_2=match(lg$s2chr,k_locu$Chr)

lg %>% filter(scol=='S') %>% select(Species_1,s1gp,Species_2,s2gp,clg) %>% rename(Start_1=s1gp,Start_2=s2gp) %>%
  mutate(End_1=Start_1+1000,End_2=Start_2+1000,fill=colors[match(clg,rownames(colors),),]) %>%
  relocate(End_1,.after=Start_1) %>% relocate(End_2,.after=Start_2) %>% select(-clg) -> lgt

ideogram(karyotype = chroms, synteny = lgt)




##### Leri vs Ggal #####

lg=read.table("Leri-Ggal_syntr.txt",h=T)
head(lg)

lg %>% mutate(s2chr=paste('GGA',s2chr,sep="")) %>%
  mutate(s1chr=gsub('Leri_','LER',gsub('C','',s1chr))) %>%
  mutate(clg=gsub(pattern = "[1-3]", replacement = "", clg)) -> lg

lg %>% filter(s1chr %notin% c('LER50','LER45','LER44')) -> lg 



lg %>%  group_by(s1chr) %>%   summarise(End=max(s1gp)) %>% arrange(-End) %>% rename(Chr=s1chr) %>% 
  mutate(Start=1,species='Leucoraja',fill='Leucoraja',size=12,color=252525) %>% relocate(Start,.after=Chr)-> k_ler

lg  %>% group_by(s2chr) %>%   summarise(End=max(s2gp)) %>% arrange(-End) %>% rename(Chr=s2chr) %>% 
  mutate(Start=1,species='Gallus',fill='Gallus',size=12,color=252525) %>% relocate(Start,.after=Chr)-> k_gal

k_gal_ord=c("GGA4","GGAZ","GGA2","GGA3","GGA5","GGA9","GGA31","GGA1","GGA10","GGA6","GGA12","GGA8","GGA30","GGA18","GGA23","GGA13","GGA7","GGA14","GGA11","GGA15","GGA19","GGA20","GGA17","GGA21","GGA28","GGAW","GGA26","GGA25","GGA22","GGA32","GGA27","GGA24","GGA33")
k_ler_ord=c("LER12","LER1","LER3","LER2","LER4","LER8","LER5","LER18","LER9","LER14","LER37","LER19","LER22","LER13","LER6","LER33","LER36","LER34","LER15","LER16","LER10","LER39","LER23","LER26","LER11","LER7","LER20","LER17","LER25","LER28","LER21","LER31","LER30","LER29","LER24","LER35","LER27","LER38","LER32","LER40")

k_ler=k_ler[order(match(k_ler$Chr,k_ler_ord)),]
k_gal=k_gal[order(match(k_gal$Chr,k_gal_ord)),]

rbind(k_ler,k_gal) -> chroms

lg$Species_1=match(lg$s1chr,k_ler$Chr)
lg$Species_2=match(lg$s2chr,k_gal$Chr)

lg %>% filter(scol=='S') %>% select(Species_1,s1gp,Species_2,s2gp,clg) %>% rename(Start_1=s1gp,Start_2=s2gp) %>%
  mutate(End_1=Start_1+1000,End_2=Start_2+1000,fill=colors[match(clg,rownames(colors),),]) %>%
  relocate(End_1,.after=Start_1) %>% relocate(End_2,.after=Start_2) %>% select(-clg) -> lgt

to_invert=c('LER1','LER10','LER2','LER9','LER22')
lg %>% mutate(inv=ifelse(s1chr %in% to_invert,'Y','N')) -> lg
lg %>% filter(scol=='S') %>% select(Species_1,s1gp,Species_2,s2gp,clg,s1chr,inv) %>% rename(Start_1=s1gp,Start_2=s2gp) %>%
  group_by(s1chr) %>% mutate(Start_1=ifelse(inv=='Y',max(Start_1)-Start_1,Start_1)) %>% ungroup() %>% select(-s1chr,-inv) %>%
  mutate(End_1=Start_1+1000,End_2=Start_2+1000,fill=colors[match(clg,rownames(colors),),]) %>%
  relocate(End_1,.after=Start_1) %>% relocate(End_2,.after=Start_2) %>% select(-clg) %>%
  arrange(Species_1,Start_1,Species_2,Start_2) -> lgt

ideogram(karyotype = chroms, synteny = as.data.frame(lgt))




##### Locu vs Ggal #####

lg=read.table("Locu-Ggal_syntra.txt",h=T)
head(lg)
#lg %>% filter(s1chr %notin% c('LER50','LER45','LER44')) %>%
#  filter(s2chr %notin% c('LOC24','LOC28')) -> lg 

lg %>% mutate(s2chr=paste('GGA',s2chr,sep="")) %>%
       mutate(s1chr=gsub('LG','LOC',s1chr)) %>%
       mutate(clg=gsub(pattern = "[1-3]", replacement = "", clg)) %>%
       filter(s2chr %notin% c('LOC28')) -> lg


lg %>%  group_by(s1chr) %>%   summarise(End=max(s1gp)) %>% arrange(-End) %>% rename(Chr=s1chr) %>% 
  mutate(Start=1,species='Lepisosteus',fill='Lepisosteus',size=12,color=252525) %>% relocate(Start,.after=Chr)-> k_loc

lg  %>% group_by(s2chr) %>%   summarise(End=max(s2gp)) %>% arrange(-End) %>% rename(Chr=s2chr) %>% 
  mutate(Start=1,species='Gallus',fill='Gallus',size=12,color=252525) %>% relocate(Start,.after=Chr)-> k_gal

k_gal_ord=c("GGA5","GGA33","GGA4","GGAZ","GGA32","GGA2","GGA22","GGA3","GGA9","GGA1","GGA26","GGA10","GGA6","GGA12","GGA8","GGA18","GGA23","GGA13","GGA7","GGA14","GGA11","GGA15","GGA19","GGA20","GGA17","GGA21","GGA28","GGA27","GGA24","GGA25","GGA16","GGAW")
k_loc_ord=c('LOC27','LOC7','LOC4','LOC2','LOC11','LOC9','LOC16','LOC1','LOC14','LOC8','LOC17','LOC3','LOC5','LOC10','LOC6','LOC12','LOC13','LOC23','LOC20','LOC22','LOC18','LOC21','LOC25','LOC19','LOC15','LOC26','LOC24','LOC28')



k_gal=k_gal[order(match(k_gal$Chr,k_gal_ord)),]
k_loc=k_loc[order(match(k_loc$Chr,k_loc_ord)),]

rbind(k_loc,k_gal) -> chroms

lg$Species_1=match(lg$s1chr,k_loc$Chr)
lg$Species_2=match(lg$s2chr,k_gal$Chr)

lg %>% filter(scol=='S') %>% select(Species_1,s1gp,Species_2,s2gp,clg) %>% rename(Start_1=s1gp,Start_2=s2gp) %>%
  mutate(End_1=Start_1+1000,End_2=Start_2+1000,fill=colors[match(clg,rownames(colors),),]) %>%
  relocate(End_1,.after=Start_1) %>% relocate(End_2,.after=Start_2) %>% select(-clg) -> lgt

ideogram(karyotype = chroms, synteny = lgt)

lg %>% select(clg)

