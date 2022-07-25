setwd("/home/aa/lyrataWGD/popGenParam/BPM/")

#BPM
pop<-c('MAY','KAR','LOM','NKM','SPN','SUN','STU','VLH','MAL','OSL','STR','PLE','STD','JOH','MOD','FOA','WAC','PEK','TEM','BZD','SCT','GYE')
con<-c('STDSTU','STDSTR','SPNWAC','SPNVLH','SPNTEM','SPNSUN','SPNSTU','SPNSTR','SPNSTD','SCTWAC','SCTVLH','SCTTEM','SCTSUN','SCTSTU','SCTSTR','SCTSTD','SCTSPN','PLEWAC','PLEVLH','PLETEM','PLESUN','PLESTU','PLESTR','PLESTD','PLESPN','PLESCT','PEKWAC','PEKVLH','PEKTEM','PEKSUN','PEKSTU','PEKSTR','PEKSTD','PEKSPN','PEKSCT','PEKPLE','OSLWAC','OSLVLH','OSLTEM','OSLSUN','OSLSTU','OSLSTR','OSLSTD','OSLSPN','OSLSCT','OSLPLE','OSLPEK','NKMWAC','NKMVLH','NKMTEM','NKMSUN','NKMSTU','NKMSTR','NKMSTD','NKMSPN','NKMSCT','NKMPLE','NKMPEK','NKMOSL','MODWAC','MODVLH','MODTEM','MODSUN','MODSTU','MODSTR','MODSTD','MODSPN','MODSCT','MODPLE','MODPEK','MODOSL','MODNKM','MAYWAC','MAYVLH','MAYTEM','MAYSUN','MAYSTU','MAYSTR','MAYSTD','MAYSPN','MAYSCT','MAYPLE','MAYPEK','MAYOSL','MAYNKM','MAYMOD','MALWAC','MALVLH','MALTEM','MALSUN','MALSTU','MALSTR','MALSTD','MALSPN','MALSCT','MALPLE','MALPEK','MALOSL','MALNKM','MALMOD','MALMAY','LOMWAC','LOMVLH','LOMTEM','LOMSUN','LOMSTU','LOMSTR','LOMSTD','LOMSPN','LOMSCT','LOMPLE','LOMPEK','LOMOSL','LOMNKM','LOMMOD','LOMMAY','LOMMAL','KARWAC','KARVLH','KARTEM','KARSUN','KARSTU','KARSTR','KARSTD','KARSPN','KARSCT','KARPLE','KARPEK','KAROSL','KARNKM','KARMOD','KARMAY','KARMAL','KARLOM','JOHWAC','JOHVLH','JOHTEM','JOHSUN','JOHSTU','JOHSTR','JOHSTD','JOHSPN','JOHSCT','JOHPLE','JOHPEK','JOHOSL','JOHNKM','JOHMOD','JOHMAY','JOHMAL','JOHLOM','JOHKAR','GYEWAC','GYEVLH','GYETEM','GYESUN','GYESTU','GYESTR','GYESTD','GYESPN','GYESCT','GYEPLE','GYEPEK','GYEOSL','GYENKM','GYEMOD','GYEMAY','GYEMAL','GYELOM','GYEKAR','GYEJOH','FOAWAC','FOAVLH','FOATEM','FOASUN','FOASTU','FOASTR','FOASTD','FOASPN','FOASCT','FOAPLE','FOAPEK','FOAOSL','FOANKM','FOAMOD','FOAMAY','FOAMAL','FOALOM','FOAKAR','FOAJOH','FOAGYE','BZDWAC','BZDVLH','BZDTEM','BZDSUN','BZDSTU','BZDSTR','BZDSTD','BZDSPN','BZDSCT','BZDPLE','BZDPEK','BZDOSL','BZDNKM','BZDMOD','BZDMAY','BZDMAL','BZDLOM','BZDKAR','BZDJOH','BZDGYE','BZDFOA')
msum<-matrix(nrow = sum(table(pop)), ncol = sum(table(pop)),dimnames =list(c(pop),c(pop)))
stats<-c("FstH","Rho") #"FstWC","dxy",
for(s in stats){
  for(contr in  con){ # contr="HRAKAM"
    p1<-substr(contr,1,3)
    p2<-substr(contr,4,7)
    p1i<-which(pop %in% paste(p1))
    p2i<-which(pop %in% paste(p2))
    d<-read.table(paste(contr,"_WS50000_MS1_BPM.txt",sep=""),h=T)
    ds<-which(colnames(d) %in% paste(s))
    msum[p1i,p2i]<-d[nrow(d),ds]
    msum[p2i,p1i]<-d[nrow(d),ds]
  }
  write.table(msum,paste("genome_",s,"_WS50000_MS1_BPM.txt",sep=""),quote = F)}




### TREEMIX
setwd("/home/aa/alpine/treemix/lyrataWGDFinal/")
library(data.table)
pop<-c("GUN","BAB","BRE","CHO","HOC","JOH","MOD","OSL","PEK","PER","SCT","STD","SUB","SWA","TEM","VLA","VLH","ZEP")
for (p in pop) { # p="WAC"
  aa<-fread(paste(p,"_tm.table",sep=""))
  aaa<-subset(aa,!aa$V3 %like% ",")
  write.table(aaa,paste("new/",p,"_tm.table",sep=""),quote = F,col.names = F,row.names = F)
}











################# Fst scan ######################
setwd("/home/aa/lyrataWGD/fstScan/")
library(data.table)
library(stringr)
library(dplyr)
library("SuperExactTest")

treshold = 0.99
tresholdGene = 0.25
genes<-fread("/home/aa/Desktop/references/lyrataV2/genesIDsLyV2.gff",h=T)
genes$ID<-substr(genes$ID,4,12)
pdf(paste("WS1/distr_SNPs",treshold,"_Genes",tresholdGene,".pdf",sep=''),width = 9,height = 10)

pops<-c("OSLTEM","STDSCT","VLHMOD", "DIPTET")


for (p in pops) { #  p = "PANWCA"
print(p)
all<-fread(paste("WS1/",p,"_WS1_MS1_BPM.txt",sep = ""),h=T)

### 1. identify 1% outlier SNPs
#all<-all[ order(all[,13],decreasing = T), ]
outl<-subset(all,all$FstH >= quantile(all$FstH,probs = treshold,na.rm = T))
write.table(x=outl,file = paste('WS1/outSNPs_',p,"_",treshold,'.txt',sep=''),append = F,quote = F,col.names = T,row.names = F)

### 2. annotate to genes 
setkey(genes, scaff, start, end)
annot<-foverlaps(x = outl, y = genes, type="within")

### 3. distr of N outlier SNPs/gene
nOutl<-as.data.frame(table(annot$ID))
#  hist(nOutl$Freq,breaks = 100)
#  summary(nOutl$Freq)
og<-subset(genes,genes$ID %in% nOutl$Var1)
og$nsnps<-nOutl$Freq
og$length<-og$end-og$start
#  hist(og$length,breaks = 100)
#  summary(og$length)
og$density<-as.numeric(og$nsnps)/as.numeric(og$length)
# hist(og$density,breaks = 100)
# summary(og$density)
og<-data.frame(og)

### 4. outliers from there - maybe even 10%?
og<-og[ order(og[,12],decreasing = T), ]
outlG<-og[1:(nrow(og)*tresholdGene),]
outlG<-subset(outlG,nsnps>2) 
print(nrow(outlG))
print(summary(outlG[,10:12]))
par(mfrow=c(2,1))
hist(og$length,breaks = 100,xlim = c(0,max(og$length)),main = p, xlab = "lenght of gene with ANY proportion of outlier SNPs")
hist(outlG$length,breaks = 100,xlim = c(0,max(og$length)), xlab = "lenght of gene with OUTLIER prop. of outlier SNPs")
par(mfrow=c(1,1))
hist(outlG$density,breaks = 100,main = p, xlab = "OUTLIER proportion of outlier SNPs")
write.table(x=outlG,file = paste('WS1/outGenes_',p,"_SNP",treshold,"_Gene",tresholdGene,'.txt',sep=''),append = F,quote = F,col.names = T,row.names = F)
}
dev.off()

### 5. overlap to see what's sensible
file.remove(paste('outGenes_SNP',treshold,"_Gene",tresholdGene,'.txt',sep=''))
for (p in pops) { #  p = "ZEPSUB"
  o<-read.table(paste('WS1/outGenes_',p,"_SNP",treshold,"_Gene",tresholdGene,'.txt',sep=''),h=T)
  write.table(t(c(p,as.character(o$ID))),paste('outGenes_SNP',treshold,"_Gene",tresholdGene,'.txt',sep=''),quote = F,col.names = F,row.names = F,sep = " ",append = T)
}


### 6. SuperExactTest
d <- scan(paste('outGenes_SNP',treshold,"_Gene",tresholdGene,'.txt',sep=''), what="", sep="\n")
d <- strsplit(d, "[[:space:]]+")
names(d) <- c("OSLTEM","STDSCT","VLHMOD", "DIPTET")
d <- lapply(d, `[`, -1)
total=34051
res=supertest(d, n=total)
res$overlap.sizes
pdf(paste('outGenes_SNP',treshold,"_Gene",tresholdGene,'.pdf',sep=''),width = 11,height = 11,pointsize = 24)
#png(paste(dat,".png",sep=""),width = 1200,height = 960,pointsize = 24)
plot(res, Layout="landscape", sort.by="size", margin=c(0.5,5,1,2),keep.empty.intersections=F,min.intersection.size = 1, degree=1:4,show.overlap.size = T,color.on="black")
dev.off()
write.csv(summary(res)$Table, file=paste('outGenes_SNP',treshold,"_Gene",tresholdGene,'.csv',sep=''), row.names=FALSE)

  
########### FUNCTIONAL FOLLOW-UP 
#INDEPENDENT FROM THE ABOVE - ASSUMES ONE HAVE A GOOD CANDIDATE LIST

setwd("/home/aa/lyrataWGD/fstScan/")
library(data.table)

### 1. Candidate gene annotation TAIR

dict<-fread("/home/aa/Desktop/references/lyrataV2/functions/ALATdict.txt")
s<-read.csv("outGenes_SNP0.99_Gene0.1.csv",h=T)

# s<-read.table("WS1/outGenes_PANWCA_SNP0.99_Gene0.25.txt",h=T)

ann<-fread("/home/aa/Desktop/references/lyrataV2/LyV2_TAIR10orth_des_20150927.txt",h=T,quote="")
ann2018<-fread("/home/aa/Desktop/references/lyrataV2/LyV2_TAIR11orth_des_20181231.txt",h=T,quote="")
s11<-droplevels(subset(s,s$Degree == 1))
s<-subset(s,s$Degree == 1) ####Can change to 2, 3,...
s1<-toString(s$Elements)
ss<-unlist(strsplit(s1, ", "))
parCand<-subset(dict,dict$AL %in% ss) 
#  parCand<-subset(dict,dict$AL %in% s$ID) 

parCand<-parCand[!duplicated(parCand[,c('AL')]),] 
sel<-substr(parCand$AT,1,9)
sumPar<-matrix(nrow = nrow(parCand),ncol = 8,dimnames = list(c(),c("AL","AT","lineages","T15_name","T15_descr","T18_descrShort","T18_descrCurat","T18_descrComp")))
for (i in 1:nrow(parCand)) { # i=1
  g<-unlist(parCand[i,1])
  ann1<-subset(ann, ann$`Version-2` %in% g)
  ann2018_1<-subset(ann2018, ann2018$AL %in% g)
pair<-droplevels(subset(s11,s11$Elements %like% g))
pairc<-paste(pair$Intersections, sep="", collapse=', ')
sumPar[i,1]<-unlist(parCand[i,1])
sumPar[i,2]<-unlist(parCand[i,2])
sumPar[i,3]<-pairc
sumPar[i,4]<-unlist(ann1[1,7])
sumPar[i,5]<-unlist(ann1[1,4])
sumPar[i,6]<-unlist(ann2018_1[1,4])
sumPar[i,7]<-unlist(ann2018_1[1,5])
sumPar[i,8]<-unlist(ann2018_1[1,6])
}

#cyclins<-subset(sumPar, sumPar[,2] %in% cand$AT)
#write.table(x=cyclins,file = paste('cyclins_SNP0.99_Gene0.1.ALAT.ann.any.txt',sep=''),append = F,quote = F,col.names = T,row.names = F,sep='\t')


write.table(x=sumPar,file = paste('outGenes_SNP0.99_Gene0.25.ALAT.ann.tripleAAAL.txt',sep=''),append = F,quote = F,col.names = T,row.names = F,sep='\t')
write.table(x=sumPar[,1:3],file = paste('outGenes_SNP0.99_Gene0.1.ALAT.tripleAAAL.txt',sep=''),append = F,quote = F,col.names = T,row.names = F,sep='\t')



### 2. topGO
library("biomaRt")
library(topGO)

mart <- biomaRt::useMart(biomart = "plants_mart",dataset = "athaliana_eg_gene", host = 'plants.ensembl.org')
GTOGO <- biomaRt::getBM(attributes = c( "ensembl_gene_id", "go_id"), mart = mart)
GTOGO <- GTOGO[GTOGO$go_id != '',]
geneID2GO <- by(GTOGO$go_id,GTOGO$ensembl_gene_id,function(x) as.character(x))
all.genes <- sort(unique(as.character(GTOGO$ensembl_gene_id)))
int.genes <- factor(as.integer(all.genes %in% sel))
names(int.genes) = all.genes
go.obj <- new("topGOdata", ontology='BP', allGenes = int.genes, annot = annFUN.gene2GO, gene2GO = geneID2GO,nodeSize=10) ## 
resultsFe <- runTest(go.obj, algorithm = "elim", statistic = "fisher") #More conservative
resultsFc <- runTest(go.obj, algorithm = "classic", statistic = "fisher")
allRes <- GenTable(go.obj, classicFisher = resultsFc, elimFisher = resultsFe, orderBy = "elimFisher", ranksOf = "elimFisher", topNodes = 500)

allRes$genes <- sapply(allRes$GO.ID, function(x)
{genes<-genesInTerm(go.obj, x) 
genes[[1]][genes[[1]] %in% sel]})

a<-as.data.frame(subset(allRes, allRes$Annotated>=10 & allRes$Annotated <= 500 & allRes$elimFisher <= 0.05 & allRes$Significant > 1))

a$Genes<-""
for (i in 1:length(a$genes)) {
  a$Genes[i]<-paste(unlist(a$genes[i]), sep="", collapse=', ')}
a<-a[,-8]

write.table(a,file=paste("topgo_outGenes_SNP0.99_Gene0.1_min10max500.txt",sep=""),sep="\t",row.names=F)


### 3. visualize - dotplots
setwd("/home/aa/lyrataWGD/fstScan/")
library(data.table)
library(dplyr)
library(gplots)
library(RColorBrewer)
library(scales)
#1. nacist vsechny pary, gff
STDSCT<-fread("WS1/STDSCT_WS1_MS1_BPM.txt")
OSLTEM<-fread("WS1/OSLTEM_WS1_MS1_BPM.txt")
VLHMOD<-fread("WS1/VLHMOD_WS1_MS1_BPM.txt")
DIPTET<-fread("WS1/DIPTET_WS1_MS1_BPM.txt")

genes<-fread("/home/aa/Desktop/references/lyrataV2/genesIDsLyV2.gff",h=T)
genes$ID<-as.character(substr(genes$ID,4,12))
cand<-read.table("outGenes_SNP0.99_Gene0.1.ALAT.tripleAAAL.txt",sep="\t",h=T)
cand<-read.table("17tetraploidadaptive.txt",sep="\t",h=T)

flank<-15000

pdf(paste("Figures/all_",flank,".17tetraploidadaptive.pdf",sep=''),width = 6,height = 13,pointsize = 17)

for (i in 1:nrow(cand)) { # i=1
gen<-cand[i,1]  #as.character(droplevels(
int<-subset(genes,genes$ID %in% gen)
first<-int$start-flank
last<-int$end+flank
STDSCTgene<-subset(STDSCT, STDSCT$end > first & STDSCT$end<last)
OSLTEMgene<-subset(OSLTEM, OSLTEM$end > first & OSLTEM$end<last)
VLHMODgene<-subset(VLHMOD, VLHMOD$end > first & VLHMOD$end<last)
DIPTETgene<-subset(DIPTET, DIPTET$end > first & DIPTET$end<last)

minpos<-min(c(STDSCTgene$end,OSLTEMgene$end,VLHMODgene$end,DIPTETgene$end))
maxpos<-max(c(STDSCTgene$end,OSLTEMgene$end,VLHMODgene$end,DIPTETgene$end))
maxFst<-max(c(STDSCTgene$FstH,OSLTEMgene$FstH,VLHMODgene$FstH,DIPTETgene$FstH))


#pdf(paste("Figures/",gen,"_",flank,".pdf",sep=''),width = 6,height = 20,pointsize = 17)
par(mfrow=c(4,1))
par(mar=c(0.1,2.7,0.1,1), mgp=c(1, 0.5, 0))
par(bg=NA) 
plot(DIPTETgene$FstH~DIPTETgene$end,ylab = "",col = alpha("grey25", 0.8),ylim=c(0.01,maxFst+0.05),xlab="",pch=19,cex=1,xlim=c(minpos,maxpos),xaxt='n', ann=FALSE)
title(ylab="Fst DIPTET arenosa", line=1.4, cex.lab=1.2)
if (int$orientation %in% "+") {
  arrows(x0 = int$start,y0 = maxFst+0.04,x1 = int$end,y1 =maxFst+0.04,col = "red",lwd = 4,code = 2)
} else {arrows(x0 = int$end,y0 = maxFst+0.04,x1 = int$start,y1 =maxFst+0.04,col = "red",lwd =4,code=2)} 


plot(VLHMODgene$FstH~VLHMODgene$end,ylab = "",col = alpha("grey25", 0.8),ylim=c(0.01,maxFst+0.05),xlab="",pch=19,cex=1,xlim=c(minpos,maxpos),xaxt='n', ann=FALSE)
title(ylab="Fst VLHMOD", line=1.4, cex.lab=1.2)
segments(x0 = int$start,y0 = maxFst+0.04,x1 = int$end,y1 =maxFst+0.04,col = "grey50",lwd = 4)

plot(OSLTEMgene$FstH~OSLTEMgene$end,ylab = "",col = alpha("grey25", 0.8),ylim=c(0.01,maxFst+0.05),xlab="",pch=19,cex=1,xlim=c(minpos,maxpos),xaxt='n', ann=FALSE)
title(ylab="Fst OSLTEM", line=1.4, cex.lab=1.2)
segments(x0 = int$start,y0 = maxFst+0.04,x1 = int$end,y1 =maxFst+0.04,col = "grey50",lwd = 4)

par(mar=c(2.5,2.7,0.1,1))
plot(STDSCTgene$FstH~STDSCTgene$end,ylab = "",col = alpha("grey25", 0.8),ylim=c(0.01,maxFst+0.05),xlab="",pch=19,cex=1,xlim=c(minpos,maxpos))
title(ylab="Fst STDSCT", line=1.4, cex.lab=1.2)
title(xlab=paste(gen,cand[i,2],cand[i,3],sep= " - "), line=1.4, cex.lab=1.2)
segments(x0 = int$start,y0 = maxFst+0.04,x1 = int$end,y1 =maxFst+0.04,col = "grey50",lwd = 4)
#dev.off()
}

dev.off()


#Visualize - heatmap - .table.recode.txt 








### 3. visualize sanity check scaffs
VLHMODs1<-subset(VLHMOD,scaff%in% "scaffold_1")
VLHMODs2<-subset(VLHMOD,scaff%in% "scaffold_2")
VLHMODs3<-subset(VLHMOD,scaff%in% "scaffold_3")
VLHMODs4<-subset(VLHMOD,scaff%in% "scaffold_4")
VLHMODs5<-subset(VLHMOD,scaff%in% "scaffold_5")
VLHMODs6<-subset(VLHMOD,scaff%in% "scaffold_6")
VLHMODs7<-subset(VLHMOD,scaff%in% "scaffold_7")
VLHMODs8<-subset(VLHMOD,scaff%in% "scaffold_8")
png('VLHMOD.png',width = 2400,height = 3000)
par(mfrow=c(8,1))
plot(VLHMODs1$FstH~VLHMODs1$end)
plot(VLHMODs2$FstH~VLHMODs2$end)
plot(VLHMODs3$FstH~VLHMODs3$end)
plot(VLHMODs4$FstH~VLHMODs4$end)
plot(VLHMODs5$FstH~VLHMODs5$end)
plot(VLHMODs6$FstH~VLHMODs6$end)
plot(VLHMODs7$FstH~VLHMODs7$end)
plot(VLHMODs8$FstH~VLHMODs8$end)
dev.off()
  
STDSCTs1<-subset(STDSCT,scaff%in% "scaffold_1")
STDSCTs2<-subset(STDSCT,scaff%in% "scaffold_2")
STDSCTs3<-subset(STDSCT,scaff%in% "scaffold_3")
STDSCTs4<-subset(STDSCT,scaff%in% "scaffold_4")
STDSCTs5<-subset(STDSCT,scaff%in% "scaffold_5")
STDSCTs6<-subset(STDSCT,scaff%in% "scaffold_6")
STDSCTs7<-subset(STDSCT,scaff%in% "scaffold_7")
STDSCTs8<-subset(STDSCT,scaff%in% "scaffold_8")
png('STDSCT.png',width = 2400,height = 3000)
par(mfrow=c(8,1))
plot(STDSCTs1$FstH~STDSCTs1$end)
plot(STDSCTs2$FstH~STDSCTs2$end)
plot(STDSCTs3$FstH~STDSCTs3$end)
plot(STDSCTs4$FstH~STDSCTs4$end)
plot(STDSCTs5$FstH~STDSCTs5$end)
plot(STDSCTs6$FstH~STDSCTs6$end)
plot(STDSCTs7$FstH~STDSCTs7$end)
plot(STDSCTs8$FstH~STDSCTs8$end)
dev.off()

OSLTEMs1<-subset(OSLTEM,scaff%in% "scaffold_1")
OSLTEMs2<-subset(OSLTEM,scaff%in% "scaffold_2")
OSLTEMs3<-subset(OSLTEM,scaff%in% "scaffold_3")
OSLTEMs4<-subset(OSLTEM,scaff%in% "scaffold_4")
OSLTEMs5<-subset(OSLTEM,scaff%in% "scaffold_5")
OSLTEMs6<-subset(OSLTEM,scaff%in% "scaffold_6")
OSLTEMs7<-subset(OSLTEM,scaff%in% "scaffold_7")
OSLTEMs8<-subset(OSLTEM,scaff%in% "scaffold_8")
png('OSLTEM.png',width = 2400,height = 3000)
par(mfrow=c(8,1))
plot(OSLTEMs1$FstH~OSLTEMs1$end)
plot(OSLTEMs2$FstH~OSLTEMs2$end)
plot(OSLTEMs3$FstH~OSLTEMs3$end)
plot(OSLTEMs4$FstH~OSLTEMs4$end)
plot(OSLTEMs5$FstH~OSLTEMs5$end)
plot(OSLTEMs6$FstH~OSLTEMs6$end)
plot(OSLTEMs7$FstH~OSLTEMs7$end)
plot(OSLTEMs8$FstH~OSLTEMs8$end)
dev.off()


#### CYCA pathway
aaa<-read.table("file:///home/aa/lyrataWGD/fstScan/string_CYCA23_physical.tsv",h=T)
aaa$x<-substr(aaa$node1_string_id,6,14)
aaa$y<-substr(aaa$node2_string_id,6,14)
ss<-rbind(aaa$x,aaa$y)
aa<-ss[order(ss,decreasing = F)]
ss<-data.frame(aa[!duplicated(aa)])
###visualize - dotplots
setwd("/home/aa/lyrataWGD/fstScan/")
library(data.table)
library(dplyr)
library(gplots)
library(RColorBrewer)
library(scales)
#1. nacist vsechny pary, gff
STDSCT<-fread("WS1/STDSCT_WS1_MS1_BPM.txt")
OSLTEM<-fread("WS1/OSLTEM_WS1_MS1_BPM.txt")
VLHMOD<-fread("WS1/VLHMOD_WS1_MS1_BPM.txt")
DIPTET<-fread("WS1/DIPTET_WS1_MS1_BPM.txt")

genes<-fread("/home/aa/Desktop/references/lyrataV2/genesIDsLyV2.gff",h=T)
genes$ID<-as.character(substr(genes$ID,4,12))
dict<-fread("/home/aa/Desktop/references/lyrataV2/functions/ALATdict.txt")
cand<-subset(dict, substr(dict$AT,1,9) %in% ss$aa..duplicated.aa..)
flank<-15000

pdf(paste("Figures/all_",flank,".CYCAphysicalInt.pdf",sep=''),width = 6,height = 13,pointsize = 17)

for (i in 1:nrow(cand)) { # i=1
  gen<-cand[i,1]  #as.character(droplevels(
  int<-subset(genes,genes$ID %in% gen)
  first<-int$start-flank
  last<-int$end+flank
  STDSCTgene<-subset(STDSCT, STDSCT$end > first & STDSCT$end<last)
  OSLTEMgene<-subset(OSLTEM, OSLTEM$end > first & OSLTEM$end<last)
  VLHMODgene<-subset(VLHMOD, VLHMOD$end > first & VLHMOD$end<last)
  DIPTETgene<-subset(DIPTET, DIPTET$end > first & DIPTET$end<last)
  
  minpos<-min(c(STDSCTgene$end,OSLTEMgene$end,VLHMODgene$end,DIPTETgene$end))
  maxpos<-max(c(STDSCTgene$end,OSLTEMgene$end,VLHMODgene$end,DIPTETgene$end))
  maxFst<-max(c(STDSCTgene$FstH,OSLTEMgene$FstH,VLHMODgene$FstH,DIPTETgene$FstH))
  
  
  #pdf(paste("Figures/",gen,"_",flank,".pdf",sep=''),width = 6,height = 20,pointsize = 17)
  par(mfrow=c(4,1))
  par(mar=c(0.1,2.7,0.1,1), mgp=c(1, 0.5, 0))
  par(bg=NA) 
  plot(DIPTETgene$FstH~DIPTETgene$end,ylab = "",col = alpha("grey25", 0.8),ylim=c(0.01,maxFst+0.05),xlab="",pch=19,cex=1,xlim=c(minpos,maxpos),xaxt='n', ann=FALSE)
  title(ylab="Fst DIPTET arenosa", line=1.4, cex.lab=1.2)
  if (int$orientation %in% "+") {
    arrows(x0 = int$start,y0 = maxFst+0.04,x1 = int$end,y1 =maxFst+0.04,col = "red",lwd = 4,code = 2)
  } else {arrows(x0 = int$end,y0 = maxFst+0.04,x1 = int$start,y1 =maxFst+0.04,col = "red",lwd =4,code=2)} 
  
  
  plot(VLHMODgene$FstH~VLHMODgene$end,ylab = "",col = alpha("grey25", 0.8),ylim=c(0.01,maxFst+0.05),xlab="",pch=19,cex=1,xlim=c(minpos,maxpos),xaxt='n', ann=FALSE)
  title(ylab="Fst VLHMOD", line=1.4, cex.lab=1.2)
  segments(x0 = int$start,y0 = maxFst+0.04,x1 = int$end,y1 =maxFst+0.04,col = "grey50",lwd = 4)
  
  plot(OSLTEMgene$FstH~OSLTEMgene$end,ylab = "",col = alpha("grey25", 0.8),ylim=c(0.01,maxFst+0.05),xlab="",pch=19,cex=1,xlim=c(minpos,maxpos),xaxt='n', ann=FALSE)
  title(ylab="Fst OSLTEM", line=1.4, cex.lab=1.2)
  segments(x0 = int$start,y0 = maxFst+0.04,x1 = int$end,y1 =maxFst+0.04,col = "grey50",lwd = 4)
  
  par(mar=c(2.5,2.7,0.1,1))
  plot(STDSCTgene$FstH~STDSCTgene$end,ylab = "",col = alpha("grey25", 0.8),ylim=c(0.01,maxFst+0.05),xlab="",pch=19,cex=1,xlim=c(minpos,maxpos))
  title(ylab="Fst STDSCT", line=1.4, cex.lab=1.2)
  title(xlab=paste(gen,cand[i,2],sep= " - "), line=1.4, cex.lab=1.2)
  segments(x0 = int$start,y0 = maxFst+0.04,x1 = int$end,y1 =maxFst+0.04,col = "grey50",lwd = 4)
  #dev.off()
}

dev.off()





### ABBA-BABA ###
setwd("/home/aa/lyrataWGD/ABBABABA/")
abba = function(p1, p2, p3) (1 - p1) * p2 * p3
baba = function(p1, p2, p3) p1 * (1 - p2) * p3
D.stat = function(dataframe) (sum(dataframe$ABBA) - sum(dataframe$BABA)) / (sum(dataframe$ABBA) + sum(dataframe$BABA))

freq_table = read.table("data/fourfold.dp8nc.m0.5.all_var_bial_singl_outGUN.geno.tsv.gz", header=T, as.is=T)
#freq_table = read.table("hel92.DP8MP4BIMAC2HET75dist250.derFreq.tsv.gz", header=T, as.is=T)

freq_table<- na.omit(freq_table)
nrow(freq_table)

P1 = "STD"
P2 = "SCT"
P3 = "SWA"




ABBA = abba(freq_table[,P1], freq_table[,P2], freq_table[,P3])
BABA = baba(freq_table[,P1], freq_table[,P2], freq_table[,P3])
ABBA_BABA_df = as.data.frame(cbind(ABBA,BABA))
#ABBA_BABA_df1<- na.omit(ABBA_BABA_df)
D = D.stat(ABBA_BABA_df)
sum(ABBA_BABA_df$ABBA)
sum(ABBA_BABA_df$BABA)

P1 = "OSL"
P2 = "TEM"
P3 = "TET"

P1 = "LOM"
P2 = "TEM"
P3 = "TET"
ABBA = abba(freq_table[,P1], freq_table[,P2], freq_table[,P3])
BABA = baba(freq_table[,P1], freq_table[,P2], freq_table[,P3])
ABBA_BABA_df = as.data.frame(cbind(ABBA,BABA))
ABBA_BABA_df1<- na.omit(ABBA_BABA_df)
D = D.stat(ABBA_BABA_df1)
sum(ABBA_BABA_df1$ABBA)
sum(ABBA_BABA_df1$BABA)

P1 = "VLH"
P2 = "MOD"
P3 = "TET"

P1 = "LOM"
P2 = "MOD"
P3 = "TET"


P1 = "BDO"
P2 = "TET"
P3 = "MOD"

ABBA = abba(freq_table[,P1], freq_table[,P2], freq_table[,P3])
BABA = baba(freq_table[,P1], freq_table[,P2], freq_table[,P3])
ABBA_BABA_df = as.data.frame(cbind(ABBA,BABA))
ABBA_BABA_df1<- na.omit(ABBA_BABA_df)
D = D.stat(ABBA_BABA_df1)
sum(ABBA_BABA_df1$ABBA)
sum(ABBA_BABA_df1$BABA)

#Jackknife
source("jackknife_Verca.R")
chrom_table = read.table("data/scaffold_lengths.txt")
chrom_lengths = chrom_table[,2]
names(chrom_lengths) = chrom_table[,1]


block_indices <- get_block_indices(block_size=1e6,
                                   positions=freq_table$position,
                                   chromosomes=freq_table$scaffold)

n_blocks <- length(block_indices)
print(paste("Genome divided into", n_blocks, "blocks."))


aaa<-block_indices[lapply(block_indices,length)>0]

D_sd <- get_jackknife_sd(block_indices=aaa,
                         FUN=D.stat,
                         freq_table[,P1], freq_table[,P2], freq_table[,P3])
print(paste("D standard deviation = ", round(D_sd,4)))

block_jackknife(block_indices=aaa,
                FUN=D.stat,
                freq_table[,P1], freq_table[,P2], freq_table[,P3])









### HEATMAPS ###
### ID list in i2 <- tt2$Category[2:nrow(tt2)]
### ID list in i2 <- c("AL5G34820"),"AL5G34820")
#
##17 tetraploid adaptive
i2<-c("AL1G27690","AL1G35730","AL1G56960","AL2G25920","AL2G37810","AL3G24370","AL3G29960","AL4G29630","AL4G46460","AL5G32860","AL7G13140","AL7G35790","AL8G25600","AL8G44080","AL1G26770","AL1G62040","AL6G15380")

i2<-c("AL7G15530")

library(data.table)
library(dplyr)
library(gplots)
library(RColorBrewer)
setwd("/home/aa/lyrataWGD/fstScan/heatmap/")
#1.extract genes from file
ar<-fread(paste('ann/ALL.table.recode.txt',sep=''),h=F,na.strings = "-9",nThread = 3)
nnn<-fread("PopKey_lyrenosa.heatmap.csv",h=F)
names<-c(nnn[2:nrow(nnn),1])[[1]]
colnames(ar) <- c("pop","ploidy","scaff","start","AN","DP","ID","ann","aas",names,"nan")
ar1<-subset(ar,ar$ID %in% i2)[,1:(ncol(ar)-1)]
#2. filter too low freq
i<-which(colnames(ar1) %like% "_")[1]
y<-which(colnames(ar1) %like% "_")[length(which(colnames(ar1) %like% "_"))]
afh1<-ar1[,i:as.numeric(y)]
ar1$ACh1<-rowSums(afh1,na.rm = T)
ar1$NAh1<-apply(is.na(afh1), 1, sum)
ar1$ANh1<-486 ####hard calculated!
ar1$tot<-ar1$ACh1/(ar1$ANh1-ar1$NAh1)
ar2<-subset(ar1,ar1$tot > 16/486)
ar3<-subset(ar2,ar2$tot < as.numeric(1-(16/486)))
#3. calculate AF for each lineage TODO 
popnames<-unique(substr(names,1,3))
for (pop in popnames) { #  pop="BAB"
  afh1 <- ar3 %>% dplyr:: select(starts_with(pop))
  ar3$ACh1<-rowSums(afh1,na.rm = T)
  ar3$NAh1<-apply(is.na(afh1), 1, sum)
  if (subset(nnn,nnn$V1%like%pop)[1,3]==0) 
  {ar3$ANh1<-(ncol(afh1)-ar3$NAh1)*2} else {ar3$ANh1<-(ncol(afh1)-ar3$NAh1)*4}
  ar3$POP<-ar3$ACh1/ar3$ANh1
  colnames(ar3)[which(colnames(ar3) %in% "POP")]<-pop
}
#     write.table(ar3,"CAPD3.txt",row.names = F,quote = F,sep="\t")
s4<-dplyr::select(ar3,ID,ann,aas,ZEP, BAB, SUB, HOC, PER, BRE, VLA, CHO, SWA, JOH, MOD,TEM,PEK,SCT,VLH,OSL,STD,tot)
### MORE POSSIBILITIES ###
#s<-subset(s4,!s4$ann %in% "intragenic_variant" & !s4$ann %in% "downstream_gene_variant" & !s4$ann %in% "upstream_gene_variant")
s<-subset(s4,!s4$ann %in% "intragenic_variant")
#s<-subset(s4,s4$ann %like% "missense_variant")
#repolarise
#for (i in  1:nrow(s)){ # i=1
#  if (s$tot[i]>0.5)
#  {s[i,4:(length(s[i,])-1)]<-1-s[i,4:(length(s[i,])-1)]
#  } else {}}
#repolarise
popnames<-colnames(s4)[4:(length(colnames(s4))-1)]
for (i in  1:nrow(s)){ # i=1
  d<-c()
  t<-c()
  for (pop in popnames) { #  pop="BAB"
    ppp <- s %>% dplyr:: select(starts_with(pop))
    if (subset(nnn,nnn$V1%like%pop)[1,3]==0) 
    {d<-c(d,as.numeric(as.character((ppp[i,1])[[1]])))} else {t<-c(t,as.numeric(as.character((ppp[i,1])[[1]])))}
    }
  meandip<-mean(d,na.rm = T)
  meantet<-mean(t,na.rm = T)
  if (meandip>meantet)
  {s[i,4:(length(s[i,])-1)]<-(1-s[i,4:(length(s[i,])-1)])
  print(i)
  } else {}}
#plot
ann<-fread("/home/aa/Desktop/references/lyrataV2/LyV2_TAIR11orth_des_20171231.txt")
pops<-c(colnames(s)[4:(length(colnames(s))-1)])
my_palette <- colorRampPalette(c("khaki1", "green2", "blue3"))(n = 100)
for (id in i2){ #   id = "AL1G48330"
  if (nrow(subset(s,s$ID %in% id))>1)
  {s1<-subset(s,s$ID %in% id)
  ann1<-subset(ann,ann$AL %in% id)
  if (nrow(s1)>110)
  {p=0.75
  } else {p=1}
  s5<-dplyr::select(s1,SUB, HOC, PER, JOH, MOD, VLH, BAB, BRE, VLA, TEM, PEK, OSL, ZEP, CHO, SWA, SCT, STD)
  df<-as.matrix(s1[,4:(length(s1[i,])-1)],rownames = s1$aas)
  df1<-as.matrix(s5,rownames = s1$aas)
  #parallel
  pdf(paste("heatmap_",id,".pdf",sep=""),height = 11,width = 9.3)
  heatmap.2(x = df,dendrogram = "none",Colv="NA", Rowv="NA",key = F,col=my_palette,colsep= c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,14,17,17,17,17,17,17,17,17,17),sepcolor= c("black"),sepwidth = c(0.02),trace="none",ColSideColors = c("red","red","red","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","red","red","red"),labCol=pops,colCol= c("tomato","tomato","tomato","skyblue2","skyblue2","skyblue2","skyblue2","skyblue2","skyblue2","blue","blue","blue","blue","blue","red","red","red"),offsetRow = c(0.05),offsetCol= c(0.05),cexRow = c(p),cexCol = c(2),margins = c(7,7),lhei = c(0.2,9),lwid = c(0.03,0.8),xlab = paste(id,ann1$AT,ann1$annShort,sep=" - ")) 
  #independent
  heatmap.2(x = df1,dendrogram = "none",Colv="NA", Rowv="NA",key = F,col=my_palette,colsep= c(0,0,6,12,17),sepcolor= c("black"),sepwidth = c(0.02),trace="none",ColSideColors = c("red","blue","blue","blue","blue","red","red","blue","blue","blue","blue","red","red","blue","blue","blue","red"),labCol=colnames(s5),colCol= c("tomato","skyblue2","skyblue2","blue","blue","red","tomato","skyblue2","skyblue2","blue","blue","red","tomato","blue","blue",'skyblue2',"red"),offsetRow = c(0.05),offsetCol= c(0.05),cexRow = c(p),cexCol = c(2),margins = c(7,7),lhei = c(0.2,9),lwid = c(0.03,0.8),xlab = paste(id,ann1$AT,ann1$annShort,sep=" - ")) 
  dev.off()
  } else {} 
}



#####################################
#########HaplotypeAF#################
#####################################

#### 0. get intervals - once for candidate process/function
library(data.table)
genes<-fread("/home/aa/Desktop/references/lyrataV2/genesIDsLyV2.gff",h=T)
genes$ID<-substr(genes$ID,4,12)
int<-c('AL1G27690','AL1G35730','AL1G56960','AL2G25920','AL2G37810','AL3G24370','AL3G29960','AL4G29630','AL4G46460','AL5G32860','AL7G13140','AL7G35790','AL8G25600','AL8G44080','AL1G26770','AL1G62040','AL6G15380')
l<-subset(genes, genes$ID %in% int)
write.table(paste("-L ",l$scaff,":",l$start-5000,"-",l$end+5000,sep=""),"/home/aa/lyrataWGD/fstScan/haplotypeMaps/17WGDcand_intervals10May.txt",col.names = F,row.names = F,quote = F) 
#this is then used in script fish://holcovam@tilia.metacentrum.cz/auto/pruhonice1-ibot/home/holcovam/ScanTools/lyrenosa.merged.annotated/filterLyrAreAnn17cand.sh - serves as info for the haplotype maps (ALL.ann)

### 0b. To get candidate HISs
library(data.table)
library(dplyr)
library(data.table)
library(dplyr)
library(ggmap)
library(argosfilter)
library(mapplots)
library(RgoogleMaps)
library(RCurl) 
setwd("/home/aa/lyrataWGD/fstScan/heatmap/")
namesAll<-read.table("file:///home/aa/lyrataWGD/fstScan/haplotypeMaps/17WGDcand_names.txt",h=T)
allall<-fread(paste('ann/ALL.table.recode.txt',sep=''),h=F,na.strings = "-9",nThread = 3)
nnn<-fread("PopKey_lyrenosa.heatmap.csv",h=F)
names<-c(nnn[2:nrow(nnn),1])[[1]]
colnames(allall) <- c("pop","ploidy","scaff","start","AN","DP","ID","ann","aas",names,"nan")
thresh<-0.6

for (i in 1:nrow(namesAll)) { #  i=1 name<-"SDS" i2<-"AL1G26770"
  i2<-namesAll[i,1]
  name<-namesAll[i,2]
#plotHaplotypeAFs <- function(i2,name) {
#  tryCatch({
#0.identify HIS
ar1<-subset(allall,allall$ID %in% i2)[,1:(ncol(allall)-1)]
i<-which(colnames(ar1) %like% "_")[1]
y<-which(colnames(ar1) %like% "_")[length(which(colnames(ar1) %like% "_"))]
afh1<-ar1[,i:as.numeric(y)]
ar1$ACh1<-rowSums(afh1,na.rm = T)
ar1$NAh1<-apply(is.na(afh1), 1, sum)
ar1$ANh1<-486 ####hard calculated!
ar1$tot<-ar1$ACh1/(ar1$ANh1-ar1$NAh1)
ar2<-subset(ar1,ar1$tot > 16/486)
ar3<-subset(ar2,ar2$tot < as.numeric(1-(16/486)))
popnames<-unique(substr(names,1,3))
for (pop in popnames) { #  pop="BAB"
  afh1 <- ar3 %>% dplyr:: select(starts_with(pop))
  ar3$ACh1<-rowSums(afh1,na.rm = T)
  ar3$NAh1<-apply(is.na(afh1), 1, sum)
  if (subset(nnn,nnn$V1%like%pop)[1,3]==0) 
  {ar3$ANh1<-(ncol(afh1)-ar3$NAh1)*2} else {ar3$ANh1<-(ncol(afh1)-ar3$NAh1)*4}
  ar3$POP<-ar3$ACh1/ar3$ANh1
  colnames(ar3)[which(colnames(ar3) %in% "POP")]<-pop
}
#ar4<-subset(ar3, ar3$ann %like% "missense")
#ar4<-subset(ar3, !ar3$ann %in% "synonymous_variant")
ar4<-ar3

ar5<-dplyr::select(ar4,scaff,start,ID,ann,aas,ZEP, BAB, SUB, HOC, PER, BRE, VLA, CHO, SWA, JOH, MOD,TEM,PEK,SCT,VLH,OSL,STD,tot)
dipa<-dplyr::select(ar5,ZEP, BAB, SUB)
teta<-dplyr::select(ar5,HOC, PER, BRE, VLA, CHO, SWA)
ar5$meandipa<-apply(dipa,MARGIN = 1,mean)
ar5$meanteta<-apply(teta,MARGIN = 1,mean)
ar5$afda<-abs(ar5$meandipa-ar5$meanteta)
ar6a<-subset(ar5,ar5$afda > thresh)
print(nrow(ar6a))
dipl<-dplyr::select(ar5,VLH,OSL,STD)
tetl<-dplyr::select(ar5, JOH, MOD,TEM,PEK,SCT)
ar5$meandipl<-apply(dipl,MARGIN = 1,mean)
ar5$meantetl<-apply(tetl,MARGIN = 1,mean)
ar5$afdl<-abs(ar5$meandipl-ar5$meantetl)
ar6l<-subset(ar5,ar5$afdl > thresh)
print(nrow(ar6l))
overl<-subset(ar6a, ar6a$scaff %in% ar6l$scaff & ar6a$start %in% ar6l$start)
print(nrow(overl))
write.table(ar6a,paste("../haplotypeMaps/",name,".HIS.arenosa.txt",sep = ""),row.names = F,quote = F,sep="\t")
write.table(ar6l,paste("../haplotypeMaps/",name,".HIS.lyrata.txt",sep = ""),row.names = F,quote = F,sep="\t")
write.table(overl,paste("../haplotypeMaps/",name,".HIS.both.txt",sep = ""),row.names = F,quote = F,sep="\t")
#union
both<-rbind(ar6a[,1:5],ar6l[,1:5])
both<-both[ order(both[,2]), ]
bo<-both[!duplicated(both[,c('scaff','start')]),]
write.table(bo,paste("../haplotypeMaps/",name,".HIS.union.txt",sep = ""),row.names = F,quote = F,sep="\t")
}

####1a. extract HISs from the tolat dataset - arenosa
setwd("/home/aa/lyrataWGD/fstScan/haplotypeMaps/")
ar<-fread(paste('ALL.ann/ALL.table.recode.txt',sep=''),h=F,na.strings = "-9",nThread = 3)
popkey<-fread("PopKeyAll1000.csv",h=T)
namesAll<-read.table("file:///home/aa/lyrataWGD/fstScan/haplotypeMaps/17WGDcand_names.txt",h=T)
namesAll<-namesAll[10,]

for (i in 1:nrow(namesAll)) { #  i=1
  i2<-namesAll[i,1]
  name<-namesAll[i,2]
hifs<-fread(paste(name,".HIS.arenosa.txt",sep = ""),h=T)
prot <- name
names<-c(popkey[1:nrow(popkey),1])[[1]]
colnames(ar) <- c("pop","ploidy","scaff","start","AN","DP","ID","ann","AAS",names,"nan")
ar1<-subset(ar,ar$scaff %in% hifs$scaff & ar$start %in% hifs$start)
ar2<-as.data.frame(ar1[,c(3,4)])
popnames<-unique(substr(names,1,3))
popnames<-popnames[!popnames %in% "XXX"]
ploidy <- data.frame(pops = popnames, ploidy.lev = NA)
for(i in 1:length(popnames)){  # i = 1
  afh1 <- ar1 %>% dplyr:: select(starts_with(popnames[i]))
  ar1$ACh1 <- rowSums(afh1,na.rm = T)
  ar1$NAh1 <- apply(is.na(afh1), 1, sum)
  if (subset(popkey,popkey$Samples %like% popnames[i])[1,3]==0) 
  {ar1$ANh1<-(ncol(afh1)-ar1$NAh1)*2
  ploidy[i,2] <- 2
  } else {ar1$ANh1<-(ncol(afh1)-ar1$NAh1)*4
  ploidy[i,2] <- 4}
  ar2$POP<-ar1$ACh1/ar1$ANh1
  colnames(ar2)[which(colnames(ar2) %in% "POP")]<-popnames[i]
}
ploidy<-ploidy[ order(ploidy[,1]), ]
#repolarise - #TODO
for (i in  1:nrow(ar2)){ # i=1
  d<-c()
  t<-c()
  for (pop in popnames) { #  pop="BAB"
    ppp <- ar2 %>% dplyr:: select(starts_with(pop))
    if (subset(popkey,popkey$Samples%like%pop)[1,3]==0) 
    {d<-c(d,as.numeric(as.character((ppp[i,1])[[1]])))} else {t<-c(t,as.numeric(as.character((ppp[i,1])[[1]])))}
  }
  meandip<-mean(d,na.rm = T)
  meantet<-mean(t,na.rm = T)
  if (meandip>meantet)
  {ar2[i,4:(length(ar2[i,])-1)]<-(1-ar2[i,4:(length(ar2[i,])-1)])
  print(i)
  } else {}}
#calculate haplotype allele frequencies
#df2<-read.table("AFsperLin.txt",h=T)
g<-read.table("PopGeogrAll10001.txt",h=T) ###Doplnit barvu, offset, repolarizovat??
g<-g[ order(g[,1]), ]
g$col <- ploidy$ploidy.lev
df2.transp1<-as.data.frame(t(ar2[,3:(ncol(ar2))]))
df2.transp<-df2.transp1[ order(row.names(df2.transp1)), ]
#reference
g$c1<-1-(apply(df2.transp, 1, max))
#alternative (selected)
g$c2<-apply(df2.transp, 1, min)
#recombined
g$oth<-abs(1-g$c1-g$c2)
g$nAA<-nrow(ar2)
g$length<-abs(max(ar2$start)-min(ar2$start))
g<-subset(g,!g$pop %in% "PAR" & !g$pop %in% "SUK")
write.table(g,paste(name,".hapl.arenosa.txt",sep=""),quote = F,sep = "\t",row.names = F)
#MAP
#list<-list.files(path = ".", pattern = "*.hapl.txt")
#hapl<-read.table(list,h=T)
hapl<-read.table(paste(name,".hapl.arenosa.txt",sep=""),h=T)
contr<-c("diploid","tetraploid","recombined")
colors<-c("tomato", "skyblue3", "grey")
### Central Europe
#for (l in list){ # l="AL8G38150_SMC6B.hapl.txt"
#  df3<-read.table(paste(l), header=T)
haplCeur<-subset(hapl,hapl$lon < 28 & hapl$lon > 2 & hapl$lat < 56)
lat <- c(min(haplCeur$lat)-1,max(haplCeur$lat)+1) #define our map's ylim
lon <- c(min(haplCeur$lon)-1,max(haplCeur$lon)+1) #define our map's xlim
#  name<-strsplit(l, split='[.]')[[1]][1]
# prot<-strsplit(name, split='[_]')[[1]][2]
terrmap<-GetMap.bbox(latR=c(min(haplCeur$lat)-1,max(haplCeur$lat)+1),lonR=c(min(haplCeur$lon)-1,max(haplCeur$lon)+1), maptype = "terrain", destfile = "terrain1.png",path = "&style=feature:all|element:labels|visibility:off", MINIMUMSIZE = T,GRAYSCALE=T)
size<-terrmap$size
pdf(file=paste(name,".CeurA.hapl.pdf", sep=""),width = size[1]/72, height = size[2]/72,pointsize = 10)
PlotOnStaticMap(terrmap)
df3<-na.omit(haplCeur)
for (i in 1:nrow(df3)){
  d1<-distance(df3$lat[i], df3$lat[i], lon[1], lon[2]) 
  d2<-distance(df3$lat[i], df3$lat[i], lon[1], df3$lon[i])
  a1<-size[1]/d1
  a2<-d2*a1
  a3<- -(size[1]/2)+a2+df3$offsetx[i]
  b1<-size[2]/(lat[2]-lat[1])
  b2<-(df3$lat[i]-lat[1])*b1
  b3<- -(size[2]/2)+b2+df3$offsety[i]
  terrmap<-add.pie(z=c(df3$c1[i],df3$c2[i],df3$oth[i]), x=a3, y=b3, radius=5.5, col = colors, labels=NA)
  terrmap<-text(x=a3, y=b3, labels = df3$pop[i], pos = df3$pos[i], col = df3$colorText[i], offset = 0.63,cex = 0.7)}
terrmap<-legend("topleft", inset=.005, title=paste(prot," haplotypes arenosa",sep=""),contr,fill = colors, horiz=F, cex=1,bg = 'white')
terrmap<-legend("bottomleft", inset=.005, legend = c(paste("Length = ",df3$length[1]," bp",sep=""),paste("Markers = ",df3$nAA[1]," AAs",sep="")), horiz=F, cex=1,bg = 'white')
dev.off()
### whole Europe
hapleur<-subset(hapl,hapl$lon < 70 & hapl$lon > -8)
lat <- c(min(hapleur$lat)-1,max(hapleur$lat)+1) #define our map's ylim
lon <- c(min(hapleur$lon)-1,max(hapleur$lon)+1) #define our map's xlim
terrmap<-GetMap.bbox(latR=c(min(hapleur$lat)-1,max(hapleur$lat)+1),lonR=c(min(hapleur$lon)-1,max(hapleur$lon)+1), maptype = "terrain", destfile = "terrain1.png",path = "&style=feature:all|element:labels|visibility:off", MINIMUMSIZE = T,GRAYSCALE=T)
size<-terrmap$size
pdf(file=paste(name,"EurA.hapl.pdf", sep=""),width = size[1]/72, height = size[2]/72,pointsize = 10)
PlotOnStaticMap(terrmap)
df3<-na.omit(hapleur)
for (i in 1:nrow(df3)){
  d1<-distance(df3$lat[i], df3$lat[i], lon[1], lon[2]) 
  d2<-distance(df3$lat[i], df3$lat[i], lon[1], df3$lon[i])
  a1<-size[1]/d1
  a2<-d2*a1
  a3<- -(size[1]/2)+a2+df3$offsetx[i]
  b1<-size[2]/(lat[2]-lat[1])
  b2<-(df3$lat[i]-lat[1])*b1
  b3<- -(size[2]/2)+b2+df3$offsety[i]
  terrmap<-add.pie(z=c(df3$c1[i],df3$c2[i],df3$oth[i]), x=a3, y=b3, radius=5.5, col = colors, labels=NA)
  terrmap<-text(x=a3, y=b3, labels = df3$pop[i], pos = df3$pos[i], col = df3$colorText[i], offset = 0.63,cex = 0.7)}
terrmap<-legend("bottomright", inset=.005, title=paste(prot," haplotypes arenosa",sep=""),contr,fill = colors, horiz=F, cex=1,bg = 'white')
terrmap<-legend("bottomleft", inset=.005, legend = c(paste("Length = ",df3$length[1]," bp",sep=""),paste("Markers = ",df3$nAA[1]," AAs",sep="")), horiz=F, cex=1,bg = 'white')
dev.off()
####Siberia
hapleur<-subset(hapl, hapl$lon > 70)
lat <- c(min(hapleur$lat)-1,max(hapleur$lat)+1) #define our map's ylim
lon <- c(min(hapleur$lon)-2,max(hapleur$lon)+4) #define our map's xlim
terrmap<-GetMap.bbox(latR=c(min(hapleur$lat),max(hapleur$lat)),lonR=c(min(hapleur$lon)-2,max(hapleur$lon)+4), maptype = "terrain", destfile = "terrain1.png",path = "&style=feature:all|element:labels|visibility:off", MINIMUMSIZE = T,GRAYSCALE=T)
size<-terrmap$size
pdf(file=paste(name,".SiberA.hapl.pdf", sep=""),width = size[1]/72, height = size[2]/72,pointsize = 10)
PlotOnStaticMap(terrmap)
df3<-na.omit(hapleur)
for (i in 1:nrow(df3)){
  d1<-distance(df3$lat[i], df3$lat[i], lon[1], lon[2]) 
  d2<-distance(df3$lat[i], df3$lat[i], lon[1], df3$lon[i])
  a1<-size[1]/d1
  a2<-d2*a1
  a3<- -(size[1]/2)+a2+df3$offsetx[i]
  b1<-size[2]/(lat[2]-lat[1])
  b2<-(df3$lat[i]-lat[1])*b1
  b3<- -(size[2]/2)+b2+df3$offsety[i]
  terrmap<-add.pie(z=c(df3$c1[i],df3$c2[i],df3$oth[i]), x=a3, y=b3, radius=5.5, col = colors, labels=NA)
  terrmap<-text(x=a3, y=b3, labels = df3$pop[i], pos = df3$pos[i], col = df3$colorText[i], offset = 0.63,cex = 0.7)}
terrmap<-legend("topright", inset=.005, title=paste(prot," arenosa",sep=""),contr,fill = colors, horiz=F, cex=0.7,bg = 'white')
dev.off()
####America
hapleur<-subset(hapl, hapl$lon < -7)
lat <- c(min(hapleur$lat)-1,max(hapleur$lat)+1) #define our map's ylim
lon <- c(min(hapleur$lon)-3,max(hapleur$lon)+3) #define our map's xlim
terrmap<-GetMap.bbox(latR=c(min(hapleur$lat)-1,max(hapleur$lat)+1),lonR=c(min(hapleur$lon)-3,max(hapleur$lon)+3), maptype = "terrain", destfile = "terrain1.png",path = "&style=feature:all|element:labels|visibility:off", MINIMUMSIZE = T,GRAYSCALE=T)
size<-terrmap$size
pdf(file=paste(name,".ZAmerA.hapl.pdf", sep=""),width = size[1]/72, height = size[2]/72,pointsize = 10)
PlotOnStaticMap(terrmap)
df3<-na.omit(hapleur)
for (i in 1:nrow(df3)){
  d1<-distance(df3$lat[i], df3$lat[i], lon[1], lon[2]) 
  d2<-distance(df3$lat[i], df3$lat[i], lon[1], df3$lon[i])
  a1<-size[1]/d1
  a2<-d2*a1
  a3<- -(size[1]/2)+a2+df3$offsetx[i]
  b1<-size[2]/(lat[2]-lat[1])
  b2<-(df3$lat[i]-lat[1])*b1
  b3<- -(size[2]/2)+b2+df3$offsety[i]
  terrmap<-add.pie(z=c(df3$c1[i],df3$c2[i],df3$oth[i]), x=a3, y=b3, radius=7, col = colors, labels=NA)
  terrmap<-text(x=a3, y=b3, labels = df3$pop[i], pos = df3$pos[i], col = df3$colorText[i], offset = 0.63,cex = 1)}
terrmap<-legend("topleft", inset=.005, title=paste(prot," arenosa",sep=""),contr,fill = colors, horiz=F, cex=0.7,bg = 'white')
dev.off()

####1b. get HISs - lyrata
setwd("/home/aa/lyrataWGD/fstScan/haplotypeMaps/")
hifs<-fread(paste(name,".HIS.lyrata.txt",sep = ""),h=T)
ar<-fread(paste('ALL.ann/ALL.table.recode.txt',sep=''),h=F,na.strings = "-9",nThread = 3)
popkey<-fread("PopKeyAll1000.csv",h=T)
prot <- name
names<-c(popkey[1:nrow(popkey),1])[[1]]
colnames(ar) <- c("pop","ploidy","scaff","start","AN","DP","ID","ann","AAS",names,"nan")
ar1<-subset(ar,ar$scaff %in% hifs$scaff & ar$start %in% hifs$start)
ar2<-as.data.frame(ar1[,c(3,4)])
popnames<-unique(substr(names,1,3))
popnames<-popnames[!popnames %in% "XXX"]
ploidy <- data.frame(pops = popnames, ploidy.lev = NA)
for(i in 1:length(popnames)){  # i = 109
  afh1 <- ar1 %>% dplyr:: select(starts_with(popnames[i]))
  ar1$ACh1 <- rowSums(afh1,na.rm = T)
  ar1$NAh1 <- apply(is.na(afh1), 1, sum)
  if (subset(popkey,popkey$Samples %like% popnames[i])[1,3]==0) 
  {ar1$ANh1<-(ncol(afh1)-ar1$NAh1)*2
  ploidy[i,2] <- 2
  } else {ar1$ANh1<-(ncol(afh1)-ar1$NAh1)*4
  ploidy[i,2] <- 4}
  ar2$POP<-ar1$ACh1/ar1$ANh1
  colnames(ar2)[which(colnames(ar2) %in% "POP")]<-popnames[i]
}
ploidy<-ploidy[ order(ploidy[,1]), ]
#repolarise - #TODO
for (i in  1:nrow(ar2)){ # i=1
  d<-c()
  t<-c()
  for (pop in popnames) { #  pop="BAB"
    ppp <- ar2 %>% dplyr:: select(starts_with(pop))
    if (subset(popkey,popkey$Samples%like%pop)[1,3]==0) 
    {d<-c(d,as.numeric(as.character((ppp[i,1])[[1]])))} else {t<-c(t,as.numeric(as.character((ppp[i,1])[[1]])))}
  }
  meandip<-mean(d,na.rm = T)
  meantet<-mean(t,na.rm = T)
  if (meandip>meantet)
  {ar2[i,4:(length(ar2[i,])-1)]<-(1-ar2[i,4:(length(ar2[i,])-1)])
  print(i)
  } else {}}
#calculate haplotype allele frequencies
#df2<-read.table("AFsperLin.txt",h=T)
g<-read.table("PopGeogrAll10001.txt",h=T) ###Doplnit barvu, offset, repolarizovat??
g<-g[ order(g[,1]), ]
g$col <- ploidy$ploidy.lev
df2.transp1<-as.data.frame(t(ar2[,3:(ncol(ar2))]))
df2.transp<-df2.transp1[ order(row.names(df2.transp1)), ]
#reference
g$c1<-1-(apply(df2.transp, 1, max))
#alternative (selected)
g$c2<-apply(df2.transp, 1, min)
#recombined
g$oth<-abs(1-g$c1-g$c2)
g$nAA<-nrow(ar2)
g$length<-abs(max(ar2$start)-min(ar2$start))
g<-subset(g,!g$pop %in% "PAR" & !g$pop %in% "SUK")
write.table(g,paste(name,".hapl.lyrata.txt",sep=""),quote = F,sep = "\t",row.names = F)
#MAP
#list<-list.files(path = ".", pattern = "*.hapl.txt")
#hapl<-read.table(list,h=T)
hapl<-read.table(paste(name,".hapl.lyrata.txt",sep=""),h=T)
contr<-c("diploid","tetraploid","recombined")
colors<-c("red", "blue", "grey")
### Central Europe
#for (l in list){ # l="AL8G38150_SMC6B.hapl.txt"
#  df3<-read.table(paste(l), header=T)
haplCeur<-subset(hapl,hapl$lon < 28 & hapl$lon > 2 & hapl$lat < 56)
lat <- c(min(haplCeur$lat)-1,max(haplCeur$lat)+1) #define our map's ylim
lon <- c(min(haplCeur$lon)-1,max(haplCeur$lon)+1) #define our map's xlim
#  name<-strsplit(l, split='[.]')[[1]][1]
# prot<-strsplit(name, split='[_]')[[1]][2]
terrmap<-GetMap.bbox(latR=c(min(haplCeur$lat)-1,max(haplCeur$lat)+1),lonR=c(min(haplCeur$lon)-1,max(haplCeur$lon)+1), maptype = "terrain", destfile = "terrain1.png",path = "&style=feature:all|element:labels|visibility:off", MINIMUMSIZE = T,GRAYSCALE=T)
size<-terrmap$size
pdf(file=paste(name,".CeurL.hapl.pdf", sep=""),width = size[1]/72, height = size[2]/72,pointsize = 10)
PlotOnStaticMap(terrmap)
df3<-na.omit(haplCeur)
for (i in 1:nrow(df3)){
  d1<-distance(df3$lat[i], df3$lat[i], lon[1], lon[2]) 
  d2<-distance(df3$lat[i], df3$lat[i], lon[1], df3$lon[i])
  a1<-size[1]/d1
  a2<-d2*a1
  a3<- -(size[1]/2)+a2+df3$offsetx[i]
  b1<-size[2]/(lat[2]-lat[1])
  b2<-(df3$lat[i]-lat[1])*b1
  b3<- -(size[2]/2)+b2+df3$offsety[i]
  terrmap<-add.pie(z=c(df3$c1[i],df3$c2[i],df3$oth[i]), x=a3, y=b3, radius=5.5, col = colors, labels=NA)
  terrmap<-text(x=a3, y=b3, labels = df3$pop[i], pos = df3$pos[i], col = df3$colorText[i], offset = 0.63,cex = 0.7)}
terrmap<-legend("topleft", inset=.005, title=paste(prot," haplotypes lyrata",sep=""),contr,fill = colors, horiz=F, cex=1,bg = 'white')
terrmap<-legend("bottomleft", inset=.005, legend = c(paste("Length = ",df3$length[1]," bp",sep=""),paste("Markers = ",df3$nAA[1]," AAs",sep="")), horiz=F, cex=1,bg = 'white')
dev.off()
### whole Europe
hapleur<-subset(hapl,hapl$lon < 70 & hapl$lon > -8)
lat <- c(min(hapleur$lat)-1,max(hapleur$lat)+1) #define our map's ylim
lon <- c(min(hapleur$lon)-1,max(hapleur$lon)+1) #define our map's xlim
terrmap<-GetMap.bbox(latR=c(min(hapleur$lat)-1,max(hapleur$lat)+1),lonR=c(min(hapleur$lon)-1,max(hapleur$lon)+1), maptype = "terrain", destfile = "terrain1.png",path = "&style=feature:all|element:labels|visibility:off", MINIMUMSIZE = T,GRAYSCALE=T)
size<-terrmap$size
pdf(file=paste(name,"EurL.hapl.pdf", sep=""),width = size[1]/72, height = size[2]/72,pointsize = 10)
PlotOnStaticMap(terrmap)
df3<-na.omit(hapleur)
for (i in 1:nrow(df3)){
  d1<-distance(df3$lat[i], df3$lat[i], lon[1], lon[2]) 
  d2<-distance(df3$lat[i], df3$lat[i], lon[1], df3$lon[i])
  a1<-size[1]/d1
  a2<-d2*a1
  a3<- -(size[1]/2)+a2+df3$offsetx[i]
  b1<-size[2]/(lat[2]-lat[1])
  b2<-(df3$lat[i]-lat[1])*b1
  b3<- -(size[2]/2)+b2+df3$offsety[i]
  terrmap<-add.pie(z=c(df3$c1[i],df3$c2[i],df3$oth[i]), x=a3, y=b3, radius=5.5, col = colors, labels=NA)
  terrmap<-text(x=a3, y=b3, labels = df3$pop[i], pos = df3$pos[i], col = df3$colorText[i], offset = 0.63,cex = 0.7)}
terrmap<-legend("topleft", inset=.005, title=paste(prot," haplotypes lyrata",sep=""),contr,fill = colors, horiz=F, cex=1,bg = 'white')
terrmap<-legend("bottomleft", inset=.005, legend = c(paste("Length = ",df3$length[1]," bp",sep=""),paste("Markers = ",df3$nAA[1]," AAs",sep="")), horiz=F, cex=1,bg = 'white')
dev.off()
####Siberia
hapleur<-subset(hapl, hapl$lon > 70)
lat <- c(min(hapleur$lat)-1,max(hapleur$lat)+1) #define our map's ylim
lon <- c(min(hapleur$lon)-2,max(hapleur$lon)+4) #define our map's xlim
terrmap<-GetMap.bbox(latR=c(min(hapleur$lat),max(hapleur$lat)),lonR=c(min(hapleur$lon)-2,max(hapleur$lon)+4), maptype = "terrain", destfile = "terrain1.png",path = "&style=feature:all|element:labels|visibility:off", MINIMUMSIZE = T,GRAYSCALE=T)
size<-terrmap$size
pdf(file=paste(name,".SiberL.hapl.pdf", sep=""),width = size[1]/72, height = size[2]/72,pointsize = 10)
PlotOnStaticMap(terrmap)
df3<-na.omit(hapleur)
for (i in 1:nrow(df3)){
  d1<-distance(df3$lat[i], df3$lat[i], lon[1], lon[2]) 
  d2<-distance(df3$lat[i], df3$lat[i], lon[1], df3$lon[i])
  a1<-size[1]/d1
  a2<-d2*a1
  a3<- -(size[1]/2)+a2+df3$offsetx[i]
  b1<-size[2]/(lat[2]-lat[1])
  b2<-(df3$lat[i]-lat[1])*b1
  b3<- -(size[2]/2)+b2+df3$offsety[i]
  terrmap<-add.pie(z=c(df3$c1[i],df3$c2[i],df3$oth[i]), x=a3, y=b3, radius=5.5, col = colors, labels=NA)
  terrmap<-text(x=a3, y=b3, labels = df3$pop[i], pos = df3$pos[i], col = df3$colorText[i], offset = 0.63,cex = 0.7)}
terrmap<-legend("bottomright", inset=.005, title=paste(prot," lyrata",sep=""),contr,fill = colors, horiz=F, cex=0.7,bg = 'white')
dev.off()
####America
hapleur<-subset(hapl, hapl$lon < -7)
lat <- c(min(hapleur$lat)-1,max(hapleur$lat)+1) #define our map's ylim
lon <- c(min(hapleur$lon)-3,max(hapleur$lon)+3) #define our map's xlim
terrmap<-GetMap.bbox(latR=c(min(hapleur$lat)-1,max(hapleur$lat)+1),lonR=c(min(hapleur$lon)-3,max(hapleur$lon)+3), maptype = "terrain", destfile = "terrain1.png",path = "&style=feature:all|element:labels|visibility:off", MINIMUMSIZE = T,GRAYSCALE=T)
size<-terrmap$size
pdf(file=paste(name,".ZAmerL.hapl.pdf", sep=""),width = size[1]/72, height = size[2]/72,pointsize = 10)
PlotOnStaticMap(terrmap)
df3<-na.omit(hapleur)
for (i in 1:nrow(df3)){
  d1<-distance(df3$lat[i], df3$lat[i], lon[1], lon[2]) 
  d2<-distance(df3$lat[i], df3$lat[i], lon[1], df3$lon[i])
  a1<-size[1]/d1
  a2<-d2*a1
  a3<- -(size[1]/2)+a2+df3$offsetx[i]
  b1<-size[2]/(lat[2]-lat[1])
  b2<-(df3$lat[i]-lat[1])*b1
  b3<- -(size[2]/2)+b2+df3$offsety[i]
  terrmap<-add.pie(z=c(df3$c1[i],df3$c2[i],df3$oth[i]), x=a3, y=b3, radius=7, col = colors, labels=NA)
  terrmap<-text(x=a3, y=b3, labels = df3$pop[i], pos = df3$pos[i], col = df3$colorText[i], offset = 0.63,cex = 1)}
terrmap<-legend("topleft", inset=.005, title=paste(prot," lyrata",sep=""),contr,fill = colors, horiz=F, cex=0.7,bg = 'white')
dev.off()

####1c. get HISs - union
setwd("/home/aa/lyrataWGD/fstScan/haplotypeMaps/")
hifs<-fread(paste(name,".HIS.union.txt",sep = ""),h=T)
ar<-fread(paste('ALL.ann/ALL.table.recode.txt',sep=''),h=F,na.strings = "-9",nThread = 3)
popkey<-fread("PopKeyAll1000.csv",h=T)
prot <- name
names<-c(popkey[1:nrow(popkey),1])[[1]]
colnames(ar) <- c("pop","ploidy","scaff","start","AN","DP","ID","ann","AAS",names,"nan")
ar1<-subset(ar,ar$scaff %in% hifs$scaff & ar$start %in% hifs$start)
ar2<-as.data.frame(ar1[,c(3,4)])
popnames<-unique(substr(names,1,3))
popnames<-popnames[!popnames %in% "XXX"]
ploidy <- data.frame(pops = popnames, ploidy.lev = NA)
for(i in 1:length(popnames)){  # i = 121
  afh1 <- ar1 %>% dplyr:: select(starts_with(popnames[i]))
  ar1$ACh1 <- rowSums(afh1,na.rm = T)
  ar1$NAh1 <- apply(is.na(afh1), 1, sum)
  if (subset(popkey,popkey$Samples %like% popnames[i])[1,3]==0) 
  {ar1$ANh1<-(ncol(afh1)-ar1$NAh1)*2
  ploidy[i,2] <- 2
  } else {ar1$ANh1<-(ncol(afh1)-ar1$NAh1)*4
  ploidy[i,2] <- 4}
  ar2$POP<-ar1$ACh1/ar1$ANh1
  colnames(ar2)[which(colnames(ar2) %in% "POP")]<-popnames[i]
}
ploidy<-ploidy[ order(ploidy[,1]), ]
#repolarise - #TODO
for (i in  1:nrow(ar2)){ # i=1
  d<-c()
  t<-c()
  for (pop in popnames) { #  pop="BAB"
    ppp <- ar2 %>% dplyr:: select(starts_with(pop))
    if (subset(popkey,popkey$Samples%like%pop)[1,3]==0) 
    {d<-c(d,as.numeric(as.character((ppp[i,1])[[1]])))} else {t<-c(t,as.numeric(as.character((ppp[i,1])[[1]])))}
  }
  meandip<-mean(d,na.rm = T)
  meantet<-mean(t,na.rm = T)
  if (meandip>meantet)
  {ar2[i,4:(length(ar2[i,])-1)]<-(1-ar2[i,4:(length(ar2[i,])-1)])
  print(i)
  } else {}}
#calculate haplotype allele frequencies
#df2<-read.table("AFsperLin.txt",h=T)
g<-read.table("PopGeogrAll10001.txt",h=T) ###Doplnit barvu, offset, repolarizovat??
g<-g[ order(g[,1]), ]
g$col <- ploidy$ploidy.lev
df2.transp1<-as.data.frame(t(ar2[,3:(ncol(ar2))]))
df2.transp<-df2.transp1[ order(row.names(df2.transp1)), ]
#reference
g$c1<-1-(apply(df2.transp, 1, max))
#alternative (selected)
g$c2<-apply(df2.transp, 1, min)
#recombined
g$oth<-abs(1-g$c1-g$c2)
g$nAA<-nrow(ar2)
g$length<-abs(max(ar2$start)-min(ar2$start))
g<-subset(g,!g$pop %in% "PAR" & !g$pop %in% "SUK")
write.table(g,paste(name,".hapl.union.txt",sep=""),quote = F,sep = "\t",row.names = F)
#MAP
#list<-list.files(path = ".", pattern = "*.hapl.txt")
#hapl<-read.table(list,h=T)
hapl<-read.table(paste(name,".hapl.union.txt",sep=""),h=T)
contr<-c("diploid","tetraploid","recombined")
colors<-c("red3", "blue3", "grey")
### Central Europe
#for (l in list){ # l="AL8G38150_SMC6B.hapl.txt"
#  df3<-read.table(paste(l), header=T)
haplCeur<-subset(hapl,hapl$lon < 28 & hapl$lon > 2 & hapl$lat < 56)
lat <- c(min(haplCeur$lat)-1,max(haplCeur$lat)+1) #define our map's ylim
lon <- c(min(haplCeur$lon)-1,max(haplCeur$lon)+1) #define our map's xlim
#  name<-strsplit(l, split='[.]')[[1]][1]
# prot<-strsplit(name, split='[_]')[[1]][2]
terrmap<-GetMap.bbox(latR=c(min(haplCeur$lat)-1,max(haplCeur$lat)+1),lonR=c(min(haplCeur$lon)-1,max(haplCeur$lon)+1), maptype = "terrain", destfile = "terrain1.png",path = "&style=feature:all|element:labels|visibility:off", MINIMUMSIZE = T,GRAYSCALE=T)
size<-terrmap$size
pdf(file=paste(name,".CeurS.hapl.pdf", sep=""),width = size[1]/72, height = size[2]/72,pointsize = 10)
PlotOnStaticMap(terrmap)
df3<-na.omit(haplCeur)
for (i in 1:nrow(df3)){
  d1<-distance(df3$lat[i], df3$lat[i], lon[1], lon[2]) 
  d2<-distance(df3$lat[i], df3$lat[i], lon[1], df3$lon[i])
  a1<-size[1]/d1
  a2<-d2*a1
  a3<- -(size[1]/2)+a2+df3$offsetx[i]
  b1<-size[2]/(lat[2]-lat[1])
  b2<-(df3$lat[i]-lat[1])*b1
  b3<- -(size[2]/2)+b2+df3$offsety[i]
  terrmap<-add.pie(z=c(df3$c1[i],df3$c2[i],df3$oth[i]), x=a3, y=b3, radius=5.5, col = colors, labels=NA)
  terrmap<-text(x=a3, y=b3, labels = df3$pop[i], pos = df3$pos[i], col = df3$colorText[i], offset = 0.63,cex = 0.7)}
terrmap<-legend("topleft", inset=.005, title=paste(prot," haplotypes - union",sep=""),contr,fill = colors, horiz=F, cex=1,bg = 'white')
terrmap<-legend("bottomleft", inset=.005, legend = c(paste("Length = ",df3$length[1]," bp",sep=""),paste("Markers = ",df3$nAA[1]," AAs",sep="")), horiz=F, cex=1,bg = 'white')
dev.off()
### whole Europe
hapleur<-subset(hapl,hapl$lon < 70 & hapl$lon > -8)
lat <- c(min(hapleur$lat)-1,max(hapleur$lat)+1) #define our map's ylim
lon <- c(min(hapleur$lon)-1,max(hapleur$lon)+1) #define our map's xlim
terrmap<-GetMap.bbox(latR=c(min(hapleur$lat)-1,max(hapleur$lat)+1),lonR=c(min(hapleur$lon)-1,max(hapleur$lon)+1), maptype = "terrain", destfile = "terrain1.png",path = "&style=feature:all|element:labels|visibility:off", MINIMUMSIZE = T,GRAYSCALE=T)
size<-terrmap$size
pdf(file=paste(name,"EurS.hapl.pdf", sep=""),width = size[1]/72, height = size[2]/72,pointsize = 10)
PlotOnStaticMap(terrmap)
df3<-na.omit(hapleur)
for (i in 1:nrow(df3)){
  d1<-distance(df3$lat[i], df3$lat[i], lon[1], lon[2]) 
  d2<-distance(df3$lat[i], df3$lat[i], lon[1], df3$lon[i])
  a1<-size[1]/d1
  a2<-d2*a1
  a3<- -(size[1]/2)+a2+df3$offsetx[i]
  b1<-size[2]/(lat[2]-lat[1])
  b2<-(df3$lat[i]-lat[1])*b1
  b3<- -(size[2]/2)+b2+df3$offsety[i]
  terrmap<-add.pie(z=c(df3$c1[i],df3$c2[i],df3$oth[i]), x=a3, y=b3, radius=5.5, col = colors, labels=NA)
  terrmap<-text(x=a3, y=b3, labels = df3$pop[i], pos = df3$pos[i], col = df3$colorText[i], offset = 0.63,cex = 0.7)}
terrmap<-legend("topleft", inset=.005, title=paste(prot," haplotypes - union",sep=""),contr,fill = colors, horiz=F, cex=1,bg = 'white')
terrmap<-legend("bottomleft", inset=.005, legend = c(paste("Length = ",df3$length[1]," bp",sep=""),paste("Markers = ",df3$nAA[1]," AAs",sep="")), horiz=F, cex=1,bg = 'white')
dev.off()
####Siberia
hapleur<-subset(hapl, hapl$lon > 70)
lat <- c(min(hapleur$lat)-1,max(hapleur$lat)+1) #define our map's ylim
lon <- c(min(hapleur$lon)-2,max(hapleur$lon)+4) #define our map's xlim
terrmap<-GetMap.bbox(latR=c(min(hapleur$lat),max(hapleur$lat)),lonR=c(min(hapleur$lon)-2,max(hapleur$lon)+4), maptype = "terrain", destfile = "terrain1.png",path = "&style=feature:all|element:labels|visibility:off", MINIMUMSIZE = T,GRAYSCALE=T)
size<-terrmap$size
pdf(file=paste(name,".SiberS.hapl.pdf", sep=""),width = size[1]/72, height = size[2]/72,pointsize = 10)
PlotOnStaticMap(terrmap)
df3<-na.omit(hapleur)
for (i in 1:nrow(df3)){
  d1<-distance(df3$lat[i], df3$lat[i], lon[1], lon[2]) 
  d2<-distance(df3$lat[i], df3$lat[i], lon[1], df3$lon[i])
  a1<-size[1]/d1
  a2<-d2*a1
  a3<- -(size[1]/2)+a2+df3$offsetx[i]
  b1<-size[2]/(lat[2]-lat[1])
  b2<-(df3$lat[i]-lat[1])*b1
  b3<- -(size[2]/2)+b2+df3$offsety[i]
  terrmap<-add.pie(z=c(df3$c1[i],df3$c2[i],df3$oth[i]), x=a3, y=b3, radius=5.5, col = colors, labels=NA)
  terrmap<-text(x=a3, y=b3, labels = df3$pop[i], pos = df3$pos[i], col = df3$colorText[i], offset = 0.63,cex = 0.7)}
terrmap<-legend("bottomright", inset=.005, title=paste(prot," union",sep=""),contr,fill = colors, horiz=F, cex=0.7,bg = 'white')
dev.off()
####America
hapleur<-subset(hapl, hapl$lon < -7)
lat <- c(min(hapleur$lat)-1,max(hapleur$lat)+1) #define our map's ylim
lon <- c(min(hapleur$lon)-3,max(hapleur$lon)+3) #define our map's xlim
terrmap<-GetMap.bbox(latR=c(min(hapleur$lat)-1,max(hapleur$lat)+1),lonR=c(min(hapleur$lon)-3,max(hapleur$lon)+3), maptype = "terrain", destfile = "terrain1.png",path = "&style=feature:all|element:labels|visibility:off", MINIMUMSIZE = T,GRAYSCALE=T)
size<-terrmap$size
pdf(file=paste(name,".ZAmerS.hapl.pdf", sep=""),width = size[1]/72, height = size[2]/72,pointsize = 10)
PlotOnStaticMap(terrmap)
df3<-na.omit(hapleur)
for (i in 1:nrow(df3)){
  d1<-distance(df3$lat[i], df3$lat[i], lon[1], lon[2]) 
  d2<-distance(df3$lat[i], df3$lat[i], lon[1], df3$lon[i])
  a1<-size[1]/d1
  a2<-d2*a1
  a3<- -(size[1]/2)+a2+df3$offsetx[i]
  b1<-size[2]/(lat[2]-lat[1])
  b2<-(df3$lat[i]-lat[1])*b1
  b3<- -(size[2]/2)+b2+df3$offsety[i]
  terrmap<-add.pie(z=c(df3$c1[i],df3$c2[i],df3$oth[i]), x=a3, y=b3, radius=7, col = colors, labels=NA)
  terrmap<-text(x=a3, y=b3, labels = df3$pop[i], pos = df3$pos[i], col = df3$colorText[i], offset = 0.63,cex = 1)}
terrmap<-legend("topleft", inset=.005, title=paste(prot," union",sep=""),contr,fill = colors, horiz=F, cex=0.7,bg = 'white')
dev.off()
###combine together
infiles <- Sys.glob("*.hapl.pdf") 
outfile <- paste("haplotypes.",name,".pdf",sep="") 
system(paste("pdftk", paste(infiles, collapse = " "), "cat output", outfile))  
junk <- dir(path=".",  pattern="*.hapl.pdf") 
file.remove(junk) }
#############THE END  
#})}
  
######### trace HIS through phypogeny
#PRIOR: blast gene, take 15 best hits of silimar length, add to cart, download fasta, geneious aligh + set reference + delete indels, save as alignment fasta

library(seqinr)
library(readr)
library(data.table)
library(dplyr)
library(gplots)

setwd("/home/aa/lyrataWGD/fstScan/ProtHistory/")
name<-"CYCD32"
aa<-read.table(paste("../haplotypeMaps/",name,".HIS.arenosa.txt",sep = ""),h=T,sep="\t")
al<-read.table(paste("../haplotypeMaps/",name,".HIS.lyrata.txt",sep = ""),h=T,sep="\t")
both<-rbind(aa[,1:5],al[,1:5])
both<-both[ order(both[,2]), ]
bo<-both[!duplicated(both[,c('scaff','start')]),]
bo$pos<-parse_number(bo$aas)
aa<-read.alignment(file = paste("/home/aa/lyrataWGD/fstScan/ProtHistory/",name,"alignment.fasta",sep=""),format = "fasta")
aaa<-as.matrix(aa$seq)
aaa1<-as.data.frame(aaa)
aaa2<-as.data.frame(strsplit(as.character(aa$seq), ""))
colnames(aaa2)<-substr(aa$nam,1,10)
aaa3<-aaa2[bo$pos,]
fin<-cbind(bo,aaa3)
write.table(fin,paste(name,".HIS.brassicaceae.txt",sep=""),row.names = F,quote = F,sep = "\t")




###Heatmap candidate only
library(data.table)
library(dplyr)
library(gplots)

setwd("/home/aa/lyrataWGD/fstScan/ProtHistory/")
allall<-fread(paste('../heatmap/ann/ALL.table.recode.txt',sep=''),h=F,na.strings = "-9",nThread = 3)
nnn<-fread("../heatmap/PopKey_lyrenosa.heatmap.csv",h=F)
names<-c(nnn[2:nrow(nnn),1])[[1]]
colnames(allall) <- c("pop","ploidy","scaff","start","AN","DP","ID","ann","aas",names,"nan")
namesAll<-read.table("file:///home/aa/lyrataWGD/fstScan/haplotypeMaps/17WGDcand_names.txt",h=T)
for (i in 1:nrow(namesAll)) { #  i=1
  i2<-namesAll[i,1]
  name<-namesAll[i,2]
  bo<-fread(paste("../haplotypeMaps/",name,".HIS.union.txt",sep = ""),h=T)
ar1<-subset(allall,allall$ID %in% bo$ID & allall$start %in% bo$start)[,1:(ncol(allall)-1)]
#2. filter too low freq
i<-which(colnames(ar1) %like% "_")[1]
y<-which(colnames(ar1) %like% "_")[length(which(colnames(ar1) %like% "_"))]
afh1<-ar1[,i:as.numeric(y)]
ar1$ACh1<-rowSums(afh1,na.rm = T)
ar1$NAh1<-apply(is.na(afh1), 1, sum)
ar1$ANh1<-486 ####hard calculated!
ar1$tot<-ar1$ACh1/(ar1$ANh1-ar1$NAh1)
ar2<-subset(ar1,ar1$tot > 16/486)
ar3<-subset(ar2,ar2$tot < as.numeric(1-(16/486)))
#3. calculate AF for each lineage TODO 
popnames<-unique(substr(names,1,3))
for (pop in popnames) { #  pop="BAB"
  afh1 <- ar3 %>% dplyr:: select(starts_with(pop))
  ar3$ACh1<-rowSums(afh1,na.rm = T)
  ar3$NAh1<-apply(is.na(afh1), 1, sum)
  if (subset(nnn,nnn$V1%like%pop)[1,3]==0) 
  {ar3$ANh1<-(ncol(afh1)-ar3$NAh1)*2} else {ar3$ANh1<-(ncol(afh1)-ar3$NAh1)*4}
  ar3$POP<-ar3$ACh1/ar3$ANh1
  colnames(ar3)[which(colnames(ar3) %in% "POP")]<-pop
}
#     write.table(ar3,"CAPD3.txt",row.names = F,quote = F,sep="\t")
s<-dplyr::select(ar3,ID,ann,aas,ZEP, BAB, SUB, HOC, PER, BRE, VLA, CHO, SWA, JOH, MOD,TEM,PEK,SCT,VLH,OSL,STD,tot)
#repolarise
popnames<-colnames(s)[4:(length(colnames(s))-1)]
for (i in  1:nrow(s)){ # i=1
  d<-c()
  t<-c()
  for (pop in popnames) { #  pop="BAB"
    ppp <- s %>% dplyr:: select(starts_with(pop))
    if (subset(nnn,nnn$V1%like%pop)[1,3]==0) 
    {d<-c(d,as.numeric(as.character((ppp[i,1])[[1]])))} else {t<-c(t,as.numeric(as.character((ppp[i,1])[[1]])))}
  }
  meandip<-mean(d,na.rm = T)
  meantet<-mean(t,na.rm = T)
  if (meandip>meantet)
  {s[i,4:(length(s[i,])-1)]<-(1-s[i,4:(length(s[i,])-1)])
  print(i)
  } else {}}
#plot
ann<-fread("/home/aa/Desktop/references/lyrataV2/LyV2_TAIR11orth_des_20171231.txt")
pops<-c(colnames(s)[4:(length(colnames(s))-1)])
my_palette <- colorRampPalette(c("tomato", "blue"))(n = 100)
id<-ar1$ID[1]
  if (nrow(subset(s,s$ID %in% id))>1)
  {s1<-subset(s,s$ID %in% id)
  ann1<-subset(ann,ann$AL %in% id)
  if (nrow(s1)>110)
  {p=0.75
  } else {p=1}
  s5<-dplyr::select(s1,SUB, HOC, PER, JOH, MOD, VLH, BAB, BRE, VLA, TEM, PEK, OSL, ZEP, CHO, SWA, SCT, STD)
  df<-as.matrix(s1[,4:(length(s1[i,])-1)],rownames = s1$aas)
  df1<-as.matrix(s5,rownames = s1$aas)
  #parallel
  pdf(paste("heatmaps/heatmap_",id,"_",name,".pdf",sep=""),height = 11,width = 9.3)
  heatmap.2(x = df,dendrogram = "none",Colv="NA", Rowv="NA",key = F,col=my_palette,colsep= c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,14,17,17,17,17,17,17,17,17,17),sepcolor= c("black"),sepwidth = c(0.02),trace="none",ColSideColors = c("red","red","red","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","red","red","red"),labCol=pops,colCol= c("tomato","tomato","tomato","skyblue2","skyblue2","skyblue2","skyblue2","skyblue2","skyblue2","blue","blue","blue","blue","blue","red","red","red"),offsetRow = c(0.05),offsetCol= c(0.05),cexRow = c(p),cexCol = c(2),margins = c(7,7),lhei = c(0.2,9),lwid = c(0.03,0.8),xlab = paste(id,ann1$AT,ann1$annShort,sep=" - ")) 
  #independent
  heatmap.2(x = df1,dendrogram = "none",Colv="NA", Rowv="NA",key = F,col=my_palette,colsep= c(0,0,6,12,17),sepcolor= c("black"),sepwidth = c(0.02),trace="none",ColSideColors = c("red","blue","blue","blue","blue","red","red","blue","blue","blue","blue","red","red","blue","blue","blue","red"),labCol=colnames(s5),colCol= c("tomato","skyblue2","skyblue2","blue","blue","red","tomato","skyblue2","skyblue2","blue","blue","red","tomato","blue","blue",'skyblue2',"red"),offsetRow = c(0.05),offsetCol= c(0.05),cexRow = c(p),cexCol = c(2),margins = c(7,7),lhei = c(0.2,9),lwid = c(0.03,0.8),xlab = paste(id,ann1$AT,ann1$annShort,sep=" - ")) 
  dev.off()
  } else {} 
}


####Europe and 
###Heatmap candidate only
library(data.table)
library(dplyr)
library(gplots)

setwd("/home/aa/lyrataWGD/fstScan/ProtHistory/")
allall<-fread(paste('../heatmap/ann/ALL.table.recode.txt',sep=''),h=F,na.strings = "-9",nThread = 3)
nnn<-fread("../heatmap/PopKey_lyrenosa.heatmap.csv",h=F)
names<-c(nnn[2:nrow(nnn),1])[[1]]
colnames(allall) <- c("pop","ploidy","scaff","start","AN","DP","ID","ann","aas",names,"nan")
namesAll<-read.table("file:///home/aa/lyrataWGD/fstScan/haplotypeMaps/17WGDcand_names.txt",h=T)
for (i in 1:nrow(namesAll)) { #  i=1
  i2<-namesAll[i,1]
  name<-namesAll[i,2]
  bo<-fread(paste("../haplotypeMaps/",name,".HIS.union.txt",sep = ""),h=T)
  ar1<-subset(allall,allall$ID %in% bo$ID & allall$start %in% bo$start)[,1:(ncol(allall)-1)]
  #2. filter too low freq
  i<-which(colnames(ar1) %like% "_")[1]
  y<-which(colnames(ar1) %like% "_")[length(which(colnames(ar1) %like% "_"))]
  afh1<-ar1[,i:as.numeric(y)]
  ar1$ACh1<-rowSums(afh1,na.rm = T)
  ar1$NAh1<-apply(is.na(afh1), 1, sum)
  ar1$ANh1<-486 ####hard calculated!
  ar1$tot<-ar1$ACh1/(ar1$ANh1-ar1$NAh1)
  ar2<-subset(ar1,ar1$tot > 16/486)
  ar3<-subset(ar2,ar2$tot < as.numeric(1-(16/486)))
  #3. calculate AF for each lineage TODO 
  popnames<-unique(substr(names,1,3))
  for (pop in popnames) { #  pop="BAB"
    afh1 <- ar3 %>% dplyr:: select(starts_with(pop))
    ar3$ACh1<-rowSums(afh1,na.rm = T)
    ar3$NAh1<-apply(is.na(afh1), 1, sum)
    if (subset(nnn,nnn$V1%like%pop)[1,3]==0) 
    {ar3$ANh1<-(ncol(afh1)-ar3$NAh1)*2} else {ar3$ANh1<-(ncol(afh1)-ar3$NAh1)*4}
    ar3$POP<-ar3$ACh1/ar3$ANh1
    colnames(ar3)[which(colnames(ar3) %in% "POP")]<-pop
  }
  #     write.table(ar3,"CAPD3.txt",row.names = F,quote = F,sep="\t")
  s<-dplyr::select(ar3,ID,ann,aas,ZEP, BAB, SUB, HOC, PER, BRE, VLA, CHO, SWA, JOH, MOD,TEM,PEK,SCT,VLH,OSL,STD,tot)
  #repolarise
  popnames<-colnames(s)[4:(length(colnames(s))-1)]
  for (i in  1:nrow(s)){ # i=1
    d<-c()
    t<-c()
    for (pop in popnames) { #  pop="BAB"
      ppp <- s %>% dplyr:: select(starts_with(pop))
      if (subset(nnn,nnn$V1%like%pop)[1,3]==0) 
      {d<-c(d,as.numeric(as.character((ppp[i,1])[[1]])))} else {t<-c(t,as.numeric(as.character((ppp[i,1])[[1]])))}
    }
    meandip<-mean(d,na.rm = T)
    meantet<-mean(t,na.rm = T)
    if (meandip>meantet)
    {s[i,4:(length(s[i,])-1)]<-(1-s[i,4:(length(s[i,])-1)])
    print(i)
    } else {}}
  #plot
  ann<-fread("/home/aa/Desktop/references/lyrataV2/LyV2_TAIR11orth_des_20171231.txt")
  pops<-c(colnames(s)[4:(length(colnames(s))-1)])
  my_palette <- colorRampPalette(c("tomato", "blue"))(n = 100)
  id<-ar1$ID[1]
  if (nrow(subset(s,s$ID %in% id))>1)
  {s1<-subset(s,s$ID %in% id)
  ann1<-subset(ann,ann$AL %in% id)
  if (nrow(s1)>110)
  {p=0.75
  } else {p=1}
  s5<-dplyr::select(s1,SUB, HOC, PER, JOH, MOD, VLH, BAB, BRE, VLA, TEM, PEK, OSL, ZEP, CHO, SWA, SCT, STD)
  df<-as.matrix(s1[,4:(length(s1[i,])-1)],rownames = s1$aas)
  df1<-as.matrix(s5,rownames = s1$aas)
  #parallel
  pdf(paste("heatmaps/heatmap_",id,"_",name,".pdf",sep=""),height = 11,width = 9.3)
  heatmap.2(x = df,dendrogram = "none",Colv="NA", Rowv="NA",key = F,col=my_palette,colsep= c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,14,17,17,17,17,17,17,17,17,17),sepcolor= c("black"),sepwidth = c(0.02),trace="none",ColSideColors = c("red","red","red","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","red","red","red"),labCol=pops,colCol= c("tomato","tomato","tomato","skyblue2","skyblue2","skyblue2","skyblue2","skyblue2","skyblue2","blue","blue","blue","blue","blue","red","red","red"),offsetRow = c(0.05),offsetCol= c(0.05),cexRow = c(p),cexCol = c(2),margins = c(7,7),lhei = c(0.2,9),lwid = c(0.03,0.8),xlab = paste(id,ann1$AT,ann1$annShort,sep=" - ")) 
  #independent
  heatmap.2(x = df1,dendrogram = "none",Colv="NA", Rowv="NA",key = F,col=my_palette,colsep= c(0,0,6,12,17),sepcolor= c("black"),sepwidth = c(0.02),trace="none",ColSideColors = c("red","blue","blue","blue","blue","red","red","blue","blue","blue","blue","red","red","blue","blue","blue","red"),labCol=colnames(s5),colCol= c("tomato","skyblue2","skyblue2","blue","blue","red","tomato","skyblue2","skyblue2","blue","blue","red","tomato","blue","blue",'skyblue2',"red"),offsetRow = c(0.05),offsetCol= c(0.05),cexRow = c(p),cexCol = c(2),margins = c(7,7),lhei = c(0.2,9),lwid = c(0.03,0.8),xlab = paste(id,ann1$AT,ann1$annShort,sep=" - ")) 
  dev.off()
  } else {} 
  
}





  
####Interval list for PCA etc
library(data.table)
setwd<-("/home/aa/lyrataWGD/fstScan/haplotypeMaps/")
namesAll<-read.table("file:///home/aa/lyrataWGD/fstScan/haplotypeMaps/17WGDcand_names.txt",h=T)
for (i in 1:nrow(namesAll)) { #  i=1
  i2<-namesAll[i,1]
  name<-namesAll[i,2]
aa<-fread(paste("/home/aa/lyrataWGD/fstScan/haplotypeMaps/",name,".HIS.union.txt",sep=""),h=T)
aa$int<-paste(aa$scaff,":",aa$start,sep="")
write.table(file = "/home/aa/lyrataWGD/fstScan/haplotypeMaps/HIS0.7Syn17May.intervals",aa$int,row.names = F,col.names = F,quote = F,append = T)
}
#1. extract using filterHIS script (in /auto/pruhonice1-ibot/home/holcovam/ScanTools/17WGDcandidatesAdaptiveSNPs/)

  
  
  
  
############ ENDOREDUPLICATION ###########
  library(ggpubr)
  setwd("/home/aa/lyrataWGD/endoredupl/")
  df<-read.table("arabidopsis_lyrata_arenosa_endoreduplication_MajdaModif.csv",h=T,sep = ",")
  
  df$Eprop<-rowSums(df[,7:12])/rowSums(df[,6:12])
  df$ENum<-rowSums(df[,7:12] > 4)
  df$EInt<-df$X4C_COUNT*4+df$X8C_COUNT*8+df$X16_COUNT*16+df$X32C_COUNT*32+df$X64_COUNT*64+df$X128_COUNT*128
  
  
  mean(subset(df,df$PLOIDY %in% "2x")$Eprop)
  mean(subset(df,df$PLOIDY %in% "4x")$Eprop)
  mean(subset(df,df$PLOIDY %in% "neo-4x")$Eprop)
  
  wilcox.test(x = subset(df,df$PLOIDY %in% "2x")$Eprop,y = subset(df,df$PLOIDY %in% "4x")$Eprop)
  wilcox.test(x = subset(df,df$PLOIDY %in% "neo-4x")$Eprop,y = subset(df,df$PLOIDY %in% "4x")$Eprop)
  

  median(subset(df,df$PLOIDY %in% "2x")$ENum)
  median(subset(df,df$PLOIDY %in% "4x")$ENum)
  median(subset(df,df$PLOIDY %in% "neo-4x")$ENum)
  
  wilcox.test(x = subset(df,df$PLOIDY %in% "2x")$ENum,y = subset(df,df$PLOIDY %in% "4x")$ENum)
  wilcox.test(x = subset(df,df$PLOIDY %in% "2x")$ENum,y = subset(df,df$PLOIDY %in% "4x")$ENum)
  
  
  
  my_comparisons <- list( c("2x", "4x") )
  pdf("endorEP.pdf",width = 4,height = 4,pointsize = 24)
  #EP
  ggviolin(df, x = "PLOIDY", y = "Eprop", fill = "PLOIDY",
           palette = c("red", "blue"),
           add = "boxplot", add.params = list(fill = "white")) +
    stat_compare_means(comparisons = my_comparisons, label = "p.signif") + 
    stat_compare_means(label.y = 1.5)
  #EI
  ggviolin(df, x = "PLOIDY", y = "EInt", fill = "PLOIDY",
           palette = c("red", "blue"),
           add = "boxplot", add.params = list(fill = "white")) +
    stat_compare_means(comparisons = my_comparisons, label = "p.signif") + 
    stat_compare_means(label.y = 1.5)
  #EN
  ggviolin(df, x = "PLOIDY", y = "ENum", fill = "PLOIDY",
           palette = c("red", "blue"),
           add = "boxplot", add.params = list(fill = "white")) +
    stat_compare_means(comparisons = my_comparisons, label = "p.signif") + 
    stat_compare_means(label.y = 1.5)
  dev.off()
  
  
  df$LP<-paste(df$Lineage,df$PLOIDY,sep=".")
  df<-df[order(df$LP,decreasing = F),]
  my_comparisons <- list(c("AACeur.2x", "AACeur.4x"),c("AACeur.2x", "AACeur.neo-4x"),c("AACeur.4x", "AACeur.neo-4x"),c("Austria.2x", "Austria.4x"),c("Czechia.2x", "Czechia.4x"), c("Germany.2x", "Germany.4x"))
  pdf("endor_All4_neo.pdf",width = 6,height = 5,pointsize = 24)
  #EP
  ggviolin(df, x = "LP", y = "Eprop", fill = "LP",
           palette = c(c("red", "blue","blue"),rep(c("red", "blue"),3)),
           add = "boxplot", add.params = list(fill = "white")) + stat_compare_means(comparisons = my_comparisons, label = "p.signif") +  stat_compare_means(label.y = 1.5)
  ggviolin(df, x = "LP", y = "EInt", fill = "LP",
           palette = c(c("red", "blue","blue"),rep(c("red", "blue"),3)),
           add = "boxplot", add.params = list(fill = "white")) + stat_compare_means(comparisons = my_comparisons, label = "p.signif") +  stat_compare_means(label.y = 1.5)
  ggviolin(df, x = "LP", y = "ENum", fill = "LP",
           palette = c(c("red", "blue","blue"),rep(c("red", "blue"),3)),
           add = "boxplot", add.params = list(fill = "white"))  + stat_compare_means(comparisons = my_comparisons, label = "p.signif") +  stat_compare_means(label.y = 1.5)
  dev.off()

  ##### Neo in arenosa
  library(ggpubr)
  library(data.table)
  setwd("/home/aa/lyrataWGD/endoredupl/")
  df<-read.table("arabidopsis_lyrata_arenosa_endoreduplication_MajdaModif.csv",h=T,sep = ",")
  aa<-subset(df,df$Lineage %in% "AACeur")
  dip<-subset(aa,aa$PLOIDY %in% "2x")
  tet<-subset(aa,aa$PLOIDY %like% "4x")
  dip$dna<-dip$X2C_COUNT*2+dip$X4C_COUNT*4+dip$X8C_COUNT*8+dip$X16_COUNT*16+dip$X32C_COUNT*32+dip$X64_COUNT*64+dip$X128_COUNT*128
  tet$dna<-tet$X2C_COUNT*4+tet$X4C_COUNT*8+tet$X8C_COUNT*16+tet$X16_COUNT*32+tet$X32C_COUNT*64+tet$X64_COUNT*128+tet$X128_COUNT*256
  df1<-rbind(dip,tet)
  df1$perCell<-df1$dna/(df1$X2C_COUNT+df1$X4C_COUNT+df1$X8C_COUNT+df1$X16_COUNT+df1$X32C_COUNT+df1$X64_COUNT+df1$X128_COUNT)
  
  
  df1$Eprop<-rowSums(df1[,7:12])/rowSums(df1[,6:12])
  df1$ENum<-rowSums(df1[,7:12] > 4)
  
  mean(subset(df1,df1$PLOIDY %in% "2x")$dna)
  mean(subset(df1,df1$PLOIDY %in% "4x")$dna)
  mean(subset(df1,df1$PLOIDY %in% "neo-4x")$dna)
  mean(subset(df1,df1$PLOIDY %in% "2x")$perCell)
  mean(subset(df1,df1$PLOIDY %in% "4x")$perCell)
  mean(subset(df1,df1$PLOIDY %in% "neo-4x")$perCell)
  mean(subset(df1,df1$PLOIDY %in% "2x")$Eprop)
  mean(subset(df1,df1$PLOIDY %in% "4x")$Eprop)
  mean(subset(df1,df1$PLOIDY %in% "neo-4x")$Eprop)
  median(subset(df1,df1$PLOIDY %in% "2x")$ENum)
  median(subset(df1,df1$PLOIDY %in% "4x")$ENum)
  median(subset(df1,df1$PLOIDY %in% "neo-4x")$ENum)
  
  my_comparisons <- list( c("2x", "neo-4x", "4x") )
  pdf("endorDNA.pdf",width = 4,height = 4,pointsize = 24)
  #EP
  ggviolin(df1, x = "PLOIDY", y = "perCell", fill = "PLOIDY",
           palette = c("red", "#8800ffff","blue"),
           add = "boxplot", add.params = list(fill = "white")) +
    stat_compare_means(comparisons = my_comparisons, label = "p.signif") + 
    stat_compare_means(label.y = 1.5)
  dev.off()
  
  
  df$LP<-paste(df$Lineage,df$PLOIDY,sep=".")
  df<-df[order(df$LP,decreasing = F),]
  my_comparisons <- list(c("AACeur.2x", "AACeur.4x"),c("AACeur.2x", "AACeur.neo-4x"),c("AACeur.4x", "AACeur.neo-4x"),c("Austria.2x", "Austria.4x"),c("Czechia.2x", "Czechia.4x"), c("Germany.2x", "Germany.4x"))
  pdf("endor_All4_neo.pdf",width = 6,height = 5,pointsize = 24)
  #EP
  ggviolin(df, x = "LP", y = "Eprop", fill = "LP",
           palette = c(c("red", "blue","blue"),rep(c("red", "blue"),3)),
           add = "boxplot", add.params = list(fill = "white")) + stat_compare_means(comparisons = my_comparisons, label = "p.signif") +  stat_compare_means(label.y = 1.5)
  ggviolin(df, x = "LP", y = "EInt", fill = "LP",
           palette = c(c("red", "blue","blue"),rep(c("red", "blue"),3)),
           add = "boxplot", add.params = list(fill = "white")) + stat_compare_means(comparisons = my_comparisons, label = "p.signif") +  stat_compare_means(label.y = 1.5)
  ggviolin(df, x = "LP", y = "ENum", fill = "LP",
           palette = c(c("red", "blue","blue"),rep(c("red", "blue"),3)),
           add = "boxplot", add.params = list(fill = "white"))  + stat_compare_means(comparisons = my_comparisons, label = "p.signif") +  stat_compare_means(label.y = 1.5)
  dev.off()
  
  
  
  
  
  
  
  
  #### DIPLOID
  df<-subset(df,df$PLOIDY %in% "2x")
  df$LP<-paste(df$Lineage,df$PLOIDY,sep=".")
  df<-df[order(df$LP,decreasing = F),]
  my_comparisons <- list(c("AACeur.2x", "Austria.2x"),c("Austria.2x", "Czechia.2x"),c("Czechia.2x", "AACeur.2x"), c("Germany.2x", "Austria.2x"), c("Germany.2x", "Czechia.2x"), c("Germany.2x", "AACeur.2x"))
  pdf("endor_Alldiploids.pdf",width = 6,height = 5,pointsize = 24)
  #EP
  ggviolin(df, x = "LP", y = "Eprop", fill = "LP",
           palette = rep(c("red", "blue"),4),
           add = "boxplot", add.params = list(fill = "white")) + stat_compare_means(comparisons = my_comparisons, label = "p.signif") +  stat_compare_means(label.y = 1.5)
  ggviolin(df, x = "LP", y = "EInt", fill = "LP",
           palette = rep(c("red", "blue"),4),
           add = "boxplot", add.params = list(fill = "white")) + stat_compare_means(comparisons = my_comparisons, label = "p.signif") +  stat_compare_means(label.y = 1.5)
  ggviolin(df, x = "LP", y = "ENum", fill = "LP",
           palette = rep(c("red", "blue"),4),
           add = "boxplot", add.params = list(fill = "white"))  + stat_compare_means(comparisons = my_comparisons, label = "p.signif") +  stat_compare_means(label.y = 1.5)
  dev.off()
  
  hist(df$ENum,nclass = 100)
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  




pdf("mapsHEI10.pdf",width = 14,height = 14)
mapPoints.Ceur <- ggmap(map.Ceur) + geom_scatterpie(aes(x=lon, y=lat,r=rep(0.2,276)), data=haplCeur, cols=c("diploid","tetraploid","other"), alpha=.8) + guides(fill=guide_legend(title="Haplotype"))
mapPoints.Ceur
#+ geom_point(aes(x = lon, y = lat, colour = factor(col)), data = haplCeur)
#+ geom_scatterpie_legend(haplCeur$pop, x=-16, y=-5)
dev.off()

mapPoints.Ceur <- ggmap(map.Ceur) + geom_point(aes(x = lon, y = lat, colour = factor(col)), data = haplCeur) + scale_colour_manual(values = c("red","blue"))


library(ggplot2)
library(scatterpie)
d <- data.frame(x=rnorm(5), y=rnorm(5))
d$A <- abs(rnorm(5, sd=1))
d$B <- abs(rnorm(5, sd=2))
d$C <- abs(rnorm(5, sd=3))
pdf("mapsHEI10.pdf")
ggplot() + geom_scatterpie(aes(x=x, y=y), data=d, cols=c("A", "B", "C")) + coord_fixed()
dev.off()





hapl<-read.table(list,h=T)
haplCeur<-subset(hapl,hapl$lon < 28 & hapl$lon > 2 & hapl$lat < 56)
bbox.Ceur<-make_bbox(lon, lat, data=haplCeur)
map.Ceur <- get_map(location = bbox.Ceur, zoom = 5, source = "stamen", maptype = "terrain")
mapPoints.Ceur <- ggmap(map.Ceur) + geom_point(aes(x = lon, y = lat, colour = factor(col)), data = haplCeur) + scale_colour_manual(values = c("red","blue"))

pdf("mapsHEI10.pdf")
mapPoints.eur
mapPoints.Ceur
dev.off()


world <- map_data('world')
p <- ggplot(world, aes(long, lat)) +
  geom_map(map=world, aes(map_id=region), fill=NA, color="black") +
  coord_quickmap()
p + geom_scatterpie(aes(x=long, y=lat, group=region, r=radius),
                    data=d, cols=LETTERS[1:4], color=NA, alpha=.8) +
  geom_scatterpie_legend(d$radius, x=-160, y=-55)


bbox.eur<-make_bbox(lon, lat, data=hapleur)
map.eur <- get_map(location = bbox.eur, zoom = 5, source = "stamen", maptype = "terrain")
mapPoints.eur <- ggmap(map.eur) + geom_point(aes(x = lon, y = lat, colour = factor(col)), data = hapleur) + scale_colour_manual(values = c("red","blue"))

library(scatterpie)
hapl<-read.table(list,h=T)
haplCeur<-subset(hapl,hapl$lon < 28 & hapl$lon > 2 & hapl$lat < 56)
haplCeur1<-na.omit(haplCeur)
bbox.Ceur<-make_bbox(lon, lat, data=haplCeur)
map.Ceur <- get_map(location = bbox.Ceur, zoom = 5, source = "stamen", maptype = "terrain")























#3. calculate AF for each lineage TODO 
popnames<-unique(substr(names,1,3))
for (pop in popnames) { #  pop="BAB"
  afh1 <- ar3 %>% dplyr:: select(starts_with(pop))
  ar3$ACh1<-rowSums(afh1,na.rm = T)
  ar3$NAh1<-apply(is.na(afh1), 1, sum)
  if (subset(nnn,nnn$V1%like%pop)[1,3]==0) 
  {ar3$ANh1<-(ncol(afh1)-ar3$NAh1)*2} else {ar3$ANh1<-(ncol(afh1)-ar3$NAh1)*4}
  ar3$POP<-ar3$ACh1/ar3$ANh1
  colnames(ar3)[which(colnames(ar3) %in% "POP")]<-pop
}
s4<-dplyr::select(ar3,ID,ann,aas,ZEP, BAB, SUB, HOC, PER, BRE, VLA, CHO, SWA, JOH, MOD,TEM,PEK,SCT,VLH,OSL,STD,tot)
### MORE POSSIBILITIES ###
#s<-subset(s4,!s4$ann %in% "intragenic_variant" & !s4$ann %in% "downstream_gene_variant" & !s4$ann %in% "upstream_gene_variant")
#s<-subset(s4,!s4$ann %in% "intragenic_variant")
s<-subset(s4,s4$ann %like% "missense_variant")
#repolarise

#for (i in  1:nrow(s)){ # i=1
#  if (s$tot[i]>0.5)
#  {s[i,4:(length(s[i,])-1)]<-1-s[i,4:(length(s[i,])-1)]
#  } else {}}

#repolarise
popnames<-colnames(s4)[4:(length(colnames(s4))-1)]
for (i in  1:nrow(s)){ # i=1
  d<-c()
  t<-c()
  for (pop in popnames) { #  pop="BAB"
    ppp <- s %>% dplyr:: select(starts_with(pop))
    if (subset(nnn,nnn$V1%like%pop)[1,3]==0) 
    {d<-c(d,as.numeric(as.character((ppp[i,1])[[1]])))} else {t<-c(t,as.numeric(as.character((ppp[i,1])[[1]])))}
  }
  meandip<-mean(d,na.rm = T)
  meantet<-mean(t,na.rm = T)
  if (meandip>meantet)
  {s[i,4:(length(s[i,])-1)]<-(1-s[i,4:(length(s[i,])-1)])
  print(i)
  } else {}}



#plot
ann<-fread("/home/aa/Desktop/references/lyrataV2/LyV2_TAIR11orth_des_20171231.txt")
pops<-c(colnames(s)[4:(length(colnames(s))-1)])
my_palette <- colorRampPalette(c("khaki1", "green2", "blue3"))(n = 100)
for (id in i2){ #   id = "AL1G48330"
  if (nrow(subset(s,s$ID %in% id))>1)
  {s1<-subset(s,s$ID %in% id)
  ann1<-subset(ann,ann$AL %in% id)
  if (nrow(s1)>110)
  {p=0.75
  } else {p=1}
  s5<-dplyr::select(s1,SUB, HOC, PER, JOH, MOD, VLH, BAB, BRE, VLA, TEM, PEK, OSL, ZEP, CHO, SWA, SCT, STD)
  df<-as.matrix(s1[,4:(length(s1[i,])-1)],rownames = s1$aas)
  df1<-as.matrix(s5,rownames = s1$aas)
  
  bbox.1<-make_bbox(lon, lat, data=coords)
  map <- get_map(location = bbox.1 , zoom = 5, source = "stamen", maptype = "terrain")
  mapPoints <- ggmap(map) + geom_point(aes(x = coords$lon, y = coords$lat, colour = factor(grp.pop$value)),
                                       data = grp.pop) + scale_colour_manual(values = c("red","blue","black","green"))
  mapPoints
  
  } else {} 
}











#4. N associated GOs
# I'm using gene list with manual selection but the simple .ann.txt is also fine
library(data.table)
library(ggpubr)
setwd("/home/aa/2alpine/predictors/N_GOs")
go<-fread("file:///home/aa/Desktop/references/lyrataV2/functions/ATH_GO_GOSLIMOct2020.txt",h=T)
cand<-fread("file:///home/aa/2alpine/fstScan/outGenes_SNP0.99_Gene0.1.ALAT.ann.123.txt",h=T)
go$ATGO<-paste(go$AT,go$GOcode,sep="_")
go1<-go[ order(go[,16]), ]
go2<-go1[!duplicated(go1[,c('ATGO')]),] 
# All GOs, all genes background
genNgos<-as.data.frame(table(go2$AT))
genNgosSel<-subset(genNgos,genNgos$Var1 %in% substr(cand$AT,1,9))
genNgosNonSel<-subset(genNgos,!genNgos$Var1 %in% substr(cand$AT,1,9))
hist(genNgosSel$Freq,breaks = 50)
hist(genNgosNonSel$Freq,breaks = 50)
summary(genNgosSel)
summary(genNgosNonSel)
###visualization
sel <- cbind("Parallel",genNgosSel)
non <- cbind("Any other",genNgosNonSel)
colnames(non)<-c('"Parallel"',"Var1","Freq")
all<-as.data.frame(rbind(sel,non))
colnames(all)<-c('Gene',"Var1","Freq")
all$Freq<-as.numeric(as.character(all$Freq))
my_comparisons <- list( c("Parallel", "Any other") )
pdf("GOs_all_Allgenes.pdf",width = 4,height = 4,pointsize = 24)
ggviolin(all, x = 'Gene', y = "Freq", fill = "Gene",
         palette = c("red", "grey50"),
         add = "boxplot", add.params = list(fill = "white")) +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif") + 
  stat_compare_means(label.y = max(all$Freq)-5)
dev.off()
#all GOs, all selected background
dict<-fread("/home/aa/Desktop/references/lyrataV2/functions/ALATdict.txt")
dat<-read.table("file:///home/aa/2alpine/fstScan/outGenes_SNP0.99_Gene0.1.txt",sep="\t")
datc<-""
for (i in 1:length(dat$V1)) {
 datc<-paste(datc,unlist(dat$V1[i]), sep='', collapse=NULL)}
datc1<-as.data.frame(strsplit(datc, split=" ")[[1]][-1])
datc1$AllSelected<-substr(datc1$`strsplit(datc, split = " ")[[1]][-1]`,1,9)
setDT(datc1)
datc1[, Count := .N, by=AllSelected]
datc2<-datc1[Count==1]
datc3<-subset(dict,dict$AL %in% datc2$AllSelected)
genNgosNonSel<-subset(genNgos,genNgos$Var1 %in% substr(datc3$AT,1,9))
hist(genNgosSel$Freq,breaks = 50)
hist(genNgosNonSel$Freq,breaks = 50)
summary(genNgosSel)
summary(genNgosNonSel)
sel <- cbind("Parallel",genNgosSel)
non <- cbind("Non-parallel selected",genNgosNonSel)
colnames(non)<-c('"Parallel"',"Var1","Freq")
all<-as.data.frame(rbind(sel,non))
colnames(all)<-c('Gene',"Var1","Freq")
all$Freq<-as.numeric(as.character(all$Freq))
my_comparisons <- list( c("Parallel", "Non-parallel selected") )
pdf("GOs_all_Allselected.pdf",width = 4,height = 4,pointsize = 24)
ggviolin(all, x = 'Gene', y = "Freq", fill = "Gene",
         palette = c("red", "grey50"),
         add = "boxplot", add.params = list(fill = "white")) +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif") + 
  stat_compare_means(label.y = max(all$Freq)-5)
dev.off()

### By category - the loop does not output PDF :(
cat="C" # "F" "P"
print(cat)
print("All")
go3<-subset(go2,go2$category %in% cat)
#all genes background
genNgos<-as.data.frame(table(go3$AT))
genNgosSel<-subset(genNgos,genNgos$Var1 %in% substr(cand$AT,1,9))
genNgosNonSel<-subset(genNgos,!genNgos$Var1 %in% substr(cand$AT,1,9))
#hist(genNgosSel$Freq,breaks = 50)
#hist(genNgosNonSel$Freq,breaks = 50)
print(summary(genNgosSel))
print(summary(genNgosNonSel))
sel <- cbind("Parallel",genNgosSel)
non <- cbind("Any other",genNgosNonSel)
colnames(non)<-c('"Parallel"',"Var1","Freq")
all<-as.data.frame(rbind(sel,non))
colnames(all)<-c('Gene',"Var1","Freq")
all$Freq<-as.numeric(as.character(all$Freq))
my_comparisons <- list( c("Parallel", "Any other") )
pdf(paste("GOs_",cat,"_Allgenes.pdf",sep=""),width = 4,height = 4,pointsize = 24)
ggviolin(all, x = 'Gene', y = "Freq", fill = "Gene",
         palette = c("red", "grey50"),
         add = "boxplot", add.params = list(fill = "white")) +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif") + 
  stat_compare_means(label.y = max(all$Freq)-0.5)
dev.off()
#all selected background
print("selected")
dict<-fread("/home/aa/Desktop/references/lyrataV2/functions/ALATdict.txt")
dat<-read.table("file:///home/aa/2alpine/fstScan/outGenes_SNP0.99_Gene0.1.txt",sep="\t")
datc<-""
for (i in 1:length(dat$V1)) {
  datc<-paste(datc,unlist(dat$V1[i]), sep='', collapse=NULL)}
datc1<-as.data.frame(strsplit(datc, split=" ")[[1]][-1])
datc1$AllSelected<-substr(datc1$`strsplit(datc, split = " ")[[1]][-1]`,1,9)
setDT(datc1)
datc1[, Count := .N, by=AllSelected]
datc2<-datc1[Count==1]
datc3<-subset(dict,dict$AL %in% datc2$AllSelected)
genNgosNonSel<-subset(genNgos,genNgos$Var1 %in% substr(datc3$AT,1,9))
#hist(genNgosSel$Freq,breaks = 50)
#hist(genNgosNonSel$Freq,breaks = 50)
print(summary(genNgosSel))
print(summary(genNgosNonSel))
sel <- cbind("Parallel",genNgosSel)
non <- cbind("Non-parallel selected",genNgosNonSel)
colnames(non)<-c('"Parallel"',"Var1","Freq")
all<-as.data.frame(rbind(sel,non))
colnames(all)<-c('Gene',"Var1","Freq")
all$Freq<-as.numeric(as.character(all$Freq))
my_comparisons <- list( c("Parallel", "Non-parallel selected") )
pdf(paste("GOs_",cat,"_Allselected.pdf",sep=""),width = 4,height = 4,pointsize = 24)
ggviolin(all, x = 'Gene', y = "Freq", fill = "Gene",
         palette = c("red", "grey50"),
         add = "boxplot", add.params = list(fill = "white")) +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif") + 
  stat_compare_means(label.y = max(all$Freq)-0.5)
dev.off()
  
  
###5. N interactors STRING
library(data.table)
setwd("/home/aa/2alpine/predictors/STRING/")
s1<-read.table("3702.protein.links.v11.0.txt",h=T,comment.char = "")

score = 500 # (400-medium confidnce, 700 - high conficence)
s<-subset(s1,s1$combined_score >= score)
s$x<-substr(s$protein1 ,6,14)
nInt<-as.data.frame(table(s$x))
hist(nInt$Freq,breaks =100)
summary(nInt)

### Parallel
cand<-fread("file:///home/aa/2alpine/fstScan/outGenes_SNP0.99_Gene0.1.ALAT.ann.123.txt",h=T)
cand$AT<-substr(cand$AT,1,9)
candInt<-subset(nInt,nInt$Var1 %in% cand$AT)
hist(candInt$Freq,breaks =100)
summary(candInt)

### Selected
#To make a list of all selected genes, without parallel ones
dict<-fread("/home/aa/Desktop/references/lyrataV2/functions/ALATdict.txt")
dat<-read.table("file:///home/aa/2alpine/fstScan/outGenes_SNP0.99_Gene0.1.txt",sep="\t")
datc<-""
for (i in 1:length(dat$V1)) {
  datc<-paste(datc,unlist(dat$V1[i]), sep='', collapse=NULL)}
datc1<-as.data.frame(strsplit(datc, split=" ")[[1]][-1])
datc1$AllSelected<-substr(datc1$`strsplit(datc, split = " ")[[1]][-1]`,1,9)
setDT(datc1)
datc1[, Count := .N, by=AllSelected]
datc2<-datc1[Count==1]
datc3<-subset(dict,dict$AL %in% datc2$AllSelected)
datc3<-subset(datc3,!datc3$AT %in% "nnn")
datc3$AT<-substr(datc3$AT,1,9)
selInt<-subset(nInt,nInt$Var1 %in% datc3$AT)
hist(selInt$Freq,breaks =100)
summary(selInt)

#parallel-all
print(c(mean(candInt$Freq),mean(nInt$Freq)))
print(c(median(candInt$Freq),median(nInt$Freq)))
print(wilcox.test(candInt$Freq,nInt$Freq))
#parallel-selected
print(c(mean(candInt$Freq),mean(selInt$Freq)))
print(c(median(candInt$Freq),median(selInt$Freq)))
print(wilcox.test(candInt$Freq,selInt$Freq))

### Among low inetractors
freq = 150
#parallel-all
print(c(mean(subset(candInt, candInt$Freq <= freq)$Freq),mean(subset(nInt, nInt$Freq <= freq)$Freq)))
print(wilcox.test(subset(candInt, candInt$Freq <= freq)$Freq,subset(nInt, nInt$Freq <= freq)$Freq))
#Plot
pdf("stringInteractorsRatio_conf200.pdf")
par(mfrow=c(2,1))
#MEAN
mat<-matrix(nrow = 160,ncol = 4,dimnames = list(c(),c("max_N_int","N_int_parallel","N_int_any","N_int_parallel/N_int_any")))
for (i in 1:160) { # i=1
  seq<-seq(10,1600,10)[i]
  mat[i,1]<-seq
  mat[i,2]<-mean(subset(candInt, candInt$Freq <= seq)$Freq)
  mat[i,3]<-mean(subset(nInt, nInt$Freq <= seq)$Freq)
  mat[i,4]<-mean(subset(candInt, candInt$Freq <= seq)$Freq)/mean(subset(nInt, nInt$Freq <= seq)$Freq)}
mat1<-as.data.frame(mat)
plot(mat1$`N_int_parallel/N_int_any`~mat1$max_N_int)
abline(h=1)

#MEDIAN
mat<-matrix(nrow = 160,ncol = 4,dimnames = list(c(),c("max_N_int","N_int_parallel","N_int_any","N_int_parallel/N_int_any")))
for (i in 1:160) { # i=1
  seq<-seq(10,1600,10)[i]
  mat[i,1]<-seq
  mat[i,2]<-median(subset(candInt, candInt$Freq <= seq)$Freq)
  mat[i,3]<-median(subset(nInt, nInt$Freq <= seq)$Freq)
  mat[i,4]<-median(subset(candInt, candInt$Freq <= seq)$Freq)/median(subset(nInt, nInt$Freq <= seq)$Freq)}
mat1<-as.data.frame(mat)
plot(mat1$`N_int_parallel/N_int_any`~mat1$max_N_int)
abline(h=1)
dev.off()
#what are the high frequency?
df<-subset(candInt, candInt$Freq > 600)
df2<-subset(cand,cand$AT %in% df$Var1)

###6. N interactors THALIANA NETWORK
library(data.table)
setwd('/home/aa/2alpine/predictors/thalianaNetwork/')
load("WGCNA_network_6Cand16C.rdata")
con<-WGCNA_network_6Cand16C[[5]]
hist(con$kTotal,breaks =100)
summary(con)
### Parallel
cand<-fread("file:///home/aa/2alpine/fstScan/outGenes_SNP0.99_Gene0.1.ALAT.ann.123.txt",h=T)
cand$AT<-substr(cand$AT,1,9)
candInt<-subset(con,rownames(con) %in% cand$AT)
hist(candInt$kTotal,breaks =100)
summary(candInt)
#parallel-all
print(c(mean(candInt$kTotal),mean(con$kTotal)))
print(c(median(candInt$kTotal),median(con$kTotal)))
print(wilcox.test(candInt$kTotal,con$kTotal))
### Among low inetractors
pdf("thalianaNetworkRatio.pdf")
par(mfrow=c(2,1))
#MEAN
mat<-matrix(nrow = 60,ncol = 4,dimnames = list(c(),c("max_N_int","N_int_parallel","N_int_any","N_int_parallel/N_int_any")))
for (i in 1:60) { # i=1
  seq<-seq(5,300,5)[i]
  mat[i,1]<-seq
  mat[i,2]<-mean(subset(candInt, candInt$kTotal <= seq)$kTotal)
  mat[i,3]<-mean(subset(con, con$kTotal <= seq)$kTotal)
  mat[i,4]<-mean(subset(candInt, candInt$kTotal <= seq)$kTotal)/mean(subset(con, con$kTotal <= seq)$kTotal)}
mat1<-as.data.frame(mat)
plot(mat1$`N_int_parallel/N_int_any`~mat1$max_N_int)
abline(h=1)
#MEDIAN
mat<-matrix(nrow = 60,ncol = 4,dimnames = list(c(),c("max_N_int","N_int_parallel","N_int_any","N_int_parallel/N_int_any")))
for (i in 1:60) { # i=1
  seq<-seq(5,300,5)[i]
  mat[i,1]<-seq
  mat[i,2]<-median(subset(candInt, candInt$kTotal <= seq)$kTotal)
  mat[i,3]<-median(subset(con, con$kTotal <= seq)$kTotal)
  mat[i,4]<-median(subset(candInt, candInt$kTotal <= seq)$kTotal)/median(subset(con, con$kTotal <= seq)$kTotal)}
mat1<-as.data.frame(mat)
plot(mat1$`N_int_parallel/N_int_any`~mat1$max_N_int)
abline(h=1)
dev.off()

###7. The size of the gene
library(data.table)
library(dplyr)
setwd('/home/aa/2alpine/predictors/geneSize/')
cand<-fread("file:///home/aa/2alpine/fstScan/outGenes_SNP0.99_Gene0.1.ALAT.ann.123.txt",h=T)
g<-read.table("file:///home/aa/Desktop/references/lyrataV2/genesIDsLyV2.gff",h=T)
g$ID<-substr(g$ID,4,12)
gs<-subset(g,g$ID %in% cand$AL)
g$dist<-g$end-g$start
gs$dist<-gs$end-gs$start
q<-quantile(g$dist[g$dist>100],.999)

pdf("lengthTotal.pdf",width = 16,height =16,pointsize = 24)
par(mfrow=c(2,1))
p1 <- hist(g$dist[g$dist>100],breaks = 200,main = "Genome-wide",xlab = "Gene lenght",freq = T,xlim = c(0,q))
abline(v=median(g$dist[g$dist>100]), col="blue")
p2 <- hist(gs$dist[gs$dist>100],breaks = 50,main = "Parallel",xlab = "Gene lenght",col = "red",freq = T,xlim = c(0,q)) 
abline(v=median(gs$dist[gs$dist>100]), col="blue")
dev.off()

g1<-sample_n(tbl = g,size = nrow(gs))
print("genes")
print(nrow(gs))
print("subset")
print(wilcox.test(x = g1$dist,y = gs$dist))
print("genome")
print(wilcox.test(x = g$dist,y = gs$dist))
print("median_genes")
print(median(gs$dist,na.rm = T))
print("mean_genes")
print(mean(gs$dist,na.rm = T))
median(g$dist[g$dist>100],na.rm = T)
mean(g$dist[g$dist>10],na.rm = T)

### Length upstream region
g$upstream<-""
for (i in 2:nrow(g)) { # i=2
  if (g[i,7]=="+") {
    ups<-g[i,4]-g[i-1,5]
  } else {ups<-g[i+1,4]-g[i,5]}
  g$upstream[i]<-as.numeric(ups)
}

g$upstream<-as.numeric(g$upstream)
g$upstream[g$upstream<0] <- 0
g$upstream[g$upstream>quantile(g$upstream,.99,na.rm = T)] <- NA
hist(as.numeric(g$upstream[2:nrow(g)]),breaks = 200)
summary(as.numeric(g$upstream))
gs<-subset(g,g$ID %in% cand$AL)

pdf("lengthUpstream.pdf",width = 16,height =16,pointsize = 24)
par(mfrow=c(2,1))
p1 <- hist(g$upstream,breaks = 200,main = "Genome-wide",xlab = "Upstream region length",freq = T,xlim = c(0,q))
abline(v=median(g$upstream), col="blue")
p2 <- hist(gs$dist[gs$dist>100],breaks = 50,main = "Parallel",xlab = "Upstream region length",col = "red",freq = T,xlim = c(0,q)) 
abline(v=median(gs$dist[gs$dist>100]), col="blue")
dev.off()

g1<-sample_n(tbl = g,size = nrow(gs))
print("subset")
print(wilcox.test(x = g1$upstream,y = gs$upstream))
print("genome")
print(wilcox.test(x = g$upstream,y = gs$upstream))
print("median_genes")
print(median(gs$upstream,na.rm = T))
print("mean_genes")
print(mean(gs$upstream,na.rm = T))
median(g$upstream,na.rm = T)
mean(g$upstream,na.rm = T)

### 8. Araport11 data
# tissue specific
library(data.table)
setwd("/home/aa/2alpine/predictors/araportExpression/")
ara<-read.table("Araport11_TissueSpecific.txt")
cand<-fread("file:///home/aa/2alpine/fstScan/outGenes_SNP0.99_Gene0.1.ALAT.ann.123.txt",h=T)
cand$AT<-substr(cand$AT,1,9)
candInt<-subset(ara,ara$V1 %in% cand$AT)
fisher.test(matrix(c(nrow(candInt), nrow(cand)-nrow(candInt), nrow(ara)-nrow(candInt),34052-nrow(ara)-nrow(candInt)-nrow(cand)), nrow=2))
#Proportion of tissue specific in parallel
nrow(candInt)/nrow(cand)
#Proportion of tissue specific genome-wide
(nrow(ara)-nrow(candInt))/(34052-nrow(ara)-nrow(candInt)-nrow(cand))
##Output
CandInt2<-subset(cand,cand$AT %in% ara$V1)
tot<-cbind(candInt,CandInt2)
write.table(tot,"TissueSpecificParallel.txt",quote = F,row.names = F,sep = "\t")

# housekeeping
library(data.table)
setwd("/home/aa/2alpine/predictors/araportExpression/")
ara<-read.table(file = "Araport11_HouseKeeping.txt")
ara$V1<-substr(ara$V1,1,9)
cand<-fread("file:///home/aa/2alpine/fstScan/outGenes_SNP0.99_Gene0.1.ALAT.ann.123.txt",h=T)
cand$AT<-substr(cand$AT,1,9)
candInt<-subset(ara,ara$V1 %in% cand$AT)
fisher.test(matrix(c(nrow(candInt), nrow(cand)-nrow(candInt), nrow(ara)-nrow(candInt),34052-nrow(ara)-nrow(candInt)-nrow(cand)), nrow=2))
#Proportion of house-keeping in parallel
nrow(candInt)/nrow(cand)
#Proportion of house-keeping genome-wide
(nrow(ara)-nrow(candInt))/(34052-nrow(ara)-nrow(candInt)-nrow(cand))
CandInt2<-subset(cand,cand$AT %in% ara$V1)
tot<-cbind(candInt,CandInt2)
write.table(tot,"HouseKeepingParallel.txt",quote = F,row.names = F,sep = "\t")

# in all tissues
library(data.table)
setwd("/home/aa/2alpine/predictors/araportExpression/")
ara<-read.table(file = "Araport11_inAllTissues.txt")
ara$V1<-substr(ara$V1,1,9)
cand<-fread("file:///home/aa/2alpine/fstScan/outGenes_SNP0.99_Gene0.1.ALAT.ann.123.txt",h=T)
cand$AT<-substr(cand$AT,1,9)
candInt<-subset(ara,ara$V1 %in% cand$AT)
fisher.test(matrix(c(nrow(candInt), nrow(cand)-nrow(candInt), nrow(ara)-nrow(candInt),34052-nrow(ara)-nrow(candInt)-nrow(cand)), nrow=2))
#Proportion in parallel
nrow(candInt)/nrow(cand)
#Proportion genome-wide
(nrow(ara)-nrow(candInt))/(34052-nrow(ara)-nrow(candInt)-nrow(cand))
CandInt2<-subset(cand,cand$AT %in% ara$V1)
tot<-cbind(candInt,CandInt2)
write.table(tot,"inAllTissuesParallel.txt",quote = F,row.names = F,sep = "\t")

## Loci with associated cDNAs (i.e cDNAs which overlap the locus, November 2010)
library(data.table)
setwd("/home/aa/2alpine/predictors/araportExpression/")
ara<-read.table(file = "TAIR10_Locus_cDNA_associations.txt",sep = "\t",h=T)
ara1<-as.data.frame(table(ara$name))
cand<-fread("file:///home/aa/2alpine/fstScan/outGenes_SNP0.99_Gene0.1.ALAT.ann.123.txt",h=T)
cand$AT<-substr(cand$AT,1,9)
candInt<-subset(ara1,ara1$Var1 %in% cand$AT)
#visualize
pdf("associatedcDNAs.pdf",width = 16,height =16,pointsize = 24)
par(mfrow=c(2,1))
p1 <- hist(ara1$Freq,breaks = 100,main = "Genome-wide",xlab = "N of associated cDNAs",freq = T,xlim = c(0,quantile(ara1$Freq,.999)))
abline(v=median(ara1$Freq), col="blue")
p2 <- hist(candInt$Freq,breaks = 50,main = "Parallel",xlab = "N of associated cDNAs",freq = T,xlim = c(0,quantile(ara1$Freq,.999)),col = "red")
abline(v=median(candInt$Freq), col="blue")
dev.off()

print("genome")
print(wilcox.test(x = ara1$Freq,y =candInt$Freq))
print("median_genes")
print(median(candInt$Freq,na.rm = T))
print("mean_genes")
print(mean(candInt$Freq,na.rm = T))
median(ara1$Freq,na.rm = T)
mean(ara1$Freq,na.rm = T)


#TODO - expression level
# N of tissues it is expressed in
# araport11 - rna-seq based evidence - transcript assembly - expression level?

# Diversity at the locus - RANKS
library(data.table)
setwd("/home/aa/2alpine/predictors/diversity/")
pops<-c("LOM","SUN","HIF","HFF","DRG","GUN","SUB")
cand<-fread("file:///home/aa/2alpine/fstScan/outGenes_SNP0.99_Gene0.1.ALAT.ann.123.txt",h=T)
genes<-fread("/home/aa/Desktop/references/lyrataV2/genesIDsLyV2.gff",h=T)
genes$ID<-substr(genes$ID,4,12)
ranks<-""
totranks<-''
for (pop in pops) { # pop = "SUB"
  div<-read.table(paste("/home/aa/2alpine/PopStructure/WPM/",pop,".WS50.0k_MS1_4ind_WPM.txt",sep=""),h=T,nrow = length(readLines(paste("/home/aa/2alpine/PopStructure/WPM/",pop,".WS50.0k_MS1_4ind_WPM.txt",sep=""))) - 2)
  div1<-as.data.frame(div[ order(div[,4],div[,5]), ])
  div2<-div1[!duplicated(div1[,c(4,5)],fromLast = T),]
  div2<-setDT(div2[ order(div2[,13]), ])
     div2<-subset(div2,div2$num_sites>500) ####Could be commented OUT
  div2$rank<-seq(1:nrow(div2))
  cand1<-subset(cand,cand$lineages %like% pop)
  print(pop)
  print(nrow(cand1))
  cand2<-subset(genes,genes$ID %in% cand1$AL)
  setkey(div2, scaff, start, end)
  div3<-foverlaps(cand2, div2, type="any")
  print(nrow(div3))
  ranks<-as.numeric(c(ranks,div3$rank))
  totranks<-as.numeric(c(totranks,rep(seq(1:nrow(div2)),nrow(div3))))} 
wilcox.test(ranks,totranks)
print("median_genes")
print(median(ranks,na.rm = T))
print("mean_genes")
print(mean(ranks,na.rm = T))
median(totranks,na.rm = T)
mean(totranks,na.rm = T)
#visualize
pdf("diversityRanks_moreThan500sites.pdf",width = 16,height =16,pointsize = 24)
par(mfrow=c(2,1))
p1 <- hist(totranks,breaks = 100,main = "Genome-wide",xlab = "Rank - nucleotide diversity",freq = T,xlim = c(0,max(totranks,na.rm = T)))
abline(v=median(totranks,na.rm = T), col="blue")
p2 <- hist(ranks,breaks = 100,main = "Parallel",xlab = "Rank - nucleotide diversity",freq = T,xlim = c(0,max(totranks,na.rm = T)),col = "red")
abline(v=median(ranks,na.rm = T), col="blue")
dev.off()
  
# Diversity at the locus - DIVERSITY VALUES
diver<-""
totdiver<-''
for (pop in pops) { # pop = "SUB"
  div<-read.table(paste("/home/aa/2alpine/PopStructure/WPM/",pop,".WS50.0k_MS1_4ind_WPM.txt",sep=""),h=T,nrow = length(readLines(paste("/home/aa/2alpine/PopStructure/WPM/",pop,".WS50.0k_MS1_4ind_WPM.txt",sep=""))) - 2)
  div1<-as.data.frame(div[ order(div[,4],div[,5]), ])
  div2<-div1[!duplicated(div1[,c(4,5)],fromLast = T),]
  div2<-setDT(div2[ order(div2[,13]), ])
     div2<-subset(div2,div2$num_sites>500) ####Could be commented OUT
  cand1<-subset(cand,cand$lineages %like% pop)
  print(pop)
  print(nrow(cand1))
  cand2<-subset(genes,genes$ID %in% cand1$AL)
  setkey(div2, scaff, start, end)
  div3<-foverlaps(cand2, div2, type="any")
  print(nrow(div3))
  diver<-as.numeric(c(diver,div3$Diversity))
  totdiver<-as.numeric(c(totdiver,rep(div2$Diversity,nrow(div3))))} 
wilcox.test(diver,totdiver)
print("median_genes")
print(median(diver,na.rm = T))
print("mean_genes")
print(mean(diver,na.rm = T))
median(totdiver,na.rm = T)
mean(totdiver,na.rm = T)
#visualize
pdf("diversitydiver_moreThan500sites.pdf",width = 16,height =16,pointsize = 24) # _moreThan500sites
par(mfrow=c(2,1))
p1 <- hist(totdiver,breaks = 100,main = "Genome-wide",xlab = "Rank - nucleotide diversity",freq = T,xlim = c(0,max(totdiver,na.rm = T)))
abline(v=median(totdiver,na.rm = T), col="blue")
p2 <- hist(diver,breaks = 100,main = "Parallel",xlab = "Rank - nucleotide diversity",freq = T,xlim = c(0,max(totdiver,na.rm = T)),col = "red")
abline(v=median(diver,na.rm = T), col="blue")
dev.off()


# Recombination rate at the locus
allDataMatrix<-readRDS("/home/aa/alpine/dmc/genomeScan/BALTISZEPSUB/lyrata_geneDensityData_binned.RDS")
allDataMatrix<-readRDS("/home/aa/alpine/dmc/genomeScan/resources/recMap_fx.RDS")

head(allDataMatrix)
  summary(allDataMatrix)

  
  }

  
### 8. Copy number variation
library(data.table)
setwd("/home/aa/2alpine/predictors/CNV/")
cnv<-fread("output_genotype.vcf",h=T)
colnames(cnv)[1]<-"scaff"
cnv<-subset(cnv,as.numeric(as.character(substr(cnv$scaff,10,15)))<9)

## length
id<-as.data.frame(cnv$ID)
library(stringr)
id2<-as.data.frame(str_split_fixed(id$`cnv$ID`, "-", 2))
colnames(id2)<-c("v1","end")
id2$start<-substr(id2$v1,12,25)
cnv$start<-as.numeric(as.character(id2$start))
cnv$end<-as.numeric(as.character(id2$end))
cnv$length<-abs(cnv$end-cnv$start)
summary(cnv$length)
hist(cnv$length,breaks = 100)
head(cnv)

# extract genotypes and calculate AFs
gt<-cnv[,10:59]
head(gt)
gt1<-as.data.frame(sapply(gt, substring, 1, 3))
gt2 <- as.data.frame(sapply(gt1, function(x) {
gsub("2", "1", x) }))
gt2 <- as.data.frame(sapply(gt2, function(x) {
  gsub("3", "1", x) }))
gt3 <- as.data.frame(lapply(gt2, function(x) {
  sum(as.numeric(as.character(str_split_fixed(x, "/", 2)))) }))
gt3<-as.data.frame(apply(gt2, c(1,2), function(x) {sum(as.numeric(as.character(str_split_fixed(x, "/", 2))))}))
foot<-gt3[,1:28]
alpi<-gt3[,29:50]

cnv$AFfoot<-apply(foot,1,sum)/(2*28)
cnv$AFalp<-apply(alpi,1,sum)/(2*22)
cnv$afd<-abs(cnv$AFfoot-cnv$AFalp)
cnv$AFtot<-apply(gt3,1,sum)/(2*50)

hist(cnv$AFfoot,breaks = 100)
summary(cnv$AFfoot)
hist(cnv$AFalp,breaks = 100)
summary(cnv$AFalp)
hist(cnv$afd,breaks = 100)
summary(cnv$afd)
hist(cnv$AFtot,breaks = 100)

cnvUbiq<-subset(cnv, cnv$AFtot == 1)
write.table(cnvUbiq,"fixedNonReferenceCNVs.txt",col.names = T,row.names = F, quote = F)
cnvOut<-subset(cnv,cnv$afd>quantile(cnv$afd,0.99))
write.table(cnvOut,"Outlier99CNVs.txt",col.names = T,row.names = F, quote = F)
write.table(cnv,"AllCNVs_info.txt",col.names = T,row.names = F, quote = F)


##Fixed non-reference - likely species-wide -> predictor (could still include the variability in N of duplications)
#parallel Arabidopsis-wide candiates
cand<-fread("file:///home/aa/2alpine/fstScan/outGenes_SNP0.99_Gene0.1.ALAT.ann.123.txt",h=T)
ubiq<-fread("file:///home/aa/2alpine/predictors/CNV/fixedNonReferenceCNVs.txt",h=T)
genes<-fread("/home/aa/Desktop/references/lyrataV2/genesIDsLyV2.gff",h=T)
genes$ID<-substr(genes$ID,4,12)
cand2<-subset(genes,genes$ID %in% cand$AL)
setkey(ubiq, scaff, start, end)
ubiq1<-foverlaps(cand2, ubiq, type="any")
# aVT candiates
cand<-fread("file:///home/aa/alpine/arenosaGenome/selScans/bayPass_quartetFst/Tatry2x/genes_full.txt",h=F)
ubiq<-fread("file:///home/aa/2alpine/predictors/CNV/fixedNonReferenceCNVs.txt",h=T)
genes<-fread("/home/aa/Desktop/references/lyrataV2/genesIDsLyV2.gff",h=T)
genes$ID<-substr(genes$ID,4,12)
cand2<-subset(genes,genes$ID %in% cand$V1)
setkey(ubiq, scaff, start, end)
ubiq1<-foverlaps(cand2, ubiq, type="any")
ubiq2<-subset(ubiq1,!ubiq1$POS %in% NA)
fisher.test(matrix(c(nrow(ubiq2), nrow(cand2)-nrow(ubiq2), 3075-nrow(ubiq2),nrow(genes)-nrow(ubiq2)-nrow(cand2)-3075), nrow=2))
ubiq2$i.ID

##Outlier CNVs
out<-fread("file:///home/aa/2alpine/predictors/CNV/Outlier99CNVs.txt",h=T)
genes<-fread("/home/aa/Desktop/references/lyrataV2/genesIDsLyV2.gff",h=T)
genes$ID<-substr(genes$ID,4,12)
setkey(out, scaff, start, end)
out1<-foverlaps(genes, out, type="any")
out2<-subset(out1,!out1$POS %in% NA)
write.table(out2,"Outlier99CNVs_genes.txt",col.names = T,row.names = F, quote = F)
out3<-as.data.frame(table(out2$ID))
hist(out3$Freq,breaks = 20)
# GO enrichment
dict<-fread("/home/aa/Desktop/references/lyrataV2/functions/ALATdict.txt")
library("biomaRt")
library(topGO)
library(data.table)
mart <- biomaRt::useMart(biomart = "plants_mart",dataset = "athaliana_eg_gene", host = 'plants.ensembl.org')
GTOGO <- biomaRt::getBM(attributes = c( "ensembl_gene_id", "go_id"), mart = mart)
GTOGO <- GTOGO[GTOGO$go_id != '',]
geneID2GO <- by(GTOGO$go_id,GTOGO$ensembl_gene_id,function(x) as.character(x))
all.genes <- sort(unique(as.character(GTOGO$ensembl_gene_id)))
se<-subset(dict,dict$AL %in% out2$i.ID & !dict$AT %in% "nnn")
sel<-substr(se$AT,1,9)
int.genes <- factor(as.integer(all.genes %in% sel))
names(int.genes) = all.genes
go.obj <- new("topGOdata", ontology='BP', allGenes = int.genes, annot = annFUN.gene2GO, gene2GO = geneID2GO,nodeSize=50) ## 
resultsFe <- runTest(go.obj, algorithm = "elim", statistic = "fisher") #More conservative
allRes <- GenTable(go.obj, elimFisher = resultsFe, orderBy = "elimFisher", ranksOf = "elimFisher", topNodes = 100)
a<-subset(allRes, as.numeric(allRes$Significant) > 3  & as.numeric(allRes$elimFisher) <= 0.05) 
write.table(file=paste("genes_cnv99_topgo_node50.txt",sep=""),a,sep="\t",row.names=F)
write.table(file=paste("65ATgenes_cnv99.txt",sep=""),sel,sep="\t",row.names=F, quote = F)

# aVT candiates overlapping with CNV candidates??
cand<-fread("file:///home/aa/alpine/arenosaGenome/selScans/bayPass_quartetFst/Tatry2x/genes_genic.txt",h=F)
ubiq<-fread("/home/aa/2alpine/predictors/CNV/Outlier99CNVs.txt",h=T)
genes<-fread("/home/aa/Desktop/references/lyrataV2/genesIDsLyV2.gff",h=T)
genes$ID<-substr(genes$ID,4,12)
cand2<-subset(genes,genes$ID %in% cand$V1)
setkey(ubiq, scaff, start, end)
ubiq1<-foverlaps(cand2, ubiq, type="any")
ubiq2<-subset(ubiq1,!ubiq1$POS %in% NA)
fisher.test(matrix(c(nrow(ubiq2), nrow(cand2)-nrow(ubiq2), 3075-nrow(ubiq2),nrow(genes)-nrow(ubiq2)-nrow(cand2)-3075), nrow=2))
write.table(file=paste("9genes_SNPTatry_cnv99.txt",sep=""),ubiq2,sep="\t",row.names=F, quote = F)

#length difference candidate vs. all??
library(ggpubr)
hist(ubiq$length)
hist(cnv$length)
sel <- cbind("Candidate CNVs",ubiq$length)
non <- cbind("All CNVs",cnv$length[cnv$length<20000])
all<-as.data.frame(rbind(sel,non))
colnames(all)<-c('CNVs',"length")
all$length<-as.numeric(as.character(all$length))
my_comparisons <- list( c("Candidate CNVs", "All CNVs") )
pdf("CNVsLength.pdf",width = 4,height = 4,pointsize = 20)
ggviolin(all, x = 'CNVs', y = "length", fill = "CNVs",
         palette = c("red", "grey50"),
         add = "boxplot", add.params = list(fill = "white")) +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif") + 
  stat_compare_means(label.y = max(all$length)+50)
dev.off()

# N of 'de-novo' foothill
out<-fread("file:///home/aa/2alpine/predictors/CNV/Outlier99CNVs.txt",h=T)
ubiq<-subset(out, out$AFfoot == 1)
cnvA<-subset(cnv, cnv$AFfoot == 1)
cnvF<-subset(cnv, cnv$AFalp == 1)
nrow(cnvA)/nrow(cnv)
  

##########Details pape
a<-read.table("file:///home/aa/lyrataWGD/PopKeyAll1000.csv",h=T,sep=",")
aa<-subset(a,a$ALL == 1)
aa$pop<-substr(aa$Samples,1,3)
unique(aa$pop)

b<-read.table("/home/aa/lyrataWGD/PopKey_lyrenosa.heatmap.csv",h=T,sep=",")
bb<-subset(b,b$ALL == 1)


#########tm
setwd("/home/aa/alpine/treemix/lyrataWGDFinal/")
library(data.table)
aa<-c("ZEP_tm.table","VLH_tm.table","VLA_tm.table","TEM_tm.table","SWA_tm.table","SUB_tm.table","STD_tm.table","SCT_tm.table","PER_tm.table","PEK_tm.table","OSL_tm.table","MOD_tm.table","JOH_tm.table","HOC_tm.table","CHO_tm.table","BRE_tm.table","BAB_tm.table","GUN_tm.table")
for (tm in aa) {
  aaa<-read.table(tm,h=F)
  aaa<-subset(aaa,!aaa$V3 %like% ",")
  write.table(x = aaa,file = paste("final/",tm,sep=""),quote = F,row.names = F,col.names = F,sep="\t")
}



##### TOM PicMin
#Add annotation
setwd("/home/aa/lyrataWGD/picmin/")
library(data.table)
cds<-fread("/home/aa/Desktop/references/lyrataV2/genesIDsLyV2.gff",h=T)
a<-fread("PicMinResults.csv",h=T)
a<-a[ order(a[,2]), ]
a1<-subset(a,a$pooled_q<=0.1)
a1$scaff<-paste(a1$redundan,substr(a1$scaffold,10,10),sep="_")
a1$end<-a1$start+1000
setkey(cds, scaff, start, end)
aaa<-foverlaps(a1, cds, type="any")

##Compare to my meiosis and cyclin
my<-fread("file:///home/aa/lyrataWGD/fstScan/outGenes_SNP0.99_Gene0.25.ALAT.ann.doubleAAAL.txt",h=T)
overl<- subset(my,my$AL %in% substr(aaa$ID,4,12))
overl2<- subset(aaa,substr(aaa$ID,4,12) %in% my$AL)

overl<-overl[ order(overl[,1]), ]
overl2<-overl2[ order(overl2[,15]), ]
overl2<-overl2[ order(overl2[,9]), ]

overl2<-overl2[!duplicated(overl2[,c('ID')]),]
tot<-cbind(overl,overl2)
tot<-tot[ order(tot[,23]), ] ##18

##Compare to any candidate
my<-fread("/home/aa/lyrataWGD/fstScan/outGenes_SNP0.99_Gene0.25.ALAT.ann.doubleAAAL.txt",h=T)
overl<- subset(my,my$AL %in% substr(aaa$ID,4,12))
overl2<- subset(aaa,substr(aaa$ID,4,12) %in% my$AL)
overl<-overl[ order(overl[,1]), ]
overl2<-overl2[ order(overl2[,15]), ]
overl2<-overl2[ order(overl2[,9]), ]

overl2<-overl2[!duplicated(overl2[,c('ID')]),]
tot<-cbind(overl,overl2)
tot<-tot[ order(tot[,23]), ]

write.table(file=paste("/home/aa/lyrataWGD/picmin/25GenesS01G25PicMin0.1.txt",sep=""),tot,sep="\t",row.names=F, quote = F)













################# ADEGENET ###################
setwd("/home/aa/lyrataWGD/adegenet/adaptiveGenes/")
# a function for conversion from vcfR object to genlight in tetraploids
source("adegenet_functions.r")
#install gcc-fortran, gcc-c++, zlib-devel
library(adegenet)
library(vcfR)
library(StAMPP)
library(ggplot2)
library(gplots)
library(data.table)
# ---------------------------------------------------------
# IMPORT SNP data from VCF - necessary to change path using setwd()
#vcf <- read.vcfR("lyrenosa_snp_raw.merged.HIS0.7.vcf")    # nrows=xxx - limit
vcf <- read.vcfR("allOut_all_raw.merged.12cand.vcf.gz")    # nrows=xxx - limit

head(vcf)  
# convert to genlight 
aa.genlight <- vcfR2genlight.tetra(vcf)
locNames(aa.genlight) <- paste(vcf@fix[,1],vcf@fix[,2],sep="_") 
#indNames(aa.genlight) <- c("a4al1_ING01", "a4al1_ING02", "a4al1_ING03", "a4al1_ING04", "a4al1_ING05", "a4al1_ING06", "a4al1_ING07", "a4al1_ING08")
name<-read.table("file:///home/aa/lyrataWGD/adegenet/namesAdegenetGW.txt")
indNames(aa.genlight) <- name$V1
aa.genlight1<-aa.genlight[c(1,2,3,5,6,7,8,9,10,11,12,13,14,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,47,48,49,50,51,52,53,54,55,77,78,79,80,81,82,83,84,110,111,112,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,162,163,164,165,166,167,168,169,182,183,184,185,186,187,188,189,329,330,331,332,333,346,347,348,349,350,351,352,436,437,438,439,440,441,477,478,479,482,483,484,485,486,487,488,489,490,491,492,493,494,523,524,525,526,527,528,529,530,531,534,535,536,537,542,543,544,545,546,547,557,558,559,560,561,562,563,564,565,566,571,572,573,574,575,576,577,578,579,600,601,602,603,604,605,606,607,608,610,611,612,613,614,615,616,617,618,627,628,645,646,647,648,649,650,651,652,664,665,666,667,668,669,670,671,672,673,674,675,677,678,679,680,681,682,683,685,686,687,688,689,690,691,692,693,694,695,696,712,713,714,715,716,717,718,719,720,721,738,739,740,741,742,743,757,758,759,760,761,762,763,764,765,766,767,768,769,772,773,774,775,776,777,778,779,796,797,798,799,800,801,802,803,804,805,806,807,808,809,810,811,812,813,814,815,816,817,818,819,820,821,822,823,824,825,826,827,828,831,832,833,834,835,848,849,850,851,852,856,857,858,859,860,867,893,894,895,896,897,928,929,930,931,932,933,934,935,936,937,938,939,940,941,942,943,944,953,954,955,956,957,958,1016,1017,1018,1019,1020,1021,1022,1023,1024,1025,1026,1027,1028,1029,1030,1031,1032,1033,1035)]

aa.genlight1<-aa.genlight[c(1,2,3,5,6,7,8,9,10,11,12,13,14,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,47,48,49,50,51,52,53,54,55,77,78,79,80,81,82,83,84,110,111,112,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,162,163,164,165,166,167,168,169,182,183,184,185,186,187,188,189,329,330,331,332,333,346,347,348,349,350,351,352,357,358,359,360,361,362,363,364,436,437,438,439,440,441,477,478,479,482,483,484,485,486,487,488,489,490,491,492,493,494,503,504,505,523,524,525,526,527,528,529,530,531,534,535,536,537,542,543,544,545,546,547,557,558,559,560,561,562,563,564,565,566,567,568,569,570,571,572,573,574,575,576,577,578,579,600,601,602,603,604,605,606,607,608,610,611,612,613,614,615,616,617,618,627,628,645,646,647,648,649,650,651,652,664,665,666,667,668,669,670,671,672,673,674,675,677,678,679,680,681,682,683,685,686,687,688,689,690,691,692,693,694,695,696,712,713,714,715,716,717,718,719,720,721,731,732,733,734,735,736,737,738,739,740,741,742,743,757,758,759,760,761,762,763,764,765,766,767,768,769,772,773,774,775,776,777,778,779,796,797,798,799,800,801,802,803,804,805,806,807,808,809,810,811,812,813,814,815,816,817,818,819,820,821,822,823,824,825,826,827,828,831,832,833,834,835,836,837,838,839,840,841,842,843,844,845,848,849,850,851,852,856,857,858,859,860,867,868,869,870,871,872,873,874,875,893,894,895,896,897,928,929,930,931,932,933,934,935,936,937,938,939,940,941,942,943,944,945,946,947,948,949,950,951,952,953,954,955,956,957,958,1016,1017,1018,1019,1020,1021,1022,1023,1024,1025,1026,1027,1028,1029,1030,1031,1032,1033,1035)]

pop(aa.genlight1)<-substr(indNames(aa.genlight1),1,3)  # add pop names: here pop names are first 3 chars of ind name
indNames(aa.genlight1)
ploidy(aa.genlight1)
pop(aa.genlight1)
# save the matrix and load it back


#write.table(as.matrix(aa.genlight), file="arenosa_spn_raw.fourfold.filtered.8ChromSubset.txt", row.names=T, col.names=T)
#aa.genlight <-new ("genlight",(read.table(file="all_BI_PASS_BP_4fold_matrix.txt")))
##################
#  various checks
#################
#how many particular classes of genotypes are in my dataset (i.e. numbers 0 1 2 3 4)
#table(as.matrix(aa.genlight))
##2. Final
#popNames(aa.genlight)
#subset populations
#obj <- aa.genlight[pop=c(2,11,18,38,44,62, 82, 86,91,92,97,110,116,122,145,146,156,125,132)]
#popNames(obj)
#indNames(obj)
#subset individuals
#aa.genlight1<-aa.genlight[c(1,3:156)]
head(aa.genlight1)
popNames(aa.genlight1)
indNames(aa.genlight1)
his<-read.table("HIS0.7Syn17MayNAMES.intervals")
his<-subset(his,!his$V1 %in% c("scaffold_2:17708821","scaffold_3:7249035","scaffold_3:7249358","scaffold_4:22848787","scaffold_4:22848802","scaffold_4:22848974","scaffold_4:22849516","scaffold_5:16015861","scaffold_5:16015873","scaffold_7:10777960","scaffold_8:22400532","scaffold_8_22400861"))
prot<-unique(his$V3)
for (p in prot) { # p="AT4G18490"
    coord<-which(rowSums(his == p) > 0)
    ####Separately for each protein
    fi<-which(aa.genlight1@loc.names %like% his$V1[coord[1]])
    la<-which(aa.genlight1@loc.names %like% his$V1[coord[length(coord)]])
    
    obj <- aa.genlight1[,fi:la]
    #head(obj)
    #GL plot-takes some time
    #glPlot(aa.genlight)  # takes some time 
    glPlot(obj)  # takes some time 
    
    
    aa.D.ind <- stamppNeisD(obj, pop = FALSE)  # Nei's 1972 distance between indivs
    stamppPhylip(aa.D.ind, file=paste("Neis_distance",p,"IndivDiploids.phy.dst",sep="_")) # export matrix - for SplitsTree
    aa.D.pop <- stamppNeisD(obj, pop = TRUE)   # Nei's 1972 distance between pops
    stamppPhylip(aa.D.pop, file=paste("Neis_distance",p,"PopsDiploids.phy.dst",sep="_")) # export matrix - for SplitsTree
}
    
    
    
    #aa.genlight<-obj
    ###  1.PCA
    #this is only needed after subsetting 
    #toRemove <- is.na(glMean(obj1, alleleAsUnit = FALSE)) # TRUE where NA
    #which(toRemove) # position of entirely non-typed loci
    #b <- obj1[, !toRemove]
    #calculate PCA
    pca.1 <- glPcaFast(obj, nf=30)   # retain first 300 axes (for later use in find.clusters); faster 
    ### PLOTTING
    scatter(pca.1, posi="topright")  # plot scatter with the indiv labels
    #loadingplot(pca.1)
    # just to see pops coloured in a palette
    col <- funky(19)
    s.class(pca.1$scores, pop(obj),  xax=1, yax=2, col=transp(col,.6))
    s.class(pca.1$scores, pop(obj),  xax=1, yax=3, col=transp(col,.6))
    pca.1$eig[1]/sum(pca.1$eig)
    pca.1$eig[2]/sum(pca.1$eig)
    pca.1$eig[3]/sum(pca.1$eig)
    pca.1$eig[4]/sum(pca.1$eig)
    #pca.1$eig[5]/sum(pca.1$eig)
    # save nice figs
    pdf (paste("PCA_",p,".pdf",sep="_"), width=14, height=7, pointsize = 24)
    par(mfrow=c(1,2))
    g1 <- s.class(pca.1$scores, pop(obj),  xax=1, yax=2, col=transp(col,.6),cstar = 0,cellipse = 0,clabel = 0.8,cpoint = 4,grid = F,addaxes = T,sub = paste("Axis 1 = ",round(pca.1$eig[1]/sum(pca.1$eig),2),", axis 2 = ",round(pca.1$eig[2]/sum(pca.1$eig),2),sep=""),csub = 1)
    g2 <- s.class(pca.1$scores, pop(obj),  xax=1, yax=2, col=transp(c("red","blue","blue","blue","blue","blue","blue","red","blue","blue","red","blue","red","red","blue","blue","blue","red","red"),.6),cstar = 0,cellipse = 0,clabel = 0,cpoint = 4,grid = F,addaxes = F,sub = paste("Axis 1 = ",round(pca.1$eig[1]/sum(pca.1$eig),2),", axis 2 = ",round(pca.1$eig[2]/sum(pca.1$eig),2),sep=""),csub = 1)
    
    g1 <- s.class(pca.1$scores, pop(obj),  xax=1, yax=2, col=transp(col,.6),cstar = 0,cellipse = 0,clabel = 0.8,cpoint = 4,grid = F,addaxes = T,sub = paste("Axis 1 = ",round(pca.1$eig[1]/sum(pca.1$eig),2),", axis 2 = ",round(pca.1$eig[2]/sum(pca.1$eig),2),sep=""),csub = 1)
    g2 <- s.class(pca.1$scores, pop(obj),  xax=1, yax=2, col=transp(c("red","blue","blue","blue","blue","blue","blue","red","blue","blue","red","blue","red","red","blue","blue","blue","red","red"),.6),cstar = 0,cellipse = 0,clabel = 0.8,cpoint = 4,grid = F,addaxes = T,sub = paste("Axis 1 = ",round(pca.1$eig[1]/sum(pca.1$eig),2),", axis 2 = ",round(pca.1$eig[2]/sum(pca.1$eig),2),sep=""),csub = 1)

    g1 <- s.class(pca.1$scores, pop(obj),  xax=1, yax=3, col=transp(col,.6),cstar = 0,cellipse = 0,clabel = 0.8,cpoint = 4,grid = F,addaxes = T,sub = paste("Axis 1 = ",round(pca.1$eig[1]/sum(pca.1$eig),2),", axis 3 = ",round(pca.1$eig[3]/sum(pca.1$eig),2),sep=""),csub = 1)
    g2 <- s.class(pca.1$scores, pop(obj),  xax=1, yax=3, col=transp(c("red","blue","blue","blue","blue","blue","blue","red","blue","blue","red","blue","red","red","blue","blue","blue","red","red"),.6),cstar = 0,cellipse = 0,clabel = 0.8,cpoint = 4,grid = F,addaxes = T,sub = paste("Axis 1 = ",round(pca.1$eig[1]/sum(pca.1$eig),2),", axis 3 = ",round(pca.1$eig[3]/sum(pca.1$eig),2),sep=""),csub = 1)
    #  g1 <- s.class(pca.1$scores, pop(obj),  xax=1, yax=4, col=transp(col,.6),cstar = 0,cellipse = 0,clabel = 0.8,cpoint = 4,grid = F,addaxes = T,sub = paste("Axis 1 = ",round(pca.1$eig[1]/sum(pca.1$eig),2),", axis 4 = ",round(pca.1$eig[4]/sum(pca.1$eig),2),sep=""),csub = 1)
    #  g2 <- s.class(pca.1$scores, pop(obj),  xax=1, yax=4, col=transp(c("red","blue","blue","blue","blue","blue","blue","red","blue","blue","red","blue","red","red","blue","blue","blue","red","red"),.6),cstar = 0,cellipse = 0,clabel = 0.8,cpoint = 4,grid = F,addaxes = T,sub = paste("Axis 1 = ",round(pca.1$eig[1]/sum(pca.1$eig),2),", axis 4 = ",round(pca.1$eig[4]/sum(pca.1$eig),2),sep=""),csub = 1)
    dev.off()
  })}
    
    ### 2. NJ
    #to calculate genetic distances and visulize
    aa.D.pop <- stamppNeisD(obj, pop = TRUE)   # Nei's 1972 genetic distance between pops
    #plot tree - neighbour joining tree
    pdf (paste("NJtree_",p,".pdf",sep=""), width=7, height=7)
    plot(nj(aa.D.pop))
    dev.off()
    ### Calculate Nei's distances between individuals/pops
    # ---------------------------------------------------
    aa.D.ind <- stamppNeisD(aa.genlight, pop = FALSE)  # Nei's 1972 distance between indivs
    stamppPhylip(aa.D.ind, file=paste("Neis_distance",p,"Indiv.phy.dst",sep="_")) # export matrix - for SplitsTree
    aa.D.pop <- stamppNeisD(aa.genlight, pop = TRUE)   # Nei's 1972 distance between pops
    stamppPhylip(aa.D.pop, file=paste("Neis_distance",p,"Pops.phy.dst",sep="_")) # export matrix - for SplitsTree
    # plot heatmap of the population distance matrix
    colnames(aa.D.pop)<-rownames(aa.D.pop) # name the rows of a matrix  
    pdf(paste("HeatmapDistance_",p,".pdf",sep=""), width=10,height=10)
    heatmap.2(aa.D.pop, trace="none", cexRow=0.7, cexCol=0.7)
    dev.off()
  })}
#The end of loop




########## ADEGENET NEUTRAL GENOME ###############
setwd("/home/aa/lyrataWGD/adegenet/")
# a function for conversion from vcfR object to genlight in tetraploids
source("adaptiveGenes/adegenet_functions.r")
#install gcc-fortran, gcc-c++, zlib-devel
library(adegenet)
library(vcfR)
library(StAMPP)
library(ggplot2)
library(gplots)
# ---------------------------------------------------------
# IMPORT SNP data from VCF - necessary to change path using setwd()
vcf <- read.vcfR("arabidopsisDemo.snps.bipassed.dp8.mffg0.2.vcf.gz",nrows = 50000)    # nrows=xxx - limit
head(vcf)  
# convert to genlight 
aa.genlight <- vcfR2genlight.tetra(vcf)
#locNames(aa.genlight) <- paste(vcf@fix[,1],vcf@fix[,2],sep="_") 
name<-read.table("file:///home/aa/lyrataWGD/adegenet/namesAdegenetGW.txt")
indNames(aa.genlight) <- name$V1
#indNames(aa.genlight)
ploidy(aa.genlight)
# save the matrix and load it back
##################
#  various checks
#################
#how many particular classes of genotypes are in my dataset (i.e. numbers 0 1 2 3 4)
#table(as.matrix(aa.genlight))
##2. Final
#popNames(aa.genlight)
#subset populations
#obj <- aa.genlight[pop=c(2,11,18,38,44,62, 82, 86,91,92,97,110,116,122,145,146,156,125,132)]
#popNames(obj)
#indNames(obj)
#subset GOOD individuals
aa.genlight1<-aa.genlight[c(1,2,3,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,187,188,189,195,196,197,198,199,200,201,202,203,204,205,206,207,208,209,210,211,212,213,214,215,216,217,218,219,220,221,222,223,224,225,226,227,228,229,230,231,232,233,234,235,236,237,238,239,240,241,242,243,244,245,246,247,248,249,250,251,252,253,254,255,256,257,258,259,260,261,262,263,264,265,266,267,268,269,270,271,272,273,274,275,276,277,278,279,280,281,282,283,284,285,286,287,288,289,290,291,292,293,294,295,296,297,298,299,300,301,302,303,304,305,306,307,308,309,310,311,312,313,314,315,316,317,318,319,320,321,322,323,324,325,326,327,328,329,330,331,332,333,334,335,336,337,338,339,340,341,342,343,344,345,346,347,348,349,350,351,352,353,354,355,356,357,358,359,360,361,362,363,364,365,366,367,368,369,370,371,372,373,374,375,376,377,378,379,380,381,383,384,385,386,387,388,389,390,391,392,393,394,395,396,397,398,399,400,401,402,403,404,405,406,407,408,409,410,411,412,413,414,415,416,417,418,419,420,421,422,425,426,427,428,429,430,431,432,433,434,435,436,437,438,439,440,441,442,443,444,445,446,447,448,449,450,451,453,454,455,456,457,458,459,460,461,462,463,464,465,466,467,468,469,470,471,472,473,474,475,476,477,478,479,480,481,482,483,484,485,486,487,488,489,490,491,492,493,494,495,496,497,498,499,500,501,502,503,504,505,520,521,522,523,524,525,526,527,528,529,530,531,533,534,535,536,537,538,539,540,541,542,543,544,545,546,547,548,549,550,551,552,553,554,555,556,557,558,559,560,561,562,563,564,565,566,567,568,569,570,571,572,573,574,575,576,577,578,579,580,581,582,583,584,585,586,587,588,589,590,591,592,593,594,595,596,597,598,599,600,601,602,603,604,605,606,607,608,610,611,612,613,614,615,616,617,618,619,620,621,622,623,624,625,626,627,628,629,630,631,632,633,634,635,636,637,638,639,640,641,642,643,644,645,646,647,648,649,650,651,652,653,654,655,656,657,658,659,660,661,662,663,664,665,666,667,668,669,670,671,672,673,674,675,677,678,679,680,681,682,683,685,686,687,688,689,690,691,692,693,694,695,696,699,700,701,702,703,704,705,706,707,708,709,710,711,712,713,714,715,716,717,718,719,720,721,722,723,724,725,726,727,728,729,730,731,732,733,734,735,736,737,738,739,740,741,742,743,744,745,746,747,748,749,750,751,752,753,754,755,756,757,758,759,760,761,762,763,764,765,766,767,768,769,772,773,774,775,776,777,778,779,780,781,782,783,784,785,786,787,788,789,790,791,792,793,794,795,796,797,798,799,800,801,802,803,804,805,806,807,808,809,810,811,812,813,814,815,816,817,818,819,820,821,822,823,824,825,826,827,828,831,832,833,834,835,836,837,838,839,840,841,842,843,844,845,846,847,848,849,850,851,852,853,854,855,856,857,858,859,860,861,862,863,864,865,866,867,868,869,870,871,872,873,874,875,876,877,878,879,880,881,882,883,884,885,886,887,888,889,890,891,893,894,895,896,897,898,899,900,901,902,903,904,905,906,907,908,909,910,911,912,913,914,915,916,917,918,919,920,921,922,923,924,925,926,927,928,929,930,931,932,933,934,935,936,937,938,939,940,941,942,943,944,945,946,947,948,949,950,951,952,953,954,955,956,957,958,964,965,966,967,968,969,970,971,974,975,976,977,978,979,980,981,982,983,984,985,986,987,988,989,990,991,992,993,994,995,996,997,999,1000,1001,1002,1003,1004,1005,1006,1007,1008,1009,1010,1011,1012,1013,1014,1015,1016,1017,1018,1019,1020,1021,1022,1023,1024,1025,1026,1027,1028,1029,1030,1031,1032,1033,1035,1036,1037,1038,1039,1040,1041,1042,1043)]

#subset GOOD diploids
aa.genlight1<-aa.genlight[c(1,2,3,5,6,7,8,9,10,11,12,13,14,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,47,48,49,50,51,52,53,54,55,77,78,79,80,81,82,83,84,110,111,112,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,152,153,154,155,156,157,158,159,162,163,164,165,166,167,168,169,180,181,182,183,184,185,186,187,188,189,214,215,216,217,218,219,220,221,222,223,224,225,226,227,228,229,230,231,232,233,234,235,236,248,249,250,253,254,255,256,257,258,259,260,261,262,263,264,265,266,267,268,269,270,271,272,273,274,275,276,277,278,279,280,281,282,283,284,285,286,287,288,289,290,291,292,293,294,295,296,297,298,299,300,301,302,303,304,305,306,307,308,309,310,311,312,313,314,315,316,317,318,319,320,321,322,323,324,325,326,329,330,331,332,333,334,335,336,337,338,339,340,341,342,343,344,345,346,347,348,349,350,351,352,353,354,355,356,365,366,367,368,389,390,391,392,393,394,395,396,397,398,399,436,437,438,439,440,441,477,478,479,482,483,484,485,486,487,488,489,490,491,492,493,494,523,524,525,526,527,528,529,530,531,534,535,536,537,542,543,544,545,546,547,557,558,559,560,561,562,563,564,565,566,571,572,573,574,575,576,577,578,579,580,581,582,583,584,585,586,587,588,589,600,601,602,603,604,605,606,607,608,610,611,612,613,614,615,616,617,618,627,628,645,646,647,648,649,650,651,652,664,665,666,667,668,669,670,671,672,673,674,675,677,678,679,680,681,682,683,685,686,687,688,689,690,691,692,693,694,695,696,712,713,714,715,716,717,718,719,720,721,738,739,740,741,742,743,757,758,759,760,761,762,763,764,765,766,767,768,769,772,773,774,775,776,777,778,779,796,797,798,799,800,801,802,803,804,805,806,807,808,809,810,811,812,813,814,815,816,817,818,819,820,821,822,823,824,825,826,827,828,829,830,831,832,833,834,835,848,849,850,851,852,856,857,858,859,860,867,893,894,895,896,897,928,929,930,931,932,933,934,935,936,937,938,939,940,941,942,943,944,953,954,955,956,957,958,1007,1016,1017,1018,1019,1020,1021,1022,1023,1024,1025,1026,1027,1028,1029,1030,1031,1032,1033,1035)]

#subset GOOD diploids without halleri
aa.genlight1<-aa.genlight[c(1,2,3,5,6,7,8,9,10,11,12,13,14,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,47,48,49,50,51,52,53,54,55,77,78,79,80,81,82,83,84,110,111,112,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,162,163,164,165,166,167,168,169,182,183,184,185,186,187,188,189,329,330,331,332,333,346,347,348,349,350,351,352,436,437,438,439,440,441,477,478,479,482,483,484,485,486,487,488,489,490,491,492,493,494,523,524,525,526,527,528,529,530,531,534,535,536,537,542,543,544,545,546,547,557,558,559,560,561,562,563,564,565,566,571,572,573,574,575,576,577,578,579,600,601,602,603,604,605,606,607,608,610,611,612,613,614,615,616,617,618,627,628,645,646,647,648,649,650,651,652,664,665,666,667,668,669,670,671,672,673,674,675,677,678,679,680,681,682,683,685,686,687,688,689,690,691,692,693,694,695,696,712,713,714,715,716,717,718,719,720,721,738,739,740,741,742,743,757,758,759,760,761,762,763,764,765,766,767,768,769,772,773,774,775,776,777,778,779,796,797,798,799,800,801,802,803,804,805,806,807,808,809,810,811,812,813,814,815,816,817,818,819,820,821,822,823,824,825,826,827,828,831,832,833,834,835,848,849,850,851,852,856,857,858,859,860,867,893,894,895,896,897,928,929,930,931,932,933,934,935,936,937,938,939,940,941,942,943,944,953,954,955,956,957,958,1007,1016,1017,1018,1019,1020,1021,1022,1023,1024,1025,1026,1027,1028,1029,1030,1031,1032,1033,1035)]


head(aa.genlight1)
pop(aa.genlight1)<-substr(indNames(aa.genlight1),1,3)  # add pop names: here pop names are first 3 chars of ind name
popNames(aa.genlight1)
#indNames(aa.genlight1)
obj1<-aa.genlight1
#his<-read.table("HIS0.7Names.intervals")
#prot<-unique(his$V3)
#prot<-prot[c(1:13)]
    glPlot(obj1)  # takes some time 
    #aa.genlight<-obj
    ###  1.PCA
    #this is only needed after subsetting 
#    toRemove <- is.na(glMean(obj1, alleleAsUnit = FALSE)) # TRUE where NA
#    which(toRemove) # position of entirely non-typed loci
#    obj <- obj1[, !toRemove]
    #calculate PCA
#    pca.1 <- glPcaFast(obj, nf=30)   # retain first 300 axes (for later use in find.clusters); faster 
    ### PLOTTING
#    scatter(pca.1, posi="topright")  # plot scatter with the indiv labels
### 2. NJ
#to calculate genetic distances and visulize
#aa.D.pop <- stamppNeisD(obj, pop = TRUE)   # Nei's 1972 genetic distance between pops
#plot tree - neighbour joining tree
#pdf (paste("NJtree_GW.pdf",sep=""), width=7, height=7)
#plot(nj(aa.D.pop))
#dev.off()
### Calculate Nei's distances between individuals/pops
# ---------------------------------------------------
aa.D.pop <- stamppNeisD(obj1, pop = TRUE)   # Nei's 1972 distance between pops
stamppPhylip(aa.D.pop, file=paste("Neis_distanceGWPopsDiploidNoHalleri.phy.dst",sep="_")) # export matrix - for SplitsTree
aa.D.ind <- stamppNeisD(obj1, pop = FALSE)  # Nei's 1972 distance between indivs
stamppPhylip(aa.D.ind, file=paste("Neis_distanceGWIndivDiploidNoHalleri.phy.dst",sep="_")) # export matrix - for SplitsTree

####




##### BARPLOTS - N, S, U, R
library(data.table)
library(dplyr)
library(gplots)
library(ggplot2)
setwd("/home/aa/lyrataWGD/fstScan/ProtHistory/")
popkey<-fread("/home/aa/lyrataWGD/fstScan/haplotypeMaps/PopKeyAll1000.csv",h=T)
namesAll<-read.table("file:///home/aa/lyrataWGD/fstScan/haplotypeMaps/17WGDcand_names.txt",h=T)
for (i in 1:nrow(namesAll)) { #  i=1
  i2<-namesAll[i,1]
  name<-namesAll[i,2]
  aa<-read.table(paste("../haplotypeMaps/",name,".HIS.arenosa.txt",sep = ""),h=T,sep="\t")
al<-read.table(paste("../haplotypeMaps/",name,".HIS.lyrata.txt",sep = ""),h=T,sep="\t")
both<-rbind(aa[,1:5],al[,1:5])
both<-both[ order(both[,2]), ]
bo<-both[!duplicated(both[,c('scaff','start')]),]
s<-nrow(subset(bo,bo$ann %in% "synonymous_variant"))/nrow(bo)
n<-nrow(subset(bo,bo$ann %in% "missense_variant"))/nrow(bo)
u<-nrow(subset(bo,bo$ann %in% "upstream_gene_variant" | bo$ann %in% "5_prime_UTR_variant"))/nrow(bo)
r<-(nrow(bo)-nrow(subset(bo,bo$ann %in% "synonymous_variant"))-nrow(subset(bo,bo$ann %in% "missense_variant")))/nrow(bo)-nrow(subset(bo,bo$ann %in% "upstream_gene_variant" | bo$ann %in% "5_prime_UTR_variant"))/nrow(bo)

print(name)
#print(nrow(subset(bo,bo$ann %in% "upstream_gene_variant" | bo$ann %in% "5_prime_UTR_variant")))
print(nrow(subset(bo,bo$ann %in% "synonymous_variant")))
print(nrow(bo))
# write.table(paste(name,s,n,u,r,sep = "\t"),"NSUR_proportionsXX.txt",append = T,row.names = F,col.names = F,quote = F)
}
##Barplot
a1<-read.table("NSUR_proportions.txt",h=F)
colnames(a1)<-c("V1","Nonsyn.","Upstream","Oth. regul.","Synon.")
pdf("NSURBarplot.pdf",width = 15,height = 5)
barplot(height = t(a1[,2:5]), col = c("grey25", "grey50","grey65", "grey98"), names.arg = a1$V1,legend.text = colnames(a1)[2:5],args.legend = list(x = "topleft", inset = c(0.45, -0.25)),las=2)
dev.off()

###More regul in cyclins?
prop.test(x = c(5,4,4),n=c(23,7,7) ,p=rep(0.324,3), alternative = "g")
mean(c(5,4,4)/c(23,7,7))

###More nonsyn in meiosis?
prop.test(x = c(4,11,17,25,2),n=c(9,30,35,53,8) ,p=rep(0.278,5), alternative = "g")
mean(c(4,11,17,25,2)/c(9,30,35,53,8))

###More synon in RNA?
prop.test(x = c(3,13,9),n=c(14,14,11) ,p=rep(0.269,3), alternative = "g")
mean(c(3,13,9)/c(14,14,11))


##### BARPLOTS - Haplotypes
library(data.table)
library(dplyr)
library(gplots)
library(ggplot2)
setwd("/home/aa/lyrataWGD/fstScan/ProtHistory/")
popkey<-fread("/home/aa/lyrataWGD/fstScan/haplotypeMaps/PopKeyAll1000.csv",h=T)
namesAll<-read.table("file:///home/aa/lyrataWGD/fstScan/haplotypeMaps/17WGDcand_names.txt",h=T)
for (i in 1:nrow(namesAll)) { #  i=1
  i2<-namesAll[i,1]
  name<-namesAll[i,2]
  tet<-read.table(paste("../haplotypeMaps/",name,".hapl.union.txt",sep = ""),h=T,sep="\t")
  tet1<-subset(tet,tet$col == 4 & !substr(tet$pop,1,2) %in% "LL" & !substr(tet$pop,1,2) %in% "XX")
  tetfreq<-mean(tet1$c2,na.rm = T)
  tea<-read.table(paste("../haplotypeMaps/",name,".hapl.arenosa.txt",sep = ""),h=T,sep="\t")
  tea1<-subset(tea,tea$col == 4 & !substr(tea$pop,1,2) %in% "LL" & !substr(tea$pop,1,2) %in% "XX")
  teafreq<-mean(tea1$c1,na.rm = T)
  tel<-read.table(paste("../haplotypeMaps/",name,".hapl.lyrata.txt",sep = ""),h=T,sep="\t")
  tel1<-subset(tel,tel$col == 4 & !substr(tel$pop,1,2) %in% "LL" & !substr(tel$pop,1,2) %in% "XX")
  telfreq<-mean(tel1$c1,na.rm = T)
  rec<-1-tetfreq-teafreq-telfreq
  write.table(paste(name,tetfreq,teafreq,telfreq,rec,sep = "\t"),"HaplTetrPops_proportions_inclKRK.txt",append = T,row.names = F,col.names = F,quote = F)
}

a1<-read.table("HaplTetrPops_proportions_inclKRK.txt",h=F)
colnames(a1)<-c("V1","Tetraploid.","Diploid - A. arenosa","Diploid - A. lyrata","Recombined")
summary(a1)
pdf("HaplTetrPopsBarplot_inclKRK.pdf",width = 15,height = 5)
barplot(height = t(a1[,2:5]), col = c("Blue", "tomato", "red","grey"), names.arg = a1$V1,legend.text = colnames(a1)[2:5],args.legend = list(x = "topleft", inset = c(0.45, -0.25)))
dev.off()
  
  
  ##Diploid arenosa 
namesAll<-read.table("file:///home/aa/lyrataWGD/fstScan/haplotypeMaps/17WGDcand_names.txt",h=T)
for (i in 1:nrow(namesAll)) { #  i=1
  i2<-namesAll[i,1]
  name<-namesAll[i,2]
  
  dia<-read.table(paste("../haplotypeMaps/",name,".hapl.arenosa.txt",sep = ""),h=T,sep="\t")
  dia1<-subset(dia,dia$colorText == "tomato")
  diafreq<-mean(dia1$c1,na.rm = T)
  
  rec<-1-diafreq
  write.table(paste(name,diafreq,rec,sep = "\t"),"HaplDiAPops_proportions.txt",append = T,row.names = F,col.names = F,quote = F)
}
a1<-read.table("HaplDiAPops_proportions.txt",h=F)
colnames(a1)<-c("V1","Diploid - A. arenosa","Recombined")
pdf("HaplDiAPopsBarplot.pdf",width = 15,height = 5)
barplot(height = t(a1[,2:3]), col = c( "tomato", "grey"), names.arg = a1$V1,legend.text = colnames(a1)[2:3],args.legend = list(x = "topleft", inset = c(0.45, -0.25)))
dev.off()
  

##Diploid lyrata 
namesAll<-read.table("file:///home/aa/lyrataWGD/fstScan/haplotypeMaps/17WGDcand_names.txt",h=T)
for (i in 1:nrow(namesAll)) { #  i=1
  i2<-namesAll[i,1]
  name<-namesAll[i,2]
  dia<-read.table(paste("../haplotypeMaps/",name,".hapl.lyrata.txt",sep = ""),h=T,sep="\t")
#  dia1<-subset(dia,pop == "STD" | pop == "POT"| pop == "PLE"| pop == "STR"| pop == "MAL"| pop == "OSL"| pop == "VLH")
  dia1<-subset(dia,dia$colorText == "red" & !substr(dia$pop,1,2) %in% "LL")
  diafreq<-mean(dia1$c1,na.rm = T)
  rec<-1-diafreq
  write.table(paste(name,diafreq,rec,sep = "\t"),"HaplDiLAllPops_proportions.txt",append = T,row.names = F,col.names = F,quote = F)
}
a1<-read.table("HaplDiLAllPops_proportions.txt",h=F)
colnames(a1)<-c("V1","Diploid - A. lyrata","Recombined")
pdf("HaplDiLAllPopsBarplot.pdf",width = 15,height = 5)
barplot(height = t(a1[,2:3]), col = c( "red", "grey"), names.arg = a1$V1,legend.text = colnames(a1)[2:3],args.legend = list(x = "topleft", inset = c(0.45, -0.25)))
dev.off()






################ BARPLOTS: source of SNPs within haplotype
####
library(data.table)
library(dplyr)
setwd("/home/aa/lyrataWGD/fstScan/haplotypeMaps/")
ar<-fread(paste('ALL.ann/ALL.table.recode.txt',sep=''),h=F,na.strings = "-9",nThread = 3)
or<-fread("ALL.ann/ALO.table.recode.txt",h=F,na.strings = "-9")
popkey<-fread("PopKeyAll1000SVdistr.csv",h=T)
popkeyo<-fread("PopKeyOutgroups.csv",h=T)
popkeyall<-rbind(popkey,popkeyo,use.names=FALSE)
names<-c(popkey[1:nrow(popkey),1])[[1]]
colnames(ar) <- c("pop","ploidy","scaff","start","AN","DP","ID","ann","AAS",names,"nan")
nameso<-c(popkeyo[1:nrow(popkeyo),1])[[1]]
colnames(or) <- c("pop","ploidy","scaff","start","AN","DP",nameso,"nan")
namesAll<-read.table("file:///home/aa/lyrataWGD/fstScan/haplotypeMaps/17WGDcand_namesSV.txt",h=T)
species<-read.csv("diploidSpecies.csv",h=F)

for (i in 1:nrow(namesAll)) { #  i=1
#  i2<-namesAll[i,1]
  name<-namesAll[i,2]
  hifs<-fread(paste(name,".HIS.union.txt",sep = ""),h=T)
  #ar<-fread(paste('ALL.ann/ALL.table.recode.txt',sep=''),h=F,na.strings = "-9",nThread = 3)
  #popkey<-fread("PopKeyAll1000.csv",h=T)
  prot <- name
  ar1<-subset(ar,ar$scaff %in% hifs$scaff & ar$start %in% hifs$start)
  or1<-subset(or,or$scaff %in% hifs$scaff & or$start %in% hifs$start)
  
  ar1$end<-ar1$start  
  or1$end<-or1$start  
  setkey(or1, scaff, start, end)
  artot<-foverlaps(ar1, or1, type="any")
  #artot<-artot[,7:ncol(artot)]
  
  arAC<-as.data.frame(cbind(artot$scaff,artot$i.start))
  arAN<-as.data.frame(cbind(artot$scaff,artot$i.start))
  popnames<-unique(substr(c(names,nameso),1,3))
  popnames<-popnames[!popnames %in% "XXX"]
  ploidy <- data.frame(pops = popnames, ploidy.lev = NA)
  for(i in 1:length(popnames)){  # i = 121
    afh1 <- artot %>% dplyr:: select(starts_with(popnames[i]))
    artot$ACh1 <- rowSums(afh1,na.rm = T)
    artot$NAh1 <- apply(is.na(afh1), 1, sum)
    if (subset(popkeyall,popkeyall$Samples %like% popnames[i])[1,3]==0) 
    {artot$ANh1<-(ncol(afh1)-artot$NAh1)*2
    ploidy[i,2] <- 2
    } else {artot$ANh1<-(ncol(afh1)-artot$NAh1)*4
    ploidy[i,2] <- 4}
    arAC$POP<-artot$ACh1
    arAN$POP<-artot$ANh1
    colnames(arAC)[which(colnames(arAC) %in% "POP")]<-popnames[i]
    colnames(arAN)[which(colnames(arAN) %in% "POP")]<-popnames[i]
  }
  
  ploidy<-ploidy[ order(ploidy[,1]), ]
  #repolarise - #TODO
  for (i in  1:nrow(arAC)){ # i=1
    d<-c()
    t<-c()
    for (pop in popnames) { #  pop="BAB"
      pppac <- arAC %>% dplyr:: select(starts_with(pop))
      pppan <- arAN %>% dplyr:: select(starts_with(pop))
      
      if (subset(popkeyall,popkeyall$Samples%like%pop)[1,3]==0) 
      {d<-c(d,as.numeric(as.character((pppac[i,1])[[1]]))/as.numeric(as.character((pppan[i,1])[[1]])))} else {t<-c(t,as.numeric(as.character((pppac[i,1])[[1]]))/as.numeric(as.character((pppan[i,1])[[1]])))}
    }
    meandip<-mean(d,na.rm = T)
    meantet<-mean(t,na.rm = T)
    if (meandip>meantet)
    {arAC[i,4:(length(arAC[i,])-1)]<-(arAN[i,4:(length(arAN[i,])-1)]-arAC[i,4:(length(arAC[i,])-1)])
    print(i)
    } else {}}
  
  #### Find for each species
  ACtot<-arAC[,1:2]
  for(i in 1:length(unique(species$V2))){  # i = 1
    bbb<- as.data.frame(arAC[, subset(species,species$V2 %in% unique(species$V2)[i])$V1])
    ACtot$POP<-rowSums(bbb)
    colnames(ACtot)[which(colnames(ACtot) %in% "POP")]<-unique(species$V2)[i]
  }
  
  #### Find for each lineage
  ACtotlin<-arAC[,1:2]
  for(i in 1:length(unique(species$V3))){  # i = 1
    bbb<- as.data.frame(arAC[, subset(species,species$V3 %in% unique(species$V3)[i])$V1])
    ACtotlin$POP<-rowSums(bbb)
    colnames(ACtotlin)[which(colnames(ACtotlin) %in% "POP")]<-unique(species$V3)[i]
  }
  
  write.table(ACtot,paste(name,".SVPerSpecies.txt",sep=""),quote = F,sep = "\t",row.names = F)
  write.table(ACtotlin,paste(name,".SVPerLineage.txt",sep=""),quote = F,sep = "\t",row.names = F)
  write.table(arAC,paste(name,".SVPerPop.txt",sep=""),quote = F,sep = "\t",row.names = F)
  }

###SV Heatmap candidate only
library(data.table)
library(dplyr)
library(gplots)

setwd("/home/aa/lyrataWGD/fstScan/haplotypeMaps/")
namesAll<-read.table("file:///home/aa/lyrataWGD/fstScan/haplotypeMaps/17WGDcand_names.txt",h=T)
for (i in 1:nrow(namesAll)) { #  i=1
  i2<-namesAll[i,1]
  name<-namesAll[i,2]
  s<-fread(paste("../haplotypeMaps/",name,".SVPerSpecies.txt",sep = ""),h=T)
  popnames<-colnames(s)[3:(length(colnames(s)))]
  my_palette <- colorRampPalette(c("tomato", "blue"))(n = 1000)
  id<-i2
  s5<-dplyr::select(s,AA, AL, AH, AC, AN, AP)
  s5$AA<-s5$AA/384
  s5$AL<-s5$AL/236
  s5$AH<-s5$AH/282
  s5$AN<-s5$AN/26
  s5$AC<-s5$AC/22
  s5$AP<-s5$AP/4
  
  df<-as.matrix(s5[,1:(length(s5))],rownames = s$V2)
  
  df1<-log10(df+0.01)
  #df1<-as.matrix(s5,rownames = s1$aas)
  #parallel
  pdf(paste("SVheatmaps/SVheatmap_",id,"_",name,".pdf",sep=""),height = 11,width = 9.3)
  heatmap.2(x = df1,dendrogram = "none",Colv="NA", Rowv="NA",key = F,col=my_palette,colsep= c(0,1,2,3,4,5,6),sepcolor= c("black"),sepwidth = c(0.02),trace="none",ColSideColors = c("red","red","red","red","red","red"),labCol=colnames(df),colCol= c("red","red","red","red","red","red"),offsetRow = c(0.05),offsetCol= c(0.05),cexRow = 1,cexCol = c(2),margins = c(7,7),lhei = c(0.2,9),lwid = c(0.03,0.8),xlab = paste(id,name,sep=" - ")) 
  dev.off()
}

aa<-read.table("ALL.SVPerSpecies.txt",h=T)
aa1<- aa[, 3:8]
aa$tot<-rowSums(aa1)
pdf("AFSofAdaptiveVariationAll.pdf",width = 18,height = 8)
hist(aa$tot,breaks = c(0,(seq(0,max(aa$tot)+1,1))),main = "All diploids",freq = T)
hist(aa$AA,breaks = c(0,(seq(0,max(aa$AA),1))),main = "All arenosa",freq = T)
hist(aa$AL,breaks = c(0,(seq(0,max(aa$AL),1))),main = "All lyrata",freq = T)
hist(aa$AH,breaks = c(0,(seq(0,max(aa$AH),1))),main = "All halleri",freq = T)
hist(aa$AN,breaks = c(0,(seq(0,max(aa$AN),1))),main = "All cebennensis",freq = T)
hist(aa$AC,breaks = c(0,(seq(0,max(aa$AC),1))),main = "All croatica",freq = T)
hist(aa$AP,breaks = c(0,(seq(0,max(aa$AP),1))),main = "All pedemontana",freq = T)

dev.off()

#write.table(ACtot,paste(name,".SVPerSpecies.txt",sep=""),quote = F,sep = "\t",row.names = F)

#Overall - in needed???
arAC<-read.table("ALL.SVPerPop.txt",h=T)
dip<-subset(ploidy,ploidy$ploidy.lev %in% 2)
arACdip<- arAC[, dip$pops]
arACdip$tot<-rowSums(arACdip)
hist(arACdip$tot,nclass = 100)
table(arACdip$tot)
write.table(arACdip,paste("diploids.SVPerPop.txt",sep=""),quote = F,sep = "\t",row.names = F)

#########Barplot
##### BARPLOTS
library(data.table)
library(dplyr)
library(gplots)
library(ggplot2)
library(RColorBrewer)
#setwd("/home/aa/lyrataWGD/fstScan/ProtHistory/")
#popkey<-fread("/home/aa/lyrataWGD/fstScan/haplotypeMaps/PopKeyAll1000.csv",h=T)
#namesSNPs<-read.table("file:///home/aa/lyrataWGD/fstScan/haplotypeMaps/17WGDcand_names.txt",h=T)


his<-read.table("HIS0.7Syn17MayNAMES.intervals",h=T)
table(his$name)
sv<-read.table("ALL.SVPerSpecies.txt",h=T)
his<-cbind(his,sv)
prot<-unique(his$name)
#Incl.Singletons
for (p in prot) { # p="CYCA2_3"
    coord<-which(rowSums(his == p) > 0)
    obj <- his[coord[1]:coord[length(coord)],]
    ar<-nrow(subset(obj,obj$AA > 0 & obj$AL == 0 & obj$AH == 0 & obj$AC == 0 & obj$AP == 0 & obj$AN == 0))/nrow(obj) #arenosa only
    ly<-nrow(subset(obj,obj$AL > 0 & obj$AA == 0 & obj$AH == 0 & obj$AC == 0 & obj$AP == 0 & obj$AN == 0))/nrow(obj) #lyrata only
    nla<-nrow(subset(obj,obj$AA == 0 & obj$AL == 0 & (obj$AH > 0 | obj$AC > 0 | obj$AP > 0 | obj$AN > 0)))/nrow(obj) #outgroup only
    na<-nrow(subset(obj,obj$AA == 0 & obj$AL > 0 & (obj$AH > 0 | obj$AC > 0 | obj$AP > 0 | obj$AN > 0)))/nrow(obj) #non-arenosa TS
    nl<-nrow(subset(obj,obj$AA > 0 & obj$AL == 0 & (obj$AH > 0 | obj$AC > 0 | obj$AP > 0 | obj$AN > 0)))/nrow(obj) #non-lyrata TS
    dn<-nrow(subset(obj,(rowSums(obj[,6:11])==0)))/nrow(obj)
    tp<-(1-ar-ly-na-nla-dn-nl)
    write.table(paste(p,tp,nl,na,nla,ar,ly,dn,sep = "\t"),"SVsource_proportionsSingletons_TSLyr.txt",append = T,row.names = F,col.names = F,quote = F)}
#NOSingletons
for (p in prot) { # p="AT4G18490"
  coord<-which(rowSums(his == p) > 0)
  obj <- his[coord[1]:coord[length(coord)],]
  ar<-nrow(subset(obj,obj$AA > 1 & obj$AL <2 & obj$AH <2 & obj$AC <2 & obj$AP <2 & obj$AN <2))/nrow(obj)
  ly<-nrow(subset(obj,obj$AL > 1 & obj$AA <2 & obj$AH <2 & obj$AC <2 & obj$AP <2 & obj$AN <2))/nrow(obj)
  na<-nrow(subset(obj,obj$AA < 2 & obj$AL > 1 & (obj$AH > 1 | obj$AC > 1 | obj$AP > 1 | obj$AN > 1)))/nrow(obj)
  nl<-nrow(subset(obj,obj$AA > 1 & obj$AL < 2 & (obj$AH > 1 | obj$AC > 1 | obj$AP > 1 | obj$AN > 1)))/nrow(obj)
  nla<-nrow(subset(obj,obj$AA <2 & obj$AL <2 & (obj$AH > 1 | obj$AC > 1 | obj$AP > 1 | obj$AN > 1)))/nrow(obj)
  dn<-nrow(subset(obj,obj$AL < 2 & obj$AA < 2 & obj$AH < 2 & obj$AC < 2 & obj$AP < 2 & obj$AN < 2))/nrow(obj)
  tp<-(1-ar-ly-nla-dn-na-nl)
  write.table(paste(p,tp,nl,na,nla,ar,ly,dn,sep = "\t"),"SVsource_proportionsNOSingletons_TSLyr.txt",append = T,row.names = F,col.names = F,quote = F)}
##Barplot
a1<-read.table("SVsource_proportionsSingletons_TSLyr.txt",h=F)
ord<-read.table("17WGDcand_names.txt",h=T)
a1<-a1[order(match(a1[,1],ord[,2])),]
colnames(a1)<-c("V1","1)","2)","3)","4)","5)","6)", "7)")
pdf("SVsource_proportionsSingletons_20JunSorted.pdf",width = 15,height = 5)
barplot(height = t(a1[,2:8]), col =  c("#ffa500ff", "#dc3f00ff", "#8c0000ff","#fed976ff","#ff6347ff","#ff0000ff","#bebebeff"), names.arg = a1$V1,legend.text = colnames(a1)[2:8],args.legend = list(x = "topleft", inset = c(0.45, -0.25)))
#dev.off()
a1<-read.table("SVsource_proportionsNOSingletons_TSLyr.txt",h=F)
ord<-read.table("17WGDcand_names.txt",h=T)
a1<-a1[order(match(a1[,1],ord[,2])),]
colnames(a1)<-c("V1","1)","2)","3)","4)","5)","6)", "7)") #"arenosa only","lyrata only","TS non arenosa","TS non lyrata","trans-sp.","outgroup only", "de-novo"
#pdf("SVsource_proportionsNOSingletons_TSLyr.pdf",width = 15,height = 5)
barplot(height = t(a1[,2:8]),  col=c("#ffa500ff", "#dc3f00ff", "#8c0000ff","#fed976ff","#ff6347ff","#ff0000ff","#bebebeff"), names.arg = a1$V1,legend.text = colnames(a1)[2:8],args.legend = list(x = "topleft", inset = c(0.45, -0.25)))
dev.off()



#### Genome-wide ############
his<-read.table("HIS0.7Syn17MayNAMES.intervals",h=T)
table(his$name)
sv<-read.table("ALL.SVPerSpecies.txt",h=T)
his<-cbind(his,sv)
#Incl.Singletons
  obj <- his
  ar<-nrow(subset(obj,obj$AA > 0 & obj$AL == 0 & obj$AH == 0 & obj$AC == 0 & obj$AP == 0 & obj$AN == 0))#arenosa only
  ly<-nrow(subset(obj,obj$AL > 0 & obj$AA == 0 & obj$AH == 0 & obj$AC == 0 & obj$AP == 0 & obj$AN == 0)) #lyrata only
  nla<-nrow(subset(obj,obj$AA == 0 & obj$AL == 0 & (obj$AH > 0 | obj$AC > 0 | obj$AP > 0 | obj$AN > 0))) #outgroup only
  na<-nrow(subset(obj,obj$AA == 0 & obj$AL > 0 & (obj$AH > 0 | obj$AC > 0 | obj$AP > 0 | obj$AN > 0))) #non-arenosa TS
  nl<-nrow(subset(obj,obj$AA > 0 & obj$AL == 0 & (obj$AH > 0 | obj$AC > 0 | obj$AP > 0 | obj$AN > 0))) #non-lyrata TS
  dn<-nrow(subset(obj,(rowSums(obj[,6:11])==0)))
  tp<-(nrow(his)-ar-ly-na-nla-dn-nl)
  write.table(t(as.data.frame(c("ar","ly","na","nl","tp","nla","dn"))),"ALLSVsource_proportionsSingletons_TSLyr.txt",append = F,row.names = F,col.names = F,quote = F)
    write.table(paste(ar,ly,na,nl,tp,nla,dn,sep = "\t"),"ALLSVsource_proportionsSingletons_TSLyr.txt",append = T,row.names = F,col.names = F,quote = F)
  ar<-nrow(subset(obj,obj$AA > 1 & obj$AL <2 & obj$AH <2 & obj$AC <2 & obj$AP <2 & obj$AN <2))
  ly<-nrow(subset(obj,obj$AL > 1 & obj$AA <2 & obj$AH <2 & obj$AC <2 & obj$AP <2 & obj$AN <2))
  na<-nrow(subset(obj,obj$AA < 2 & obj$AL > 1 & (obj$AH > 1 | obj$AC > 1 | obj$AP > 1 | obj$AN > 1)))
  nl<-nrow(subset(obj,obj$AA > 1 & obj$AL < 2 & (obj$AH > 1 | obj$AC > 1 | obj$AP > 1 | obj$AN > 1)))
  nla<-nrow(subset(obj,obj$AA <2 & obj$AL <2 & (obj$AH > 1 | obj$AC > 1 | obj$AP > 1 | obj$AN > 1)))
  dn<-nrow(subset(obj,obj$AL < 2 & obj$AA < 2 & obj$AH < 2 & obj$AC < 2 & obj$AP < 2 & obj$AN < 2))
  tp<-(nrow(obj)-ar-ly-nla-dn-na-nl)
  write.table(paste(ar,ly,na,nl,tp,nla,dn,sep = "\t"),"ALLSVsource_proportionsSingletons_TSLyr.txt",append = T,row.names = F,col.names = F,quote = F)















df<-reshape2::melt(data.frame( #melt to get each variable (i.e. A, B, C) in a single row
  a1 %>% #get rid of ID
    group_by(V1) %>% #group by category
    summarise_each(funs(sum))), #get the summation for each variable
  id.vars=c("V1")) 

  ggplot(aes(x=as.factor(V1),y=value,fill=variable),data = df)

  library(ggplot2)
  # create a dataset
  specie <- c(rep("sorgho" , 3) , rep("poacee" , 3) , rep("banana" , 3) , rep("triticum" , 3) )
  condition <- rep(c("normal" , "stress" , "Nitrogen") , 4)
  value <- abs(rnorm(12 , 0 , 15))
  data <- data.frame(specie,condition,value)
  
  # Stacked + percent
  ggplot(data, aes(fill=condition, y=value, x=specie)) + 
    geom_bar(position="fill", stat="identity")
  
  
+ geom_bar(stat = "identity",position="fill") +   scale_y_continuous(labels = scales::percent)+ scale_x_discrete() + scale_fill_manual(values=c("#FEA873FF","#F97C5DFF","#1D1147FF")) + labs(x = "",y = "",title = "Type of variants", fill = "source") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text.y = element_text(size=14))


(labels=c("Synonymous" = "Synonymous", "Nonsynonymous" = "Nonsynonymous", "Regulatory" = "Regulatory"))
  
  
  
############ 2. CLUSTERING RECOMBINATION RATE PER GENE
  setwd("/home/aa/lyrataWGD/fstScan/")
  library(data.table)
  library(dplyr)
  tot<-read.table("file:///home/aa/lyrataWGD/fstScan/haplotypeMaps/HIS0.7Syn17MayNAMES.intervals",h=T)
  tot$pos<-substr(tot$int,12,40)

    ############ 
    allDataMatrix<-readRDS("/home/aa/alpine/dmc/genomeScan/BALTISZEPSUB/lyrata_geneDensityData_binned.RDS")
    nchrom=8
    #plot quintile bins across chromosomes
    toPlotPoints = lapply(1 : nchrom, function(j) lapply(1 : 5, function(i) 
      allDataMatrix[allDataMatrix[, "quintile"] == i & allDataMatrix[, "scaffold"] == j, "position"]))
    #Start plot 
    png(paste("recombBins.png",sep=""),width = 1200,height = 2400,pointsize = 24)
    par(mfrow=c(8,1)) 
    for(i in 1 : nchrom) { # i=1
      if(length(toPlotPoints[[i]][[1]]) != 0) {
        plot(x = toPlotPoints[[i]][[1]], y = rep(0.5, length(toPlotPoints[[i]][[1]])), col = "red", main = paste("scaffold", i, sep = " "), yaxt = "n", xlab = "Position", ylab = "", xlim = c(0, max(unlist(toPlotPoints[[i]]))+ 1000))
        points(x = toPlotPoints[[i]][[2]], y = rep(0.5, length(toPlotPoints[[i]][[2]])), col = "orange")
        points(x = toPlotPoints[[i]][[3]], y = rep(0.5, length(toPlotPoints[[i]][[3]])), col = "yellow")
        points(x = toPlotPoints[[i]][[4]], y = rep(0.5, length(toPlotPoints[[i]][[4]])), col = "green")
        points(x = toPlotPoints[[i]][[5]], y = rep(0.5, length(toPlotPoints[[i]][[5]])), col = "blue")
      } else if (length(toPlotPoints[[i]][[2]]) != 0) {
        plot(x = toPlotPoints[[i]][[2]], y = rep(0.5, length(toPlotPoints[[i]][[2]])), col = "orange", main = paste("scaffold", i, sep = " "), yaxt = "n", xlab = "Position", ylab = "", xlim = c(0, max(unlist(toPlotPoints[[i]]))+ 1000))
        points(x = toPlotPoints[[i]][[3]], y = rep(0.5, length(toPlotPoints[[i]][[3]])), col = "yellow")
        points(x = toPlotPoints[[i]][[4]], y = rep(0.5, length(toPlotPoints[[i]][[4]])), col = "green")
        points(x = toPlotPoints[[i]][[5]], y = rep(0.5, length(toPlotPoints[[i]][[5]])), col = "blue")
      } else  {
        plot(x = toPlotPoints[[i]][[3]], y = rep(0.5, length(toPlotPoints[[i]][[3]])), col = "orange", main = paste("scaffold", i, sep = " "), yaxt = "n", xlab = "Position", ylab = "", xlim = c(0,max(unlist(toPlotPoints[[i]]))+ 1000))
        points(x = toPlotPoints[[i]][[4]], y = rep(0.5, length(toPlotPoints[[i]][[4]])), col = "green")
        points(x = toPlotPoints[[i]][[5]], y = rep(0.5, length(toPlotPoints[[i]][[5]])), col = "blue")
      }
      
      gs1<-subset(tot,substr(tot$int,10,10) %in% i)
      gs1$mean<-gs1$pos
      if (nrow(gs1)>0) {
        for (i in 1:nrow(gs1)) {
          
          points(x = gs1[i,4],y= 0.5,col="black",pch="|", cex=1.5)
        }
      }
    }
    dev.off()
  
  
  
  
  
