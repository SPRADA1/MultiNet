###################################################################
## MultiNet algorithm 
## Part III - Chapter 8 
## Sara Prada, October 2020
###################################################################

library(minfi)
library(GEOquery)
library(igraph)
library(png)
library(grid)
library(gridExtra)
library(gapminder)
library(tidyverse)
library(lme4)
library(ggpubr)
library(np)
library(logistf)
library(sgof)
library(qvalue)
library(multtest)
library(randomForest)
library(circlize)
library(BioCircos)
library(circular)
library(OmicCircos)
library(migest)
library(ReactomePA)
library(enrichplot)
library(clusterProfiler)
library(org.Hs.eg.db)
library(qqman)
library(calibrate)
library(netbiov)
library(plot3D)
library(ggplot2)
library(gtools)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylation450kmanifest)
library(effsize)
library(scales)
library(ComplexHeatmap)
library(DOSE)
library(corrplot)
library(wesanderson)
library(RColorBrewer)

#Directory
setwd("F:/DOCTORADO/TESIS/GSE76938/test/5")

#Call supplemental scripts containing needed functions
source("F:/DOCTORADO/PROGRAMAS/TESIS PROGRAMS/MultiNet_Supplemental1.R") 
source("F:/DOCTORADO/PROGRAMAS/TESIS PROGRAMS/MultiNet_Supplemental2.R") 

#################################
#Download the database/450K array
gse1<-getGEO("GSE76938") 
gse1

#Convert the database into a data frame
all.dat0<-exprs(gse1[[1]])
dim(all.dat0)
head(all.dat0)

all.dat0<-as.data.frame(exprs(gse_h3))
head(all.dat0)
dim(all.dat0)

#Subjects metadata
all.meta1 <- pData(phenoData(gse1[[1]]))
dim(all.meta1)
head(all.meta1)

all.meta1 <- pData(phenoData(gse_h3))
head(all.meta1)

#Annotation database
ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
dim(ann450k)
head(ann450k)
annotation0<-ann450k

#Define sample groups
names(all.meta1)[34]<-"disease"
case <- all.meta1[all.meta1$disease %in% c("prostate cancer tissue"),]
dim(case)
rcase<-rownames(case)
control <- all.meta1[all.meta1$disease %in% c("prostate benign tissue"),]
dim(control) 
rcontrol<-rownames(control)

#Remove SNPs and sites related to sex chromosomes
newdata <- intersect(rownames(all.dat0), rownames(annotation0))
all.dat<-all.dat0[newdata,]
dim(all.dat)

annotation<-annotation0[rownames(all.dat),]
dim(annotation)

##Eliminamos cromosomas X e Y
nox<-annotation[annotation$chr!="chrX" & annotation$chr!="chrY",]
noxr <- rownames(nox)

##Eliminamos SNPs
snps<-annotation[is.na(annotation$Probe_rs),]
snpsr <- rownames(snps)

bothr<-intersect(noxr,snpsr)
length(bothr)

all.dat1 <- na.omit(all.dat[bothr,c(rcase,rcontrol)])
dim(all.dat1)
head(all.dat1)

#sam<-sample(40,20)
#rcase<-rownames(all.meta1)[sam]
#rcontrol<-rownames(all.meta1)[-sam]

#Filter functions
medi_global<-apply(all.dat1,1,median)
medi_global_case<-apply(all.dat1[,rcase],1,median)
medi_global_control<-apply(all.dat1[,rcontrol],1,median)

var <- apply(na.omit(all.dat1),1,var)
var_s1 <- sort(var,decreasing=TRUE)
var_case <- na.omit(apply(all.dat1[,rcase],1,var))
var_s1_case <- sort(var_case,decreasing=TRUE)
var_control <- na.omit(apply(all.dat1[,rcontrol],1,var))
var_s1_control <- sort(var_control ,decreasing=TRUE)

diff_global<-na.omit(abs(medi_global_case-medi_global_control))
diff_global1 <- sort(diff_global,decreasing=TRUE)
length(diff_global1)

##Order by maximum difference of medians
all.dat_ord_global<-na.omit(all.dat1[names(diff_global1),c(rcase,rcontrol)]) 
dim(all.dat_ord_global)

all.meta1<-all.meta1[c(rcase,rcontrol),]
dim(all.meta1)

###############
#Visualize data
#3D plot
selection_dat<-diff_global1[1:5000]
cor_global<-cor(t(all.dat_ord_global[names(selection_dat),]))
cor_mean<-apply(cor_global,1,mean)

data <- as.data.frame(cbind(var[names(selection_dat)], medi_global[names(selection_dat)],diff_global1[names(selection_dat)]))
head(data)
colnames(data)<-c("var","medi","diff")

scatter3D(data$medi,data$var, data$diff, phi = 0, bty ="g",xlab = "Median", ylab = "Variance", zlab = "Difference",ticktype = "detailed")
scatter3D(data$medi,data$var, data$diff,phi = 0, bty = "g",  type = "h", 
          ticktype = "detailed", pch = 19, cex = 0.5,xlab = "Median", ylab = "Variance", zlab = "Difference")
#scatter3d(x =data$medi, y = data$var, z =data$diff, point.col = "blue", surface=FALSE)

##By subject
all.dat_ord_global_s<-all.dat_ord_global[,c(rcase,rcontrol)]
all.meta1<-all.meta1[c(rcase,rcontrol),]

samples<-c(rep("Case",length(rcase)),rep("Control",length(rcontrol)))
names(samples)<-rownames(all.meta1[colnames(all.dat_ord_global_s),])

median_s<-apply(all.dat_ord_global_s,2,median)
variance_s<-apply(all.dat_ord_global_s,2,var)

ploti<-as.data.frame(cbind(median_s, variance_s))
head(ploti)
ploti$pop<-samples

plot12<-ggplot(ploti, aes(x=median_s, y=variance_s, color=pop, shape=pop, size=pop)) + 
  geom_point(size=4) +
  xlab("Methylation median") + ylab("Methylation variance")+
  theme(legend.direction = 'horizontal', legend.position = 'top',legend.title  =element_blank(),text = element_text(size=15))
png("Subject_Distribution_all.png", width = 10, height = 10, units = 'in', res = 600)
print(plot12)
dev.off()

##By CpG
df<-data.frame(medi_global, diff_global1[names(medi_global)])
colnames(df)<-c("Median","Median_Difference")
plot2<-ggplot(df,aes(x=Median,y=Median_Difference) ) +
  xlab("Methylation median") + ylab("Methylation difference")+
  geom_point(size=2,aes(color = Median)) + scale_color_gradient(low = "blue", high = "red")
png("Filter_Global.png", width = 10, height = 10, units = 'in', res = 600)
print(plot2)
dev.off()

##Case
df<-data.frame(medi_global_case,var_s1_case)
colnames(df)<-c("Median","Variance_case")
plot3<-ggplot(df,aes(x=Median,y=Variance_case) ) +
  xlab("Methylation median") + ylab("Methylation variance")+
  geom_point(size=2, aes(color = Median)) + scale_color_gradient(low = "blue", high = "red")
png("Filter_Case.png", width = 10, height = 10, units = 'in', res = 600)
print(plot3)
dev.off()

##Control
df<-data.frame(medi_global_control, var_s1_control)
colnames(df)<-c("Median","Variance_control")
plot4<-ggplot(df,aes(x=Median,y=Variance_control) ) +
  xlab("Methylation median") + ylab("Methylation variance")+
  geom_point(size=2,aes(color = Median)) + scale_color_gradient(low = "blue", high = "red")
png("Filter_Control.png", width = 10, height = 10, units = 'in', res = 600)
print(plot4)
dev.off()

##################
#Run the algorithm
multinet(all.dat_ord_global=all.dat_ord_global ,filter1_global= medi_global,filter2_global=diff_global1,
         filter1_control=medi_global_control,filter2_control=var_s1_control,
         filter1_case=medi_global_case,filter2_case=var_s1_case)


###############
#Local MultiNet

chr1<-as.numeric(substr(annotation[rownames(all.dat1),]$chr,4,6))
sort(chr1)
annotation_chr<-cbind(annotation[rownames(all.dat1),],chr1)

# 1, 4, 5, 7, 8, 11, 16 and 19 
chromosome=22
annotation_dat<-annotation_chr[rownames(all.dat1),]
annot_chr<-rownames(annotation_dat[annotation_dat$chr1==chromosome,])
length(annot_chr)

all.dat_ord_local1<-na.omit(all.dat1[annot_chr,])

medi_local_case<-apply(all.dat_ord_local1[,rcase],1,median)
medi_local_control<-apply(all.dat_ord_local1[,rcontrol],1,median)
diff_local<-na.omit(abs(medi_local_case-medi_local_control))
diff_local1 <- sort(diff_local,decreasing=TRUE)
length(diff_local1)

medi_local<-apply(all.dat_ord_local1,1,median)
medi_local<-sort(medi_local)

all.dat_ord_local1<-all.dat_ord_local1[names(diff_local1),]

site_local<-annotation_dat[rownames(all.dat_ord_local1),2]
names(site_local)<-rownames(all.dat_ord_local1)

var_local <- apply(na.omit(all.dat_ord_local1),1,var)
var_s1_local <- sort(var_local,decreasing=TRUE)

all.dat_ord_local1<-all.dat_ord_local1[names(site_local),]
dim(all.dat_ord_local1)

#all.meta1$disease="Children"

data <- as.data.frame(cbind(site_local[names(diff_local1)], medi_local[names(diff_local1)],diff_local1[names(diff_local1)]))
head(data)
colnames(data)<-c("var","medi","diff")

scatter3D(data$medi,data$var, data$diff, phi = 0, bty ="g",xlab = "Median", ylab = "Position", zlab = "Difference",ticktype = "detailed")

scatter3D(data$medi,data$var, data$diff,phi = 0, bty = "g",  type = "h", 
          ticktype = "detailed", pch = 19, cex = 0.5,xlab = "Median", ylab = "Variance", zlab = "Difference")
#scatter3d(x =data$medi, y = data$var, z =data$diff, point.col = "blue", surface=FALSE)

multinet(all.dat_ord_global=all.dat_ord_local1 ,medi_global= medi_local,var_s1=site_local, diff_global1=medi_local,
         medi_global_control=site_local,var_s1_control=medi_local_case,
         medi_global_case=site_local,var_s1_case=medi_local_control)
