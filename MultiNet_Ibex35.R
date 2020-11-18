###################################################################
## MultiNet algorithm - Implementation with IBEX35 Data
## Part III - Chapter 8 - 8.6
## Sara Prada, October 2020
###################################################################

#Needed packages

library(tidyquant)
library(rtsdata)
library(quantmod)
library(BatchGetSymbols)
library(BiocGenerics)
library(BBmisc)
library(gtools)

#Call supplemental scripts containing needed functions

source("F:/DOCTORADO/PROGRAMAS/TESIS PROGRAMS/MultiNet_Supplemental1.R") 
source("F:/DOCTORADO/PROGRAMAS/TESIS PROGRAMS/MultiNet_Supplemental2.R") 

#Directory

setwd("F:/DOCTORADO/TESIS/ibex/test")

#Download IBEX35 data

tickers = c("SAN.MC","ENG.MC","GRF.MC","CLNX.MC","VIS.MC","REE.MC","IBE.MC","TL5.MC","ELE.MC","NTGY.MC","BKT.MC","SGRE.MC",
            "CABK.MC", "MAP.MC","ANA.MC","FER.MC","TEF.MC","ITX.MC","AENA.MC", "MTS.MC","ENC.MC","BBVA.MC","COL.MC","AMS.MC",
            "ACX.MC","ACS.MC","MRL.MC","SAB.MC","IAG.MC")

G<-getSymbols(tickers,
           from = "2005-01-01",
           to = "2020-10-01")

chart_Series(BBVA.MC)

start_date <- "2005-01-01"
end_date <- "2020-10-01"
#end_date <- Sys.Date()-1

num_start_date<-as.numeric(as.Date(start_date))+1
num_end_date<-as.numeric(as.Date(end_date))
daysspent<-num_end_date-num_start_date+1

master_df <-matrix(NA,length(tickers),daysspent)
colnames(master_df)<-seq(num_start_date,num_end_date,1)
rownames(master_df)<-tickers

aux<-seq(1,daysspent,1)
aux1<-seq(num_start_date,num_end_date,1)
freq.data <- 'daily'

for (idx in seq(length(tickers))){
  stock_index = tickers[idx]
  l.out <- BatchGetSymbols(tickers = stock_index, 
                           first.date = start_date,
                           last.date = end_date, 
                           freq.data = freq.data) # cache in tempdir()
  
  #getSymbols(stock_index, src = "yahoo", from=start_date,to=end_date)
  temp_df = as.data.frame(l.out$df.tickers)[,c(1,2,3,4,5,6,7)]
  colnames(temp_df) = c("Open", "High", "Low", "Close", 
                        "Volume", "Adjusted","Date")
  
  temp_df$Change<-(temp_df$Close-temp_df$Open)
  
  temp_dfdat<-temp_df$Change
  names(temp_dfdat)<-as.numeric(as.Date((temp_df$Date)))

  for (j in aux){
  for (dat in aux1[j]){
     master_df[idx,j]<-temp_dfdat[as.character(get("dat"))]
  }}
}

head(master_df)
dim(master_df) #nrow is sample size, ncol is variables

#Remove missing observations because they are weekends

master_df1<-master_df[,!apply(is.na(master_df), 2, all)]
head(master_df1)
dim(master_df1)
#rownames(master_df1)<-rownames(master_df)
#colnames(master_df1)<-colnames(master_df)

master_df11<-t(master_df1)
#master_df11<-master_df11[,colSums((master_df11))==0]
master_df11[is.na(master_df11)]<-0
dim(master_df11)
head(master_df11)

for (i in 1:dim(master_df11)[2]){
  #master_df11[,i]<-abs(master_df11[,i])
  master_df11[,i]<-normalize(master_df11[,i])
}

#Select non null companies

master_df11<-master_df11[,c(1,2,5,6,7,8,9,11,12,15,16,17,18,21,22,23,25,26,28)]

#Prepare data for MultiNet

median_df<-apply(master_df11,1,median)
names(median_df)<-rownames(master_df11)
var_df<-apply(master_df11,1,var)
names(var_df)<-rownames(master_df11)

all.dat_ord_global<-na.omit(master_df11)
all.dat_ord_global[,1]
dim(all.dat_ord_global)
head(all.dat_ord_global)

median_df<-sort(apply(all.dat_ord_global,1,median),decreasing=TRUE)
names(median_df)<-rownames(all.dat_ord_global)
var_df<-apply(all.dat_ord_global,1,var)
names(var_df)<-rownames(all.dat_ord_global)

all.dat_ord_global<-na.omit(all.dat_ord_global)
dim(all.dat_ord_global)

#Run the algorithm - without case/control groups, only overall
init=1
casenet=0
controlnet=0
biological=0
maxiall=dim(all.dat_ord_global)[1]-by
by=60
max=maxiall
ov=30
int=2
ov_i=50
opt_cluster=2
pero=0.2
corm=0.8
method_map="Forgy"

multinet(all.dat_ord_global=all.dat_ord_global ,medi_global= median_df,var_s1=var_df,diff_global1=median_df)

#Post-processing

exp2<-exp1_all_global[exp1_all_global %in% c("1")]

if (length(exp2)>0){
  
exp21<-as.numeric(names(exp2))
  
points<-list()
points<-pointnode_all_global[exp21]
length(points)
  
n<-rownames(all.dat_ord_global)
  
n2<-list()
for (i in 1:length(points)){
  n2[[i]]<-n[points[[i]]]
}
length(n2)
  
n3<-as.numeric()
n3<-unique(unlist(n2))
length(n3)
}

dms_global<-n3

h1<-Heatmap(as.matrix(t(all.dat_ord_global[dms_global,colnames(all.dat_ord_global)])), name = "Price Evolution",column_order=as.character(sort(as.numeric(n3))),
            show_row_names=TRUE,show_column_names=FALSE,show_column_dend = FALSE,show_row_dend=FALSE,row_dend_reorder = TRUE,column_title = "",column_km = 3, 
            row_km=3,show_parent_dend_line = TRUE) 
length(n3)
print(h1)
h1<-draw(h1)
  
aux_dat<-all.dat_ord_global[n3,]
  
r.dend <- column_order(h1)
#c.dend <- column_dend(h1)
#c.dend$`1`
aux=r.dend[[1]]
  
test_dat<-aux_dat[aux,]
sort(as.Date(as.numeric(rownames(test_dat))))
range(sort(as.Date(as.numeric(rownames(test_dat)))))
  
dim(test_dat)
rownames(test_dat)<-as.Date(as.numeric(rownames(test_dat)))
cori<-cor(t(test_dat))
cori[abs(cori)<0.5]<-0
corrplot(cori,tl.pos='n')
corrplot(cor(t(test_dat[as.character(sort(as.numeric(rownames(test_dat)))),])))
corh<-Heatmap(cori,row_order=rownames(cori),column_order=colnames(cori))
sort(as.Date(as.numeric(rownames(cori))))
plot(apply(test_dat,1,mean),type="l")
  
test_dat<-all.dat_ord_global[aux,]
cori<-cor(t(test_dat))
cori[abs(cori)<0.8]<-0
corrplot(cori)
  
  
G<-getSymbols(tickers,
                from = "2010-01-08",
                to = "2010-07-14")
  
chart_Series(IBEX.MC)

ggscatter(as.data.frame(test_dat), x = "SAN.MC", y = "BBVA.MC",
            add = "reg.line", conf.int = TRUE,
            cor.coef = TRUE, cor.method = "pearson",
            xlab = "Santander", ylab = "BBVA")
  

  