###################################################################
## MultiNet algorithm - Suppelemental 1 - Main algorithm
## Part III - Chapter 8 
## Sara Prada, October 2020
###################################################################

multinet<-defmacro(all.dat_ord_global ,filter1_global,filter2_global,
                   filter1_control,filter2_control,
                   filter1_case,filter2_case, expr = {
  
###################################################################
#READ PARAMETERS
                     
print("Please include MultiNet parameters")

biological<-as.integer( readline(prompt="Biological results: "))                   
                     
casenet<- as.integer( readline(prompt="Case option: "))
if (is.na(casenet)) {casenet=1}
  
controlnet<- as.integer( readline(prompt="Control option: "))
if (is.na(controlnet)) {controlnet=1}
  
init <- as.integer( readline(prompt="Initial position: "))
if (is.na(init)) {init=1}
  
maxiall<-as.integer(readline(prompt="Last position: "))
if (is.na(maxiall)) {maxiall=5000}
  
by <- as.integer(readline(prompt="Window size: "))
if (is.na(by)) {by=1000}
  
ov0<-as.integer(readline(prompt="Percent window overlap: "))
if (is.na(ov0)) {ov0=50}
  
int<-as.integer(readline(prompt="Number of intervals: "))
if (is.na(int)) {int=2}
  
ov_i<-as.numeric(readline(prompt="Percent interval overlap: "))
if (is.na(ov_i)) {ov_i=50}
  
opt_cluster<-as.integer(readline(prompt="Number of clusters: "))
if (is.na(opt_cluster)) {opt_cluster=3}
  
pero<-as.numeric(readline(prompt="Points in common among nodes: "))
if (is.na(pero)) {pero=0.5}
  
corm<-as.numeric(readline(prompt="Median correlation among nodes: "))
if (is.na(corm)) {corm=0.8}
  
method_map<-readline(prompt="Cluster method: ")
if (method_map=="") {method_map="Forgy"}
  
###################################################################
##INITIALIZE PARAMETERS

graphic_global<-list()
graphic_local<-list()
graphic_control_global<-list()
graphic_control_local<-list()
graphic_case_global<-list()
graphic_case_local<-list()
mapp_global<-list()
mapp_local<-list()
mapp_case_global<-list()
mapp_case_local<-list()
mapp_control_global<-list()
mapp_control_local<-list()
links_global<-list()
links_local<-list()
nodes_global<-list()
nodes_local<-list()
links_case_global<-list()
links_case_local<-list()
nodes_case_global<-list()
nodes_case_local<-list()
links_control_global<-list()
links_control_local<-list()
nodes_control_global<-list()
nodes_control_local<-list()

max=maxiall
ov=(ov0/100)*by
#sleep_for_a_minute <- function() { Sys.sleep(60) }
start_time <- Sys.time()
#sleep_for_a_minute()

#Execution first window
print("MultiNet per window - start")

tictoc::tic("First window")
corr(init,init+by,filter1_global,filter2_global,
     filter1_control,filter2_control,
     filter1_case,filter2_case)
tictoc::toc()

#The rest of the windows
maxi0=max
a<-seq(init+by-ov,maxi0,ov) 
for (h in a){
  j=h+by
  tictoc::tic("Next window")
  corr(h,j,filter1_global,filter2_global,
       filter1_control,filter2_control,
       filter1_case,filter2_case) 
  tictoc::toc()
}
maxi=h+by
tictoc::toc()
tictoc::toc()
tictoc::toc()

print("MultiNet per window - end")

###################################################################
##Macro to reorganize vertices between nodes with points in common

print("MultiNet: joining of window networks - start")

corr1(mapp_global,graphic_global,links_global,nodes_global,all.dat_ord_global)
graphic1_global<-graphic1

if (controlnet==1){
  corr1(mapp_control_global,graphic_control_global,links_control_global,nodes_control_global, all.dat_ord_global)
  graphic1_control_global<-graphic1
}

if (casenet==1){
  corr1(mapp_case_global,graphic_case_global,links_case_global,nodes_case_global, all.dat_ord_global)
  graphic1_case_global<-graphic1
}

print("MultiNet: joining of window networks - end")

end_time <- Sys.time()

total_time<-end_time - start_time

print("MultiNet running time")
print(total_time)

##########################################################
##Graph representation

reset_par()

#global
deg_global<-igraph::degree(graphic1_global[[maxi]])

#adjacency
adja_global<-get.adjacency(graphic1_global[[maxi]])

graphic2_global<-induced_subgraph(graphic1_global[[maxi]], vids=V(graphic1_global[[maxi]]), impl = c("auto", "copy_and_delete","create_from_scratch"))

png("Global_Overall_Network.png", width = 6, height = 6, units = 'in', res = 600)
mst.plot.mod(graphic2_global,bg="black", v.size=3,e.size=.25,
             mst.e.size=1.2, sf=0,mst.edge.col="white",colors= "white", vertex.color = "green",
             v.lab=FALSE,v.lab.cex=1,v.lab.col="white")
text(0.5,1,"Global Overall Network",col="white", cex=0.9)
text(0.5,0.95,"Filter1: Median",col="white", cex=0.9)
text(0.5,0.9,"Filter2: Median difference case-control",col="white", cex=0.9)
text(0.5,0.85,paste(maxiall, "CpGs by ", by, "with", ov, "overlap", sep=" "),col="white", cex=0.9)
text(0.5,0.80,paste(int, "intervals with 50% overlap", sep=" "),col="white",cex=0.9)
text(0.5,0.75,paste(opt_cluster, "clusters", sep=" "),col="white",cex=0.9)
dev.off()

if (casenet==1){
  
  deg_global_case<-igraph::degree(graphic1_case_global[[maxi]])
  adja_global_case<-get.adjacency(graphic1_case_global[[maxi]])
  graphic2_global_case<-induced_subgraph(graphic1_case_global[[maxi]], vids=V(graphic1_case_global[[maxi]]), impl = c("auto", "copy_and_delete","create_from_scratch"))
  
  png("Global_Case_Network.png", width = 6, height = 6, units = 'in', res = 600)
  mst.plot.mod(graphic2_global_case,bg="black", v.size=5,e.size=.25,
               mst.e.size=1.2, sf=0,mst.edge.col="white",colors= "white", vertex.color = "green",
               v.lab=TRUE,v.lab.cex=0.7,v.lab.col="red")
  text(0.5,1,"Global Case Network",col="white", cex=0.9)
  text(0.5,0.95,"Filter1: Median",col="white", cex=0.9)
  text(0.5,0.9,"Filter2: Variance",col="white", cex=0.9)
  text(0.5,0.85,paste(maxiall, "CpGs by 1000 with 50% overlap", sep=" "),col="white", cex=0.9)
  text(0.5,0.80,paste(int, "intervals with 50% overlap", sep=" "),col="white", cex=0.9)
  text(0.5,0.75,paste(opt_cluster, "clusters", sep=" "),col="white", cex=0.9)
  dev.off()
}

if (controlnet==1){
  
  deg_global_control<-igraph::degree(graphic1_control_global[[maxi]])
  adja_global_control<-get.adjacency(graphic1_control_global[[maxi]])
  graphic2_global_control<-induced_subgraph(graphic1_control_global[[maxi]], vids=V(graphic1_control_global[[maxi]]), impl = c("auto", "copy_and_delete","create_from_scratch"))
  
  png("Global_Control_Network.png", width = 6, height = 6, units = 'in', res = 600)
  mst.plot.mod(graphic2_global_control,bg="black", v.size=5,e.size=.25,
               mst.e.size=1.2, sf=0,mst.edge.col="white",colors= "white", vertex.color = "green",
               v.lab=TRUE,v.lab.cex=0.7,v.lab.col="red")
  text(0.5,1,"Global Control Network",col="white", cex=0.9)
  text(0.5,0.95,"Filter1: Median",col="white", cex=0.9)
  text(0.5,0.9,"Filter2: Variance",col="white", cex=0.9)
  text(0.5,0.85,paste(maxiall, "CpGs by 1000 with 50% overlap", sep=" "),col="white", cex=0.9)
  text(0.5,0.80,paste(int, "intervals with 50% overlap", sep=" "),col="white", cex=0.9)
  text(0.5,0.75,paste(opt_cluster, "clusters", sep=" "),col="white", cex=0.9)
  dev.off()
}

###################################################################
#NETWORK MANIPULATION

#Macro to study the CpGs contained in the selected differentiated nodes  
print("MultiNet: colored networks - start")

difnode (mapp_global, graphic2_global,all.dat_ord_global,subj=colnames(all.dat_ord_global))
exp1_all_global=exp1
dif1_all_global=dif1
var_all_global=var1
cor1_med_all_global=cor1_med
pointnode_all_global=p_all
n2_all_global=n2_all

r1<-paste("<q2",round(quantile(exp)[2],2),sep=" ")
r2<-paste("q2-q4",round(quantile(exp)[2],2),round(quantile(exp)[4],2),sep=" ")
r3<-paste(">q4",round(quantile(exp)[4],2),sep=" ")

plotmap(graphic2_global,exp1_all_global,cat=c(r3,r2,r1),
        text1="Global_median_methylation.png",text="Overall by median methylation")

if (casenet==1 & controlnet==1){
r1<-paste("<q2",round(quantile(nod_dif)[2],2),sep=" ")
r2<-paste("q2-q4",round(quantile(nod_dif)[2],2),round(quantile(nod_dif)[4],2),sep=" ")
r3<-paste(">q4",round(quantile(nod_dif)[4],2),sep=" ")

plotmap(graphic2_global,dif1_all_global,cat=c(r3,r2,r1),
        text1="Global_median_difference_methylation.png",text="Overall by diff methylation")
}

r1<-paste("<q2",round(quantile(nod_var)[2],2),sep=" ")
r2<-paste("q2-q4",round(quantile(nod_var)[2],2),round(quantile(nod_var)[4],2),sep=" ")
r3<-paste(">q4",round(quantile(nod_var)[4],2),sep=" ")

plotmap(graphic2_global,var_all_global,cat=c(r3,r2,r1),
        text1="Global_median_variance.png",text="Overall by variance")

r1<-paste("<q2",round(quantile(cor_met_med)[2],2),sep=" ")
r2<-paste("q2-q4",round(quantile(cor_met_med)[2],2),round(quantile(cor_met_med)[4],2),sep=" ")
r3<-paste(">q4",round(quantile(cor_met_med)[4],2),sep=" ")

plotmap(graphic2_global,cor1_med_all_global,cat=c(r3,r2,r1),
        text1="Global_median_correlation.png",text="Overall by median correlation")

#case
if (casenet==1){
  difnode (mapp_case_global, graphic2_global_case,all.dat_ord_global,subj=rcase)#case
  exp1_case_global=exp1
  dif1_case_global=dif1
  var_case_global=var1
  cor1_med_case_global=cor1_med
  pointnode_case_global=p_all
  n2_case_global=n2_all
  
  exp_case<-exp
  nod_dif_case<-nod_dif
  cor_met_med_case<-cor_met_med
  nod_var_case<-nod_var
  
  r1<-paste("<q2",round(quantile(exp_case)[2],2),sep=" ")
  r2<-paste("q2-q4",round(quantile(exp_case)[2],2),round(quantile(exp_case)[4],2),sep=" ")
  r3<-paste(">q4",round(quantile(exp_case)[4],2),sep=" ")
  
  plotmap(graphic2_global_case,exp1_case_global,cat=c(r3,r2,r1),
          text1="Global_Case_median_methylation.png",text="Case by median methylation")
  
  r1<-paste("<q2",round(quantile(nod_dif_case)[2],2),sep=" ")
  r2<-paste("q2-q4",round(quantile(nod_dif_case)[2],2),round(quantile(nod_dif_case)[4],2),sep=" ")
  r3<-paste(">q4",round(quantile(nod_dif_case)[4],2),sep=" ")
  
  plotmap(graphic2_global_case,dif1_case_global,cat=c(r3,r2,r1),
          text1="Global_Case_median_difference_methylation.png",text="Case by diff methylation")
  
  r1<-paste("<q2",round(quantile(nod_var_case)[2],2),sep=" ")
  r2<-paste("q2-q4",round(quantile(nod_var_case)[2],2),round(quantile(nod_var_case)[4],2),sep=" ")
  r3<-paste(">q4",round(quantile(nod_var_case)[4],2),sep=" ")
  
  plotmap(graphic2_global_case,var_case_global,cat=c(r3,r2,r1),
          text1="Global_case_median_variance.png",text="Case by variance")
  
  r1<-paste("<q2",round(quantile(cor_met_med_case)[2],2),sep=" ")
  r2<-paste("q2-q4",round(quantile(cor_met_med_case)[2],2),round(quantile(cor_met_med_case)[4],2),sep=" ")
  r3<-paste(">q4",round(quantile(cor_met_med_case)[4],2),sep=" ")
  
  plotmap(graphic2_global_case,cor1_med_case_global,cat=c(r3,r2,r1),
          text1="Global_Case_median_correlation.png",text="Case by median correlation")
}

#control
if (controlnet==1){
  difnode (mapp_control_global, graphic2_global_control,all.dat_ord_global,subj=rcontrol)#control
  exp1_control_global=exp1
  dif1_control_global=dif1
  var_control_global=var1
  cor1_med_control_global=cor1_med
  pointnode_control_global=p_all
  n2_control_global=n2_all
  
  exp_control<-exp
  nod_dif_control<-nod_dif
  nod_var_control<-nod_var
  cor_met_med_control<-cor_met_med
  
  r1<-paste("<q2",round(quantile(exp_control)[2],2),sep=" ")
  r2<-paste("q2-q4",round(quantile(exp_control)[2],2),round(quantile(exp_control)[4],2),sep=" ")
  r3<-paste(">q4",round(quantile(exp_control)[4],2),sep=" ")
  
  plotmap(graphic2_global_control,exp1_control_global,cat=c(r3,r2,r1),
          text1="Global_Control_median_methylation.png",text="Control by median methylation")
  
  # #case reference
  # r1<-paste("<q2",round(quantile(exp_case)[2],2),sep=" ")
  # r2<-paste("q2-q4",round(quantile(exp_case)[2],2),round(quantile(exp_case)[4],2),sep=" ")
  # r3<-paste(">q4",round(quantile(exp_case)[4],2),sep=" ")
  # 
  # plotmap(graphic2_global_control,exp1_control_global,cat=c(r3,r2,r1),
  #         text1="Global_Control_median_methylation_oth.png",text="Older by median methylation")
  # 
  r1<-paste("<q2",round(quantile(nod_dif_control)[2],2),sep=" ")
  r2<-paste("q2-q4",round(quantile(nod_dif_control)[2],2),round(quantile(nod_dif_control)[4],2),sep=" ")
  r3<-paste(">q4",round(quantile(nod_dif_control)[4],2),sep=" ")
  
  plotmap(graphic2_global_control,dif1_control_global,cat=c(r3,r2,r1),
          text1="Global_Control_median_difference_methylation.png",text="Control by diff methylation")
  
  r1<-paste("<q2",round(quantile(nod_var_control)[2],2),sep=" ")
  r2<-paste("q2-q4",round(quantile(nod_var_control)[2],2),round(quantile(nod_var_control)[4],2),sep=" ")
  r3<-paste(">q4",round(quantile(nod_var_control)[4],2),sep=" ")
  
  plotmap(graphic2_global_control,var_control_global,cat=c(r3,r2,r1),
          text1="Global_control_median_variance.png",text="Control by variance")
  
  r1<-paste("<q2",round(quantile(cor_met_med_control)[2],2),sep=" ")
  r2<-paste("q2-q4",round(quantile(cor_met_med_control)[2],2),round(quantile(cor_met_med_control)[4],2),sep=" ")
  r3<-paste(">q4",round(quantile(cor_met_med_control)[4],2),sep=" ")
  
  plotmap(graphic2_global_control,cor1_med_control_global,cat=c(r3,r2,r1),
          text1="Global_Control_median_correlation.png",text="Control by median correlation")
}

print("MultiNet: colored networks - end")

#Histograms of the three variables distribution among case and control

if (casenet==1 & controlnet==1){
  
print("MultiNet: node histograms - start")
  
df <- data.frame(
  sample=factor(c(rep("Case",length(exp_case)),rep("Control",length(exp_control)))),
  nodes=(c(exp_case,exp_control)))
head(df)

hist1<-ggplot(df, aes(x=nodes, color=sample, fill=sample)) + 
  geom_histogram(aes(y=..density..), alpha=0.5, position="identity")+
  geom_density(alpha=.2) +xlab("Node Median methylation")+
  theme(text = element_text(size=20))
png("hist_median.png", width = 6, height =6, units = 'in', res = 600)
print(hist1)
dev.off()

#difference
df <- data.frame(
  sample=factor(c(rep("Case",length(nod_dif_case)),rep("Control",length(nod_dif_control)))),
  nodes=(c(nod_dif_case,nod_dif_control)))
head(df)

hist2<-ggplot(df, aes(x=nodes, color=sample, fill=sample)) + 
  geom_histogram(aes(y=..density..), alpha=0.5, position="identity")+
  geom_density(alpha=.2) +xlab("Node Median difference methylation")+
  theme(text = element_text(size=20))
png("hist_mediandiff.png", width = 6, height =6, units = 'in', res = 600)
print(hist2)
dev.off()

#correlation
df <- data.frame(
  sample=factor(c(rep("Case",length(cor_met_med_case)),rep("Control",length(cor_met_med_control)))),
  nodes=(c(cor_met_med_case,cor_met_med_control)))
head(df)

hist3<-ggplot(df, aes(x=nodes, color=sample, fill=sample)) + 
  geom_histogram(aes(y=..density..), alpha=0.5, position="identity")+
  geom_density(alpha=.2) +xlab("Node Median correlation")+
  theme(text = element_text(size=20))
png("hist_mediancorr.png", width = 6, height =6, units = 'in', res = 600)
print(hist3)
dev.off()

#variance
df <- data.frame(
  sample=factor(c(rep("Case",length(nod_var_case)),rep("Control",length(nod_var_control)))),
  nodes=(c(nod_var_case,nod_var_control)))
head(df)

hist4<-ggplot(df, aes(x=nodes, color=sample, fill=sample)) + 
  geom_histogram(aes(y=..density..), alpha=0.5, position="identity")+
  geom_density(alpha=.2) +xlab("Node Variance")+
  theme(text = element_text(size=20))
png("hist_var.png", width = 6, height =6, units = 'in', res = 600)
print(hist4)
dev.off()

print("MultiNet: node histograms - end")

}

#difference of networks case/control

if (casenet==1 & controlnet==1){

  print("MultiNet differences between case/control networks - start")
  
  n2_case_global_rand<-sample(n2_case_global)
  n2_control_global_rand<-sample(n2_control_global)
  
  #case
  mean_case<-list()
  var_case<-list()
  for (i in 1:length(n2_case_global)){
    if (length(n2_case_global[[i]])>1){
      mean_case[[i]]<-mean(apply(all.dat_ord_global[n2_case_global[[i]],rcase],1,mean))
      var_case[[i]]<-mean(apply(all.dat_ord_global[n2_case_global[[i]],rcase],1,var))
    } else {
      mean_case[[i]]<-mean(as.numeric(all.dat_ord_global[n2_case_global[[i]],rcase]))
      var_case[[i]]<-var(as.numeric(all.dat_ord_global[n2_case_global[[i]],rcase]))
    }
  }
  
  #random normal distribution based on means and variances
  mean_case_all<-sum(unlist(mean_case))
  var_case_all<-sum(unlist(var_case))
  
  #control
  mean_control<-list()
  var_control<-list()
  for (i in 1:length(n2_control_global)){
    if (length(n2_control_global[[i]])>1){
      mean_control[[i]]<-mean(apply(all.dat_ord_global[n2_control_global[[i]],rcontrol],1,mean))
      var_control[[i]]<-mean(apply(all.dat_ord_global[n2_control_global[[i]],rcontrol],1,var))
    } else {
      mean_control[[i]]<-mean(as.numeric(all.dat_ord_global[n2_control_global[[i]],rcontrol]))
      var_control[[i]]<-var(as.numeric(all.dat_ord_global[n2_control_global[[i]],rcontrol]))
    }
  }
  
  #random normal distribution based on means and variances
  mean_control_all<-sum(unlist(mean_control))
  var_control_all<-sum(unlist(var_control))
  
  #difference among distributions
  norm_case<-rnorm(maxi,mean_case_all,var_case_all)
  norm_control<-rnorm(maxi,mean_control_all,var_control_all)
  # 
  # hist(norm_case)
  # hist(norm_control)
  
  median(norm_case)
  median(norm_control)
  
  # reset_par()
  # plot(ecdf(norm_case), col="blue",xlim = range(c(norm_case,norm_control)),main="")
  # plot(ecdf(norm_control), add = TRUE, lty = "dashed",col="red")

  #t.test(norm_case, norm_control)
  #wilcox.test(norm_case, norm_control)
  #ks.test(norm_case,norm_control) #biased by variances, better to not use it
  
  #test also the effect size to avoid biased p-values
  print("Cliffs delta")
  print(cliff.delta(d =norm_case,f = norm_control)) #with 0 indicating stochastic equality of the two groups. 

  print("MultiNet: differences between case/control networks - end")
}

##############################################
#DMSs=Differentially Methylated Sites
#SMSs=Similarly Methylated Sites
#HCSs=Highly Correlated Sites

if(biological==1){

print("MultiNet: DMS identification - start")

color_sam<-hue_pal()(length(unique(all.meta1$disease)))
names(color_sam) <- levels(factor(all.meta1$disease)) 
col1 = list(Status = color_sam)

# Create the heatmap annotation
difcg(dif1_all_global,c("1"),pointnode_all_global,all.dat_ord_global[,rownames(all.meta1)],rownames(all.meta1),pngname="global_diff_methyl.png",
      pngname1="global_diff_methyl_chr.png",pngname2="global_diff_methyl_sample.png")
dms_genes_global<-outgenes
dms_genes_global_cg<-genescgs
dms_global<-n3
n2_global<-n2

length(dms_global)
n_dat<-all.dat_ord_global[dms_global,]
mean_dat_case<-apply(n_dat[,rcase],1,median)
mean_dat_control<-apply(n_dat[,rcontrol],1,median)

#more relevant ones
mean_dat_diff_abs<-sort(abs(mean_dat_case-mean_dat_control),decreasing=TRUE)
mean_dat_diff_abs1<-mean_dat_diff_abs

mean_dat_diff<-(mean_dat_case-mean_dat_control)[names(mean_dat_diff_abs1)]
dat_box<-as.data.frame(cbind(mean_dat_diff,rownames(mean_dat_diff)))
head(dat_box)
dat_box$id<-names(mean_dat_diff)

print("MultiNet: DMS identification - end")

#normalize and logistic
print("MultiNet: logistic regression")

all.meta2<-all.meta1[c(rcase,rcontrol),]

for (i in 1:length(c(rcase,rcontrol))){
  if (rownames(all.meta2)[i] %in% rcase){
    all.meta2$class[i]=1
  } else if (rownames(all.meta2)[i] %in% rcontrol) {all.meta2$class[i]=0}
}

logidat<-as.data.frame(na.omit(t(all.dat_ord_global[dms_global,c(rcase,rcontrol)])))

#normalize extracting the mean
#mean_logi<-apply(logidat,2,mean)
#var_logi<-apply(logidat,2,sd)
#for (i in 1:dim(logidat)[2]){
#logidat[,i]<-(logidat[,i]-mean_logi[i])/var_logi[i]
#logidat[,i]<-(logidat[,i]-mean_logi[i])
#}

class0<-all.meta2$class
names(class0)<-c(rcase,rcontrol)
logidat$class<-factor(class0)

pval<-numeric()
for (i in colnames(logidat)[1:dim(logidat)[2]-1]){
  model = glm(class ~ get(i), data=logidat,family = binomial("logit"), maxit = 100)
  #summary(model)
  pval[[i]]<-coef(summary(model))[,4][2]
}

aux_w=warnings()

head(dat_box)
head(pval)
pval<-pval[rownames(dat_box)]
dat_box$pvalue<-pval

padj<-p.adjust(pval, method = "fdr", n = length(pval))

head(padj)
dat_box$padj<-padj[rownames(dat_box)]
dat_box_sig<-dat_box[dat_box$padj<0.05,]
dim(dat_box_sig)
dim(dat_box)
head(dat_box_sig)
dms_global_sig<-rownames(dat_box_sig)
length(dms_global_sig)

if (length(aux_w)>0){
print("Complete separation warning, proceed with penalized regression & SGoF")

#penalized logistic regression
pval_pen<-numeric()
for (i in colnames(logidat)[1:dim(logidat)[2]-1]){
  model.p<-logistf(data=logidat, class~get(i),pl=FALSE)
  pval_pen[[i]]<-model.p$prob[2]
}

#glm.probs = predict(model, newdata = logidat, type = "response")
#glm.pred = ifelse(glm.probs > 0.5, "Case", "Control")
#table(glm.pred, logidat$class)

pval_pen<-pval_pen[rownames(dat_box)]
dat_box$pvalue_pen<-pval_pen

#adjust pvalues
#alternative
m41 <- SGoF(pval_pen,alpha = 0.05, gamma = 0.05)
summary(m41)
sorted<-sort(pval_pen)
padj<-m41$Adjusted.pvalues
names(padj)<-names(sorted)

head(padj)
dat_box$padj<-padj[rownames(dat_box)]
dat_box_sig<-dat_box[dat_box$padj<0.05,]
dim(dat_box_sig)
dim(dat_box)
head(dat_box_sig)
dms_global_sig<-rownames(dat_box_sig)
length(dms_global_sig)
}

# Create the heatmap annotation
col_fun <- colorRampPalette(brewer.pal(10, "RdYlBu"))(256)
ha <- HeatmapAnnotation(Status = all.meta1[c(rcase,rcontrol),]$disease,col = col1)
h1<-Heatmap(as.matrix(all.dat_ord_global[dms_global_sig,c(rcase,rcontrol)]), name = "Methylation",clustering_distance_rows = "pearson",show_row_names=FALSE,show_column_names=FALSE,show_column_dend = TRUE,row_dend_reorder = TRUE,column_title = "",col = col_fun,column_km = 2, show_parent_dend_line = FALSE, top_annotation = ha)
reset_par()
png("heatmap_sig.png", width = 10, height =10, units = 'in', res = 300)
print(h1)
dev.off()

chr1<-as.numeric(substr(annotation[dms_global_sig,]$chr,4,6))
sort(chr1)
annotation_chr<-cbind(annotation[dms_global_sig,],chr1)
chr1<-unique(sort(chr1))
annotation_chr<-annotation_chr[order(annotation_chr$chr1,annotation_chr$pos),]
head(annotation_chr)
posir<-rownames(annotation_chr)
dat_box_sig$chr<-annotation_chr[rownames(dat_box_sig),34]
dat_box_sig<-dat_box_sig[order(dat_box_sig$chr),]

rel<-ggdotchart(dat_box_sig, x = "chr", y ="mean_dat_diff",
                color = "chr", palette = "raiwnbow", size = 3, 
                add = "segment", 
                sorting = "ascending",    
                add.params = list(color = "black", size = 1.5),
                position = position_dodge(0.3),
                ggtheme = theme_gray(),
                xlab="Chromosome",
                ylab="Difference median methylation"
)+ theme(legend.position = "none",text = element_text(size=15)) 
png("dms_global_relevance.png", width = 6, height =6, units = 'in', res = 600)
print(rel)
dev.off()

logidat_ad<-logidat[,dms_global_sig]
logidat_ad$class<-factor(class0)

print("MultiNet: random forest")

set.seed(12345)
modelo <- randomForest(class ~ ., data=logidat)
modelo

relev<-as.data.frame(modelo$importance[order(modelo$importance,decreasing=TRUE),])
head(relev)
relev$id<-rownames(relev)
colnames(relev)<-c("importance","id")
relev$pvalue<-padj[rownames(relev)]
relev$diff<-mean_dat_diff[rownames(relev)]

#plot importance
# reset_par()
# varImpPlot(modelo,n.var=100)

relev_test<-relev[relev$pvalue>0.05,] #they match with low importance

relev0<-relev[relev$pvalue<0.05,]
relev1<-relev[relev$importance>=0 & relev$pvalue<0.05,]
dim(relev1)

#select the first ones to plot
relev_s<-relev1[1:50,]

rforest1<-ggplot(relev_s, aes(x=id, y=diff)) +
  geom_segment( aes(x=id, xend=id, y=0.1, yend=diff), color="black") +
  geom_point( color="blue", size=5) +
  coord_flip() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank(),
    text = element_text(size=10)
  ) +
  xlab("CGs") +
  ylab("Methylation difference")
png("dms_global_cgs_rforest.png", width = 8, height =8, units = 'in', res = 600)
print(rforest1)
dev.off()

dms_genes_global_cg<-bargen(rownames(relev_s),all.dat_ord_global,0)[[3]]
relev_s$gene<-dms_genes_global_cg[rownames(relev_s)]
head(relev_s)

rforest2<-ggplot(na.omit(relev_s), aes(x=gene, y=diff)) +
  geom_segment( aes(x=gene, xend=gene, y=0, yend=diff), color="black") +
  geom_point( color="blue", size=5) +
  coord_flip() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank(),
    text = element_text(size=10)
  ) +
  xlab("Genes") +
  ylab("Methylation difference")
png("dms_global_genes_rforest.png", width = 8, height =8, units = 'in', res = 600)
print(rforest2)
dev.off()

#predictions
#predicciones <- predict(modelo, test)
#t <- with(test, table(predicciones, class)) 
#t ; 100 * sum(diag(t)) / sum(t)

#heatmaps
mean_dat_diff_abs<-sort(mean_dat_case[dms_global_sig]-mean_dat_control[dms_global_sig],decreasing=TRUE)

dms_global_hyper<-mean_dat_diff_abs[mean_dat_diff_abs>0]
dms_global_hypo<-mean_dat_diff_abs[mean_dat_diff_abs<=0]

gene.n  <- length(mean_dat_diff_abs);
chr     <- names(mean_dat_diff_abs);
val     <- mean_dat_diff_abs;
seg.dat <- data.frame(chr=annotation_chr[dms_global_sig,]$chr, start=annotation_chr[dms_global_sig,]$pos, end=annotation_chr[dms_global_sig,]$pos+1, 
                      name=annotation_chr[dms_global_sig,]$chr1, value1=mean_dat_diff_abs,value2=mean_dat_case[dms_global_sig],value3= mean_dat_control[dms_global_sig]);
#seg.c   <- segAnglePo(seg.dat, seg=annotation_chr$chr1);
cols    <- rainbow(22, alpha=0.8);

#circle plot
#points hyper/hypo
reset_par()
png("circle.png", width = 8, height =8, units = 'in', res = 600)
par(mar = c(1, 1, 1, 1))
circos.initializeWithIdeogram(chromosome.index = paste0("chr", seq(1,22,1)))
seg.dat.1.hyper<-seg.dat[names(dms_global_hyper),c(1,2,3)]
seg.dat.1.hypo<-seg.dat[names(dms_global_hypo),c(1,2,3)]
bed_list = list(seg.dat.1.hyper, seg.dat.1.hypo)
circos.genomicRainfall(bed_list, pch = 16, cex = 0.4, col = c("#FF000080", "#0000FF80"))
circos.genomicDensity(bed_list[[1]], col = c("#FF000080"), track.height = 0.1)
circos.genomicDensity(bed_list[[2]], col = c("#0000FF80"), track.height = 0.1)
dev.off()

#heatmaps
print("MultiNet: correlation heatmaps - start")

dms_cor<-na.omit(cor(t(all.dat_ord_global[dms_global_sig,rcase])))
dim(dms_cor)

#by chromosome
chr1<-as.numeric(substr(annotation[dms_global_sig,]$chr,4,6))
sort(chr1)
annotation_chr<-cbind(annotation[dms_global_sig,],chr1)
chr1<-unique(sort(chr1))
annotation_chr<-annotation_chr[order(annotation_chr$chr1,annotation_chr$pos),]
head(annotation_chr)

dms_cor_chr<-matrix(NA, 22, 22)
for (i in 1:22){
  for (j in 1:22){
    annotation_1<-annotation_chr[annotation_chr$chr1==i,]
    length(annotation_1)
    
    annotation_2<-annotation_chr[annotation_chr$chr1==j,]
    length(annotation_2)
    
    #select correlation for this chromosome
    dms_cor1<-dms_cor[rownames(annotation_1),rownames(annotation_2)]
    dim(dms_cor1)
    
    dms_cor_chr[i,j]<-median(abs(dms_cor1))
  }
}

diag(dms_cor_chr)<-1
rownames(dms_cor_chr)<-seq(1,22,1)
col_fun <- colorRampPalette(brewer.pal(10, "RdYlBu"))(256)

reset_par()
h2<-Heatmap(dms_cor_chr, name = "Median Correlation",row_order=rownames(dms_cor_chr),column_order=colnames(dms_cor_chr),show_row_names=FALSE,show_column_names=FALSE,col=col_fun) +
  Heatmap(as.numeric(rownames(dms_cor_chr)), name = "Chromosome", width = unit(5, "mm"),
          col =  rainbow(length(rownames(dms_cor_chr))))
png("heatmap_cor_diffcg_chr_case.png", width = 10, height =10, units = 'in', res = 300)
print(h2)
dev.off()

dms_cor<-na.omit(cor(t(all.dat_ord_global[dms_global_sig,rcontrol])))
dim(dms_cor)

#by chromosome
chr1<-as.numeric(substr(annotation[dms_global_sig,]$chr,4,6))
sort(chr1)
annotation_chr<-cbind(annotation[dms_global_sig,],chr1)
chr1<-unique(sort(chr1))
annotation_chr<-annotation_chr[order(annotation_chr$chr1,annotation_chr$pos),]
head(annotation_chr)

dms_cor_chr<-matrix(NA, 22, 22)
for (i in 1:22){
  for (j in 1:22){
    annotation_1<-annotation_chr[annotation_chr$chr1==i,]
    length(annotation_1)
    
    annotation_2<-annotation_chr[annotation_chr$chr1==j,]
    length(annotation_2)
    
    #select correlation for this chromosome
    dms_cor1<-dms_cor[rownames(annotation_1),rownames(annotation_2)]
    dim(dms_cor1)
    
    dms_cor_chr[i,j]<-median(abs(dms_cor1))
  }
}

diag(dms_cor_chr)<-1
rownames(dms_cor_chr)<-seq(1,22,1)
col_fun <- colorRampPalette(brewer.pal(10, "RdYlBu"))(256)

reset_par()
h2<-Heatmap(dms_cor_chr, name = "Median Correlation",row_order=rownames(dms_cor_chr),column_order=colnames(dms_cor_chr),show_row_names=FALSE,show_column_names=FALSE,col=col_fun) +
  Heatmap(as.numeric(rownames(dms_cor_chr)), name = "Chromosome", width = unit(5, "mm"),
          col =  rainbow(length(rownames(dms_cor_chr))))
png("heatmap_cor_diffcg_chr_con.png", width = 10, height =10, units = 'in', res = 300)
print(h2)
dev.off()

dms_cor<-na.omit(cor(t(all.dat_ord_global[dms_global_sig,])))
dim(dms_cor)

#by chromosome
chr1<-as.numeric(substr(annotation[dms_global_sig,]$chr,4,6))
sort(chr1)
annotation_chr<-cbind(annotation[dms_global_sig,],chr1)
chr1<-unique(sort(chr1))
annotation_chr<-annotation_chr[order(annotation_chr$chr1,annotation_chr$pos),]
head(annotation_chr)

dms_cor_chr<-matrix(NA, 22, 22)
for (i in 1:22){
  for (j in 1:22){
    annotation_1<-annotation_chr[annotation_chr$chr1==i,]
    length(annotation_1)
    
    annotation_2<-annotation_chr[annotation_chr$chr1==j,]
    length(annotation_2)
    
    #select correlation for this chromosome
    dms_cor1<-dms_cor[rownames(annotation_1),rownames(annotation_2)]
    dim(dms_cor1)
    
    dms_cor_chr[i,j]<-median(abs(dms_cor1))
  }
}

diag(dms_cor_chr)<-1
rownames(dms_cor_chr)<-seq(1,22,1)
col_fun <- colorRampPalette(brewer.pal(10, "RdYlBu"))(256)

reset_par()
h2<-Heatmap(dms_cor_chr, name = "Median Correlation",row_order=rownames(dms_cor_chr),column_order=colnames(dms_cor_chr),show_row_names=FALSE,show_column_names=FALSE,col=col_fun) +
  Heatmap(as.numeric(rownames(dms_cor_chr)), name = "Chromosome", width = unit(5, "mm"),
          col =  rainbow(length(rownames(dms_cor_chr))))
png("heatmap_cor_diffcg_chr.png", width = 10, height =10, units = 'in', res = 300)
print(h2)
dev.off()

print("MultiNet: correlation heatmaps - end")

print("MultiNet: DMSs regions and trends - start")

png("dms_genes_global.png", width = 10, height =10, units = 'in', res = 600)
reset_par()
par(mfrow=c(2,1))
barcor<-bargen(dms_global_sig,all.dat_ord_global,0)
dev.off()

png("dms_global_site.png", width = 10, height =10, units = 'in', res = 600)
reset_par()
sitegen(names(dms_global_hyper),names(dms_global_hypo),all.dat_ord_global,0)
dev.off()

png("dms_global_island.png", width = 10, height =10, units = 'in', res = 600)
reset_par()
islandgen(names(dms_global_hyper),names(dms_global_hypo),all.dat_ord_global,0)
dev.off()

#genes
png("dms_genes_global_hyper.png", width = 10, height =10, units = 'in', res = 600)
reset_par()
par(mfrow=c(2,1))
barcor<-bargen(names(dms_global_hyper),all.dat_ord_global,0)
dev.off()

png("dms_genes_global_hypo.png", width = 10, height =10, units = 'in', res = 600)
reset_par()
par(mfrow=c(2,1))
barcor<-bargen(names(dms_global_hypo),all.dat_ord_global,0)
dev.off()

#hyper meth
difcg(exp1_all_global,c("1"),pointnode_all_global,all.dat_ord_global,c(rcase,rcontrol),pngname="global_hyper_methyl.png",
      pngname1="global_hyper_methyl_chr.png",pngname2="global_hyper_methyl_sample.png")
genes_global_all_hyper<-outgenes
cg_global_all_hyper<-n3

#hypo meth
difcg(exp1_all_global,c("3"),pointnode_all_global,all.dat_ord_global,c(rcase,rcontrol),pngname="global_hypo_methyl.png",
      pngname1="global_hypo_methyl_chr.png",pngname2="global_hypo_methyl_sample.png")
genes_global_all_hypo<-outgenes
cg_global_all_hypo<-n3

##################
#manhattan plot
sel<-dms_global_sig
database=all.dat_ord_global

anotcor<-annotation[sel,]

chr_m<-as.numeric(substr(anotcor$chr,4,6))
names(chr_m)<-sel
pos_m<-annotation[sel,2]
names(pos_m)<-sel

n_dat<-database[sel,]
mean_dat_case_sel<-mean_dat_case[sel]
mean_dat_control_sel<-mean_dat_control[sel]
mean_dat_diff_sel<-(mean_dat_case_sel-mean_dat_control_sel)

man<-as.data.frame(cbind(as.numeric(mean_dat_case_sel),as.numeric(mean_dat_control_sel),chr_m,pos_m,as.numeric(mean_dat_diff_sel)))
head(man)
man$CG<-rownames(man)
colnames(man)<-c("CASE","CONTROL","CHR","POS","DIFF","CG")

man<-man[order(man[,3],man[,4]),]
head(man)

png("manhattan_global.png", width = 8, height =8, units = 'in', res = 600)
par(mfrow=c(2,1))
manhattan_methyl(man, chr="CHR", bp="POS", p="CASE" ,snp="CG",suggestiveline = c(0.5),genomewideline =FALSE,logp = FALSE,col = rainbow(22),main="Group 1")
manhattan_methyl(man, chr="CHR", bp="POS", p="CONTROL" ,snp="CG",suggestiveline = c(0.5),genomewideline =FALSE,logp = FALSE,col = rainbow(22),main="Group 2")
dev.off()

#annotate genes names
dms_genes_global_cg<-bargen(names(dms_global_hypo),all.dat_ord_global,0)[[3]]
glist<-data.frame(dms_genes_global_cg,stringsAsFactors=FALSE)
glist$ID<-names(dms_genes_global_cg)
head(glist)

glist1<-glist[!duplicated(glist$ID),]
man1<-man[glist1$ID,]
hdata <- cbind(glist1,man1)
colnames(hdata)<-c("GENE","ID","CASE","CONTROL","CHR","POS","DIFF","CG")
glist <- as.character(hdata$dms_genes_global_cg)
hdata<-na.omit(hdata)

png("manhattan_global_diff.png", width = 6, height =6, units = 'in', res = 900)
reset_par()
manhattan_methyl(hdata, chr="CHR", bp="POS", p="DIFF" ,snp="GENE",suggestiveline = 0.0000001,
                 genomewideline =FALSE,logp = FALSE,col = rainbow(22),main="Case-Control Median - Annotated genes",
                 annotatePvalmax=0.4,annotatePvalmin=-0.4)
dev.off()

print("MultiNet: DMSs regions and trends - end")

###########
#pathways 

print("MultiNet: DMSs pathways - start")

#common marks
genes_sig<-bargen(dms_global_sig,annotation0,0)
dms_genes_global_sig<-genes_sig[[1]]
eg_dmr_global<- bitr(unique(dms_genes_global_sig), fromType="SYMBOL", toType=c("ENTREZID", "ENSEMBL"), OrgDb="org.Hs.eg.db")
head(eg_dmr_global)
genelist0_dmr_global<-as.numeric(eg_dmr_global[,2])
genelist_dmr_global<-genelist0_dmr_global
names(genelist_dmr_global) = as.character(eg_dmr_global[,1])
# 
# genes_sig<-bargen(names(dms_global_hypo),annotation0,0)
# dms_genes_global_sig<-genes_sig[[1]]
# eg_dmr_global<- bitr(unique(dms_genes_global_sig), fromType="SYMBOL", toType=c("ENTREZID", "ENSEMBL"), OrgDb="org.Hs.eg.db")
# head(eg_dmr_global)
# genelist0_dmr_global<-as.numeric(eg_dmr_global[,2])
# genelist_dmr_global2<-genelist0_dmr_global
# names(genelist_dmr_global2) = as.character(eg_dmr_global[,1])

#sample<-list(Age=genelist_dmr_global, Prostate=genelist_dmr_global2,Colorectal=genelist_dmr_global3) 

#sample<-list(Hyper=genelist_dmr_global, Hypo=genelist_dmr_global2) 

#enr10 <- compareCluster(geneCluster = sample, fun = "enrichKEGG")
enr10<-enrichKEGG(genelist_dmr_global)
if (dim(enr10)[1]>0){
enr1<-list()
enr1<-clusterProfiler::dotplot(enr10,showCategory = 20,  font.size = 8)
png("enrichKEGG.png", width = 10, height =20, units = 'in', res = 600)
print(enr1)
dev.off()
}else {print("No KEGG enrichment found")}

#enr20 <- compareCluster(geneCluster = sample, fun = "enrichDGN")
enr20<-enrichDGN(genelist_dmr_global)
if (dim(enr20)[1]>0){
enr2<-list()
enr2<-clusterProfiler::dotplot(enr20,showCategory = 20,  font.size = 8)
png("enrichDGN.png", width = 10, height =10, units = 'in', res = 600)
print(enr2)
dev.off()
}else {print("No DGN enrichment found")}

#enr30 <- compareCluster(geneCluster = sample, fun = "enrichGO",OrgDb = org.Hs.eg.db)
#head(as.data.frame(enr30))
enr30<-enrichGO(genelist_dmr_global,OrgDb = org.Hs.eg.db)
if (dim(enr30)[1]>0){
enr3<-list()
enr3<-clusterProfiler::dotplot(enr30,showCategory = 20,  font.size = 8)
png("enrichGO.png", width = 10, height =10, units = 'in', res = 600)
print(enr3)
dev.off()
} else {print("No GO enrichment found")}

#graphics of enrichment
if (dim(enr10)[1]>1 & dim(enr20)[1]>1){

edox <- setReadable(enr10, 'org.Hs.eg.db', 'ENTREZID')
edo<-setReadable(enr20, 'org.Hs.eg.db', 'ENTREZID')
p1 <- cnetplot(edox, node_label="category") 
p2 <- cnetplot(edox, node_label="gene") 
png("enrichKEGG_graph.png", width = 10, height =10, units = 'in', res = 600)
cowplot::plot_grid(p1, p2, ncol=2, labels=LETTERS[1:4])
dev.off()

p1 <- cnetplot(edo, node_label="category") 
p2 <- cnetplot(edo, node_label="gene") 
png("enrichdis_graph.png", width = 10, height =10, units = 'in', res = 600)
cowplot::plot_grid(p1, p2, ncol=2, labels=LETTERS[1:4])
dev.off()

png("enrichKEGG_heat.png", width = 10, height =10, units = 'in', res = 600)
p1<-heatplot(edox)
p2<-heatplot(edo)
cowplot::plot_grid(p1, p2,ncol=1, labels=LETTERS[1:2])
dev.off()

png("enrichKEGG_net.png", width = 10, height =10, units = 'in', res = 600)
p1<-emapplot(edox)
p2<-emapplot(edo)
cowplot::plot_grid(p1, p2,ncol=2, labels=LETTERS[1:2])
dev.off()
}

print("MultiNet: DMSs pathways - end")

if (casenet==1){
  
  # print("MultiNet: case/control correlation heatmaps - start")
  
  # col1 = list(Status = c("Newborns" = color_sam[1]))
  # 
  # #difference
  # difcg(exp1_case_global,c("1"),pointnode_case_global,all.dat_ord_global,rcase,pngname="globalcase_hyper_methyl.png",
  #       pngname1="globalcase_hyper_methyl_chr.png",pngname2="globalcase_hyper_methyl_sample.png")
  # genes_global_case_hyper<-outgenes
  # cg_global_case_hyper<-n3
  # 
  # difcg(exp1_case_global,c("3"),pointnode_case_global,all.dat_ord_global,rcase,pngname="globalcase_hypo_methyl.png",
  #       pngname1="globalcase_hypo_methyl_chr.png",pngname2="globalcase_hypo_methyl_sample.png")
  # genes_global_case_hypo<-outgenes
  # cg_global_case_hypo<-n3
  # 
  #correlation
  # difcg(cor1_med_case_global,c("1"),pointnode_case_global,all.dat_ord_global,rcase,pngname="globalcase_hyper_corr.png",
  #       pngname1="globalcase_hyper_corr_chr.png",pngname2="globalcase_hyper_corr_sample.png")
  # genes_global_case_hyper_corr<-outgenes
  # cg_global_case_hyper_corr<-n3
  # 
  # cor_global_case<-cor(t(all.dat_ord_global[cg_global_case_hyper_corr,rcase]))
  # dim(cor_global_case)
  # median(abs(cor_global_case))
  # col_fun <- colorRampPalette(brewer.pal(10, "RdYlBu"))(256)
  # png("cor_cor_case.png", width = 6, height =6, units = 'in', res = 600)
  # Heatmap(cor_global_case,name="correlation",show_row_names = FALSE,show_column_names = FALSE,col=col_fun)
  # dev.off()
  
}

if (controlnet==1){
  # col1 = list(Status = c("Nonagenarians" = color_sam[2]))
  # 
  # difcg(exp1_control_global,c("1"),pointnode_control_global,all.dat_ord_global,rcontrol,pngname="globalcontrol_hyper_methyl.png",
  #       pngname1="globalcontrol_hyper_methyl_chr.png",pngname2="globalcontrol_hyper_methyl_sample.png")
  # genes_global_control_hyper<-outgenes
  # cg_global_control_hyper<-n3
  # 
  # difcg(exp1_control_global,c("3"),pointnode_control_global,all.dat_ord_global,rcontrol,pngname="globalcontrol_hypo_methyl.png",
  #       pngname1="globalcontrol_hypo_methyl_chr.png",pngname2="globalcontrol_hypo_methyl_sample.png")
  # genes_global_control_hypo<-outgenes
  # cg_global_control_hypo<-n3
  # 
  # int_hyper<-length(intersect(cg_global_case_hyper,cg_global_control_hypo))
  # int_hypo<-length(intersect(cg_global_case_hypo,cg_global_control_hyper))
  # 
  #correlation
  # difcg(cor1_med_control_global,c("1"),pointnode_control_global,all.dat_ord_global,rcontrol,pngname="globalcontrol_hypo_corr.png",
  #       pngname1="globalcontrol_hyper_corr_chr.png",pngname2="globalcontrol_hyper_corr_sample.png")
  # genes_global_control_hyper_corr<-outgenes
  # cg_global_control_hyper_corr<-n3
  # 
  # cor_global_control<-cor(t(all.dat_ord_global[cg_global_control_hyper_corr,c(rcase,rcontrol)]))
  # dim(cor_global_control)
  # median(abs(cor_global_control))
  # col_fun <- colorRampPalette(brewer.pal(10, "RdYlBu"))(256)
  # png("cor_cor_control.png", width = 6, height =6, units = 'in', res = 600)
  # Heatmap(cor_global_control,name="correlation",show_row_names = FALSE,show_column_names = FALSE,col=col_fun)
  # dev.off()
  
  #print("MultiNet: case/control correlation heatmaps - end")
}

  #create an excel sheet with all the information
  
  print("MultiNet: create spreadsheet - start")
  
  final_data<-annotation[dms_global,]
  final_data$log_pval<-dat_box[dms_global,]$pvalue
  #final_data$gene<-dms_genes_global_cg[dms_global]
  final_data$log_adjpval<-dat_box[dms_global,]$padj
  final_data$forest<-relev[dms_global,]$importance
  head(final_data)
  dim(final_data)
  
  mean_dat_diff_abs_all<-sort(mean_dat_case[dms_global]-mean_dat_control[dms_global],decreasing=TRUE)
  length(mean_dat_diff_abs_all)
  dms_global_hyper_all<-mean_dat_diff_abs_all[mean_dat_diff_abs_all>0]
  dms_global_hypo_all<-mean_dat_diff_abs_all[mean_dat_diff_abs_all<=0]
  length(dms_global_hyper_all)
  length(dms_global_hypo_all)
  
  for (i in 1:length(dms_global)){
    
    if (final_data$log_adjpval[i] < 0.05){
      final_data$signif[i]<-"Significant"
    } else {final_data$signif[i]<-"Non-significant"}
    
    if (dms_global[i] %in% names(dms_global_hyper_all)){
      final_data$hyperhypo[i]<-"Hyper"
    } else if (dms_global[i] %in% names(dms_global_hypo_all)){
      final_data$hyperhypo[i]<-"Hypo"
    }
    
    if (dms_global[i] %in% rownames(relev_s)){
      final_data$relev[i]<-1
    } else {final_data$relev[i]<-0}
    
  }
  
  write.csv(final_data, "dms_global_dif.csv",fileEncoding = "UTF-16LE")
  
  print("MultiNet: create spreadsheet - end")
  
  print("SUMMARY")
  print("TOTAL DMSs: ")
  print(length(dms_global))
  print("TOTAL SIGNIFICANT DMSs: ")
  print(length(dms_global_sig))
  print("TOTAL SIGNIFICANT DMSs HYPERMETHYLATED: ")
  print(length(dms_global_hyper))
  print("TOTAL SIGNIFICANT DMSs HYPOMETHYLATED: ")
  print(length(dms_global_hypo))
  print("MOST RELEVANT GENES: ")
  print(dms_genes_global_cg)
}
  end_time2 <- Sys.time()
  
  total_time2<-end_time2 - start_time

  print("MultiNet TOTAL RUNNING TIME")
  print(total_time2)
  
  #remove warnings
  #assign("last.warning", NULL, envir = baseenv())
})
