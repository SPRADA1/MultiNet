###################################################################
## MultiNet algorithm - Suppelemental 1 - Other functions
## Part III - Chapter 8 
## Sara Prada, October 2020
###################################################################

####################
#Auxiliary functions

intersect_all <- function(a,b){
  all_data <- c(a,b)
  count_data<- length(list(a,b))
  #freq_dist <- count(all_data)
  freq_dist <- as.data.frame(table(all_data))
  intersect_data <- freq_dist[which(freq_dist$Freq==count_data),]
  intersect_data
}

`%not in%` <- function (x, table) is.na(match(x, table, nomatch=NA_integer_))

outersect <- function(x, y) {
  sort(c(setdiff(x, y),
         setdiff(y, x)))
}

flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}

showCols <- function(cl=colors(), bg = "grey",
                     cex = 0.75, rot = 30) {
  m <- ceiling(sqrt(n <-length(cl)))
  length(cl) <- m*m; cm <- matrix(cl, m)
  require("grid")
  grid.newpage(); vp <- viewport(w = .92, h = .92)
  grid.rect(gp=gpar(fill=bg))
  grid.text(cm, x = col(cm)/m, y = rev(row(cm))/m, rot = rot,
            vp=vp, gp=gpar(cex = cex, col = cm))
}

reset_par <- function(){
  op <- structure(list(xlog = FALSE, ylog = FALSE, adj = 0.5, ann = TRUE,
                       ask = FALSE, bg = "transparent", bty = "o", cex = 1, cex.axis = 1,
                       cex.lab = 1, cex.main = 1.2, cex.sub = 1, col = "black",
                       col.axis = "black", col.lab = "black", col.main = "black",
                       col.sub = "black", crt = 0, err = 0L, family = "", fg = "black",
                       fig = c(0, 1, 0, 1), fin = c(6.99999895833333, 6.99999895833333
                       ), font = 1L, font.axis = 1L, font.lab = 1L, font.main = 2L,
                       font.sub = 1L, lab = c(5L, 5L, 7L), las = 0L, lend = "round",
                       lheight = 1, ljoin = "round", lmitre = 10, lty = "solid",
                       lwd = 1, mai = c(1.02, 0.82, 0.82, 0.42), mar = c(5.1, 4.1,
                                                                         4.1, 2.1), mex = 1, mfcol = c(1L, 1L), mfg = c(1L, 1L, 1L,
                                                                                                                        1L), mfrow = c(1L, 1L), mgp = c(3, 1, 0), mkh = 0.001, new = FALSE,
                       oma = c(0, 0, 0, 0), omd = c(0, 1, 0, 1), omi = c(0, 0, 0,
                                                                         0), pch = 1L, pin = c(5.75999895833333, 5.15999895833333),
                       plt = c(0.117142874574832, 0.939999991071427, 0.145714307397962,
                               0.882857125425167), ps = 12L, pty = "m", smo = 1, srt = 0,
                       tck = NA_real_, tcl = -0.5, usr = c(0.568, 1.432, 0.568,
                                                           1.432), xaxp = c(0.6, 1.4, 4), xaxs = "r", xaxt = "s", xpd = FALSE,
                       yaxp = c(0.6, 1.4, 4), yaxs = "r", yaxt = "s", ylbias = 0.2), .Names = c("xlog",
                                                                                                "ylog", "adj", "ann", "ask", "bg", "bty", "cex", "cex.axis",
                                                                                                "cex.lab", "cex.main", "cex.sub", "col", "col.axis", "col.lab",
                                                                                                "col.main", "col.sub", "crt", "err", "family", "fg", "fig", "fin",
                                                                                                "font", "font.axis", "font.lab", "font.main", "font.sub", "lab",
                                                                                                "las", "lend", "lheight", "ljoin", "lmitre", "lty", "lwd", "mai",
                                                                                                "mar", "mex", "mfcol", "mfg", "mfrow", "mgp", "mkh", "new", "oma",
                                                                                                "omd", "omi", "pch", "pin", "plt", "ps", "pty", "smo", "srt",
                                                                                                "tck", "tcl", "usr", "xaxp", "xaxs", "xaxt", "xpd", "yaxp", "yaxs",
                                                                                                "yaxt", "ylbias"))
  par(op)
}

################################################################
#Apply the kind of Mapper algorithm to each window from a to a0

corr <- defmacro(a, a0, filter1_global1=filter1_global,filter2_global1=filter2_global,
                filter1_control1=filter1_control,filter2_control1=filter2_control,
                filter1_case1=filter1_case,filter2_case1=filter2_case,
                 expr = {
  
  #Data from the window selected a:a0
  all.dat.window_global<-all.dat_ord_global[a:a0,]
  tall.dat.window_global <- t(all.dat.window_global)
  dim(all.dat.window_global)
  head(all.dat.window_global)
  var_windowr_global<-rownames(all.dat.window_global)
  
  #filter1
  medi_global_w<-filter1_global1[var_windowr_global]
  
  #variance
  #var_s1_w<-var_s1[var_windowr_global]
  
  #filter2
  filter1_window_global <- sort(na.omit(filter2_global1[var_windowr_global]),decreasing=TRUE)
  
  #Correlation distance
  corwindow_global<-cor(tall.dat.window_global)
  head(corwindow_global)
  dim(corwindow_global)
  #diag(corwindow_global)<-0
  
  #Apply a kind of Mapper with our own function
  #object: methylation matrix
  #dist_object: correlation matrix
  #filter_values: filter functions
  #num_intervals: number of intervals per window
  #percent_overlap: percentage of intervals overlap
  #num_bins_when_clustering: number of clusters per interval
  #meth: clustering method
  #pero0: rate of points in common among two nodes needed to join them
  #corm0: mean correlation needed to join to nodes
  
  cg_aux<-intersect(names(sort(medi_global_w)),names(filter1_window_global))
  m_global <- mapper_corr_w(object=tall.dat.window_global,
                            dist_object = corwindow_global, 
                            filter_values = cbind(sort(medi_global_w[cg_aux]),filter1_window_global[cg_aux]),
                            #filter_values = cbind(sort(medi_global_w[cg_aux]),var_s1_w[cg_aux]),
                            num_intervals = c(int,int),
                            percent_overlap = c(ov_i,ov_i),
                            num_bins_when_clustering = opt_cluster,
                            meth=method_map,
                            pero0=pero,
                            corm0=corm) 
  
  #graphs
  graph_global <- graph.adjacency(m_global$adjacency, mode="undirected",weighted=TRUE)
  
  graphic_global[[a0]]<-graph_global
  mapp_global[[a0]]<-m_global
  
  #nodes and edges
  MapperNodes_global <- na.omit(mapperVertices(m_global, 1:length(var_windowr_global)))
  MapperLinks_global <- na.omit(mapperEdges_w(m_global))
  
  links_global[[a0]]<-MapperLinks_global
  nodes_global[[a0]]<-MapperNodes_global 
  
  #Control
  if (controlnet==1) {
    
    #filter 1
    medi_global_control_w<-filter1_control1[var_windowr_global]
    
    #filter 2
    filter2_control_window_global <- filter2_control1[var_windowr_global] 
    
    all.dat.window_control_global <- all.dat_ord_global[var_windowr_global,rcontrol]
    tall.dat.window_control_global <- t(all.dat.window_control_global)
    
    corwindow_control_global<-cor(tall.dat.window_control_global)
    head(corwindow_control_global)
    dim(corwindow_control_global)
    #diag(corwindow_control_global)<-0
    
    cg_aux<-intersect(names(na.omit(sort(medi_global_control_w))), names(na.omit(filter2_control_window_global)))
    
    m_control_global <- mapper_corr_w(object=tall.dat.window_control_global,
                                      dist_object = corwindow_control_global, 
                                      filter_values =  cbind(sort(medi_global_control_w[cg_aux]),filter2_control_window_global[cg_aux]),
                                      num_intervals = c(int,int),
                                      percent_overlap = c(ov_i,ov_i),
                                      num_bins_when_clustering = opt_cluster,
                                      meth=method_map,
                                      pero0=pero,
                                      corm0=corm) 
    
    m.graph_control_global <- graph.adjacency(m_control_global$adjacency, mode="undirected",weighted=TRUE)
    
    graphic_control_global[[a0]]<-m.graph_control_global
    mapp_control_global[[a0]]<-m_control_global
    
    #nodes and edges 
    MapperNodes_control_global <- na.omit(mapperVertices(m_control_global, 1:length(var_windowr_global)))
    MapperLinks_control_global <- na.omit(mapperEdges_w(m_control_global))
    
    links_control_global[[a0]]<-MapperLinks_control_global
    nodes_control_global[[a0]]<-MapperNodes_control_global
  }
  
  #case
  if (casenet==1) {
    
    #filter 1
    medi_global_case_w<-filter1_case1[var_windowr_global]
    
    #filter 2
    filter2_case_window_global <- filter2_case1[var_windowr_global] 
    
    all.dat.window_case_global <- all.dat_ord_global[var_windowr_global,rcase]
    tall.dat.window_case_global <- t(all.dat.window_case_global)
    
    corwindow_case_global<-cor(tall.dat.window_case_global)
    head(corwindow_case_global)
    dim(corwindow_case_global)
    #diag(corwindow_case_global)<-0
    
    cg_aux<-intersect(names(na.omit(sort(medi_global_case_w))), names(na.omit(filter2_case_window_global)))
    
    m_case_global <- mapper_corr_w(object=tall.dat.window_case_global,
                                   dist_object = corwindow_case_global, 
                                   filter_values =cbind(sort(medi_global_case_w[cg_aux]),filter2_case_window_global[cg_aux]),
                                   num_intervals = c(int,int),
                                   percent_overlap = c(ov_i,ov_i),
                                   num_bins_when_clustering = opt_cluster,
                                   meth=method_map,
                                   pero0=pero,
                                   corm0=corm) 
    
    m.graph_case_global <- graph.adjacency(m_case_global$adjacency, mode="undirected",weighted=TRUE)
    
    graphic_case_global[[a0]]<-m.graph_case_global
    mapp_case_global[[a0]]<-m_case_global
    
    #nodes and edges 
    MapperNodes_case_global <- na.omit(mapperVertices(m_case_global, 1:length(var_windowr_global)))
    MapperLinks_case_global <- na.omit(mapperEdges_w(m_case_global))
    
    links_case_global[[a0]]<-MapperLinks_case_global
    nodes_case_global[[a0]]<-MapperNodes_case_global
    
  }
  
  print(a0) 
})

#################################
#Join networks obtained with corr

corr1 <- defmacro(mapp, graphic, links, nodes, database, expr = {
  
  #new node index adapted to each window
  a1<-seq(init+by+ov,maxi,ov)
  for (i in a1){
    N=length(mapp[[i]]$points_in_vertex)
    for (j in 1:N){
      for (k in 1:length(mapp[[i]]$points_in_vertex[[j]])){
        mapp[[i]]$points_in_vertex[[j]][k] <- mapp[[i]]$points_in_vertex[[j]][k] + (i-by) - 1
      }}}
  
  N=length(mapp[[init+by]]$points_in_vertex)
  for (j in 1:N){
    for (k in 1:length(mapp[[init+by]]$points_in_vertex[[j]])){
      mapp[[init+by]]$points_in_vertex[[j]][k] <- mapp[[init+by]]$points_in_vertex[[j]][k] 
    }}
  

  a<-seq(init+by+ov,maxi,ov)
  graphic1<-list()
  mapp1<-list()
  links1<-list()
  nodes1<-list()
  links1[[init+by]]<-links[[init+by]]
  graphic1[[init+by]]<-graph.data.frame(links1[[init+by]], directed=FALSE)
  mapp1[[init+by]]<-mapp[[init+by]]
  nodes1[[init+by]]<-nodes[[init+by]]
  
  for (i in a){
    mapp1[[i]]=c(mapp[i],mapp1[i-ov])
    links1[[i]]=rbind(links[[i]],links1[[i-ov]])
    nodes1[[i]]=rbind(nodes[[i]],nodes1[[i-ov]])
  }
  
  #we take points in common among nodes in different intervals and windows
  a1<-seq(init+by+ov,maxi,ov)
  links2<-list()
  
  for (i in a1){
    
    mapp[[init+by]]<-mapp[[init+by]]
    
    N<-length(mapp[[i-ov]]$points_in_vertex)
    N1<-length(mapp[[i]]$points_in_vertex)
    a<-matrix(list(),N,N1)
    link<-vector("list", max(N1,N))
    linksource<-vector("list", max(N1,N))
    linktarget<-vector("list", max(N1,N))
    linkvalue<-vector("list", max(N1,N))
    
    for (j in 1:N){
      for (k in 1:N1){
        #if (dim(intersect_all(mapp[[i-ov]]$points_in_vertex[[j]],
        #                      mapp[[i]]$points_in_vertex[[k]]))[1]>0 ){
        ind<-intersect_all(mapp[[i-ov]]$points_in_vertex[[j]], mapp[[i]]$points_in_vertex[[k]])
        ind1<-c(mapp[[i-ov]]$points_in_vertex[[j]], mapp[[i]]$points_in_vertex[[k]])
        adja00<-abs(cor(t(database[ind1,])))
        diag(adja00)<-0
        delta<-mean(apply(adja00,1,mean))
        indmax<-max(length(mapp[[i-ov]]$points_in_vertex[[j]]), length(mapp[[i]]$points_in_vertex[[k]]))
        
        if (dim(ind)[1]>=pero*indmax | delta>=corm){
          a[[j,k]]<-1
          #dim(a)==c(N,N1)
        } else {a[[j,k]]<-0}
        
        if (a[[j,k]]==1){a[[j,k]]=delta}
        
        #if (length(a[[j,k]])>0){ 
        if (a[[j,k]]!=0){ 
          if (i==init+by+ov){
            link[[j]][[k]]<-c(j+(i-(by+ov)-1),k+(i-by)-1)
          }
          if (i > init+by+ov){
            link[[j]][[k]]<-c(j+(i-(by+ov)-1),k+(i-by)-1)
          }
          linksource[[j]][[k]]<-link[[j]][[k]][1]
          linktarget[[j]][[k]]<-link[[j]][[k]][2]
          linkvalue[[j]][[k]]<-a[[j,k]]
        }
        #}
      }	
    }
    
    link1<-na.omit(cbind(unlist(linksource),unlist(linktarget),unlist(linkvalue)))
    
    if (length(links[[i]])>0){
      Linksource1<-as.numeric(links[[i]]$Linksource) + (i-by)-1
      Linktarget1<-as.numeric(links[[i]]$Linktarget) + (i-by)-1
      Linkvalue1<-as.numeric(links[[i]]$Linkvalue)
      links2[[i]]<-cbind(Linksource1,Linktarget1,Linkvalue1)
      
      if (length(link1)>0){
        colnames(link1)<-c("Linksource1","Linktarget1","Linkvalue1")
        
        link1<-as.data.frame(link1)
        
        links2[[i]]<-rbind(links2[[i]],link1)
      } else {links2[[i]]<-links2[[i]]}
      
      graphic1[[i]]<-graph.data.frame(links2[[i]], directed=FALSE, vertices=NULL)
      
    }
  }
  
  links2[[init+by]]<-links[[init+by]]
  colnames(links2[[init+by]])<-c("Linksource1","Linktarget1","Linkvalue1")
  
  #join all links
  for (i in a1){
    links2[[i]]=rbind(links2[[i]],links2[[i-ov]])
  }
  
  for (i in a1){
    graphic1[[i]]<-graph.data.frame(links2[[i]], directed=FALSE, vertices=NULL)
  }
  
  #join all graphs
  #for (i in a1){
  #  graphic1[[i]]=graph.union(graphic1[[i]],graphic1[[i-ov]])
  #}
})

##################################################################
#Identify CpGs per node and assign a categorical value to each one 

difnode <- defmacro(map,graphic2, database, subj, expr = {
  
  #select differentiated nodes case/control
  orig3_all<-as.numeric(names(igraph::degree(graphic2)))
  
  #points within those nodes
  p_all<-vector("list", maxi)
  a<-seq(init+by+ov,maxi,ov)
  for (i in a){
    N<-length(map[[i]]$points_in_vertex)
    for (l in 1:N){
      #h=i+l
      h=l+(i-by)-1
      p_all[[h]]<-map[[i]]$points_in_vertex[[l]]
    }}
  
  N<-length(map[[init+by]]$points_in_vertex)
  for (l in 1:N){
    p_all[[l]]<-map[[init+by]]$points_in_vertex[[l]]
  }
  
  points<-list()
  points_size<-as.numeric()
  
  points<-p_all[orig3_all]
  length(points)
  
  n<-rownames(database)
  n
  length(n)
  n0<-unlist(points)
  n1<-sort(n0)
  length(n1)
  
  ###############################
  #Methylation within each node
  
  n2_all<-list()
  for (i in 1:length(points)){
    n2_all[[i]]<-n[points[[i]]]
  }
  length(n2_all)
  
  nod_met<-list()
  nod_met_case<-list()
  nod_met_control<-list()
  nod_metcase<-list()
  nod_metcontrol<-list()
  nod_dif0<-list()
  
  nod_met1<-numeric()
  nod_met1_case<-numeric()
  nod_met1_control<-numeric()
  nod_dif<-numeric()
  
  nod_max<-numeric()
  nod_max_case<-numeric()
  nod_max_control<-numeric()
  nod_max1<-numeric()
  nod_max1_case<-numeric()
  nod_max1_control<-numeric()
  
  cor_met<-numeric()
  cor_met_case<-numeric()
  cor_met_control<-numeric()
  cor_met_med<-numeric()
  chr<-numeric()
  size<-numeric()
  nod_var<-numeric()
  
  for (i in 1:length(n2_all)){
    
    database1<-na.omit(database[n2_all[[i]],subj])
    #annotation_chr<-annotation[n2_all[[i]],1]
    
    if (length(n2_all[[i]])>1) {
      cord<-abs(cor(t(database1)))
      diag(cord)<-0
      nod_met[[i]]<-apply(database1,1,median)
      nod_met1[i]<-median(nod_met[[i]]) #max

      nod_var[[i]]<-median(apply(database1,1,var))
      
      if (casenet==1 & controlnet==1){
      nod_metcase[[i]]<-apply(database[n2_all[[i]],rcase],1,median)
      nod_metcontrol[[i]]<-apply(database[n2_all[[i]],rcontrol],1,median)
      nod_dif0[[i]]<- abs(nod_metcase[[i]]-nod_metcontrol[[i]])
      nod_dif[i]<- median(nod_dif0[[i]]) #max
      }
      nod_max[i]<-max(database1)
      cor_met[i]<-max(cord)
      cor_met_med[i]<-median(cord)
      #chr[i]<-length(unique(annotation_chr))
      size[i]<-length(n2_all[[i]])
    }
    else if (length(n2_all[[i]])==1) {
      cord<-abs(cor(t(database1)))
      nod_met[[i]]<-median(as.numeric(database1))
      nod_var[[i]]<-var(as.numeric(database1))
      
      if (casenet==1 & controlnet==1){
      nod_metcase[[i]]<-median(as.numeric(database[n2_all[[i]],rcase]))
      nod_metcontrol[[i]]<-median(as.numeric(database[n2_all[[i]],rcontrol]))
      nod_dif[i]<- abs(nod_metcase[[i]]-nod_metcontrol[[i]])
      }
      
      nod_met1[i]<-nod_met[[i]]
      nod_max[i]<-max(database1)
      cor_met[i]<-1
      cor_met_med[i]<-1
      #chr[i]<-length(unique(annotation_chr))
      size[i]<-1
    }
    else if (length(n2_all[[i]])==0) {
      nod_met[[i]]<-0
      nod_met1[i]<-0
      nod_max[i]<-0
      nod_var[[i]]<-0
      cor_met[i]<-0
      cor_met_med[i]<-0
      #chr[i]<-0
      if (casenet==1 & controlnet==1){
      nod_dif[i]<-0
      }
      size[i]<-0
    }
  }	
  
  nod_met1<-nod_met1
  length(nod_met1)
  length(orig3_all)
  names(nod_met1)<-as.character(orig3_all)
  
  nod_var<-nod_var
  length(nod_var)
  length(orig3_all)
  names(nod_var)<-as.character(orig3_all)
  
  nod_max<-nod_max
  length(nod_max)
  length(orig3_all)
  names(nod_max)<-as.character(orig3_all)
  
  if (casenet==1 & controlnet==1){
  nod_dif<-nod_dif
  length(nod_dif)
  length(orig3_all)
  names(nod_dif)<-as.character(orig3_all)
  }
  
  cor_met<-cor_met[abs(cor_met)>0]
  length(cor_met)
  length(orig3_all)
  names(cor_met)<-as.character(orig3_all)
  
  cor_met_med<-cor_met_med[abs(cor_met_med)>0]
  length(cor_met_med)
  length(orig3_all)
  names(cor_met_med)<-as.character(orig3_all)
  
  #names(size)<-as.character(orig3_all)
  #names(chr)<-as.character(orig3_all)
  
  exp <- nod_met1

  # color palette
  exp1<-as.numeric()
  max1<-as.numeric()
  dif1<-as.numeric()
  var1<-as.numeric()
  
  #exp_beta<-m2beta(exp)

  for (i in 1:length(exp)){
    if (exp[i]>quantile(exp)[4]) {exp1[i]="1"
    } else if (exp[i]>=quantile(exp)[2] & exp[i]<=quantile(exp)[4]) {exp1[i]="2"
    } else if (exp[i]>=0 & exp[i]<quantile(exp)[2]) {exp1[i]="3"}
  }
  names(exp1)<-names(exp)
  # for (i in 1:length(exp)){
  #   if (exp[i]>0.8) {exp1[i]="1"
  #   } else if (exp[i]>=0.6 & exp[i]<=0.8) {exp1[i]="2"
  #   } else if (exp[i]>=0.4 & exp[i]<0.6) {exp1[i]="3"
  #   } else if (exp[i]>=0.2 & exp[i]<0.4) {exp1[i]="4"
  #   } else if (exp[i]>=0 & exp[i]<0.2) {exp1[i]="5"}
  # }
  # names(exp1)<-names(exp)
  # 
  for (i in 1:length(nod_var)){
    if (nod_var[i]>quantile(nod_var)[4]) {var1[i]="1"
    } else if (nod_var[i]>=quantile(nod_var)[2] & nod_var[i]<=quantile(nod_var)[4]) {var1[i]="2"
    } else if (nod_var[i]>=0 & nod_var[i]<quantile(nod_var)[2]) {var1[i]="3"}
  }
  names(var1)<-names(nod_var)
  # for (i in 1:length(nod_var)){
  #   if (nod_var[i]>0.8) {var1[i]="1"
  #   } else if (nod_var[i]>=0.6 & nod_var[i]<=0.8) {var1[i]="2"
  #   } else if (nod_var[i]>=0.4 & nod_var[i]<0.6) {var1[i]="3"
  #   } else if (nod_var[i]>=0.2 & nod_var[i]<0.4) {var1[i]="4"
  #   } else if (nod_var[i]>=0 & nod_var[i]<0.2) {var1[i]="5"}
  # }
  # names(var1)<-names(nod_var)
  # 
  for (i in 1:length(nod_max)){
    if (nod_max[i]>0.8) {max1[i]="1"
    } else if (nod_max[i]>=0.5 & nod_max[i]<=0.8) {max1[i]="2"
    } else if (nod_max[i]>=0 & nod_max[i]<0.5) {max1[i]="3"}
  }
  
  names(max1)<-names(nod_max)
  
  if (casenet==1 & controlnet==1){
  for (i in 1:length(nod_dif)){
    if (nod_dif[i]>quantile(nod_dif)[4]) {dif1[i]="1"
    } else if (nod_dif[i]>=quantile(nod_dif)[2] & nod_dif[i]<=quantile(nod_dif)[4]) {dif1[i]="2"
    } else if (nod_dif[i]>=0 & nod_dif[i]<quantile(nod_dif)[2]) {dif1[i]="3"}
  }

  names(dif1)<-names(nod_dif)
  }
  
  # for (i in 1:length(nod_dif)){
  #   if (nod_dif[i]>0.8) {dif1[i]="1"
  #   } else if (nod_dif[i]>=0.6 & nod_dif[i]<=0.8) {dif1[i]="2"
  #   } else if (nod_dif[i]>=0.4 & nod_dif[i]<0.6) {dif1[i]="3"
  #   } else if (nod_dif[i]>=0.2 & nod_dif[i]<0.4) {dif1[i]="4"
  #  } else if (nod_dif[i]>=0 & nod_dif[i]<0.2) {dif1[i]="5"}
  # }
  # 
  #names(dif1)<-names(nod_dif)
  
  ####################
  #correlation
  
  cor1<-numeric()
  cor1_med<-numeric()
  
  for (i in 1:length(cor_met)){
    if (cor_met[i]>0.8) {cor1[i]="1"
    } else if (cor_met[i]>=0.5 & cor_met[i]<=0.8) {cor1[i]="2"
    } else if (cor_met[i]>=0 & cor_met[i]<0.5) {cor1[i]="3"
    }
  }
  for (i in 1:length(cor_met_med)){
    if (cor_met_med[i]>quantile(cor_met_med)[4]) {cor1_med[i]="1"
    } else if (cor_met_med[i]>=quantile(cor_met_med)[2] & cor_met_med[i]<=quantile(cor_met_med)[4]) {cor1_med[i]="2"
    } else if (cor_met_med[i]>=0 & cor_met_med[i]<quantile(cor_met_med)[2]) {cor1_med[i]="3"
    }
  }
  # for (i in 1:length(cor_met_med)){
  #   if (cor_met_med[i]>0.8) {cor1_med[i]="1"
  #   } else if (cor_met_med[i]>=0.6 & cor_met_med[i]<=0.8) {cor1_med[i]="2"
  #   } else if (cor_met_med[i]>=0.4 & cor_met_med[i]<0.6) {cor1_med[i]="3"
  #   } else if (cor_met_med[i]>=0.2 & cor_met_med[i]<0.4) {cor1_med[i]="4"
  #   } else if (cor_met_med[i]>=0 & cor_met_med[i]<0.2) {cor1_med[i]="5"}
  # }
  names(cor1)<-names(cor_met)
  names(cor1_med)<-names(cor_met_med)
  
  # chr1<-numeric()
  # 
  # for (i in 1:length(chr)){
  #   if (chr[i]>10) {chr1[i]="1"
  #   } else if (chr[i]>=10 & chr[i]<=2) {chr1[i]="2"
  #   } else if (chr[i]>=0 & chr[i]<2) {chr1[i]="3"
  #   }
  # }
  
  # names(chr1)<-names(chr)
  
  # size1<-numeric()
  # 
  # for (i in 1:length(size)){
  #   if (size[i]>100) {size1[i]="1"
  #   } else if (size[i]>=10 & size[i]<=100) {size1[i]="2"
  #   } else if (size[i]>=0 & size[i]<10) {size1[i]="3"
  #   }
  # }
  # 
  # names(size1)<-names(size)
  #plot different cases
  
  # plot
  #reset_par()
  #png("Network.png", width = 4, height = 4, units = 'in', res = 600)
  #par(bg="black", mar=c(0,0,0,0))
  #plot(graphic2, 
  #     vertex.size=2,
  #     vertex.color="green", 
  #     edge.color="blue",
  #     vertex.label=NA,
  #     vertex.label.cex=0.7,
  #     vertex.label.color="white",
  #     vertex.frame.color="transparent")
  #title(main="Network",sub=text1,cex.main=1,font.main= 4,col.main="white",col.sub = "white")
  #dev.off()
})

##################################
#Identify differentiated CpG sites

difcg <- defmacro(vect,cat,pointnode,database,subj,pngname,pngname1,pngname2, expr = {
  
  #vertex differenciated
  exp2<-vect[vect %in% cat]
  
  if (length(exp2)>0){
    
    exp21<-as.numeric(names(exp2))
    
    points<-list()
    points<-pointnode[exp21]
    length(points)
    
    n<-rownames(database)
    
    n2<-list()
    for (i in 1:length(points)){
      n2[[i]]<-n[points[[i]]]
    }
    length(n2)
    
    n3<-as.numeric()
    n3<-unique(unlist(n2))
    length(n3)
    
    #if (length(n3)<1500){
    #cor_n3<-cor(t(database[n3,]))
    #dim(cor_n3)
    
    #reset_par()
    #Heatmap(cor_n3, name = "correlation")
    
    chr1<-as.numeric(substr(annotation[n3,]$chr,4,6))
    sort(chr1)
    annotation_chr<-cbind(annotation[n3,],chr1)
    chr1<-unique(sort(chr1))
    annotation_chr<-annotation_chr[order(annotation_chr$chr1,annotation_chr$pos),]
    head(annotation_chr)
    
    posir<-rownames(annotation_chr)
    
    #hyper/hypo methylated region
    # state_case<-numeric()
    # mean_posi<-apply(database[posir,rcase],1,median)
    # for (i in 1:length(posir)){
    #   if (mean_posi[i]>0.5){state_case[i]<-"Hyper"
    #   } else {state_case[i]<-"Hypo"}
    # }
    # 
    # state_control<-numeric()
    # mean_posi<-apply(database[posir,rcontrol],1,median)
    # for (i in 1:length(posir)){
    #   if (mean_posi[i]>0.5){state_control[i]<-"Hyper"
    #   } else {state_control[i]<-"Hypo"}
    # }
    #names(state)<-posir
    #ordered by position
    #cor_n31<-cor_n3[names(posi),names(posi)]
    #head(cor_n31)
    
    #all_31<-database[names(posi),]
    
    #methylation
    
    #col_fun = colorRamp2(c(0, 0.5, 1), c("white","grey","black"))
    col_fun <- colorRampPalette(brewer.pal(10, "RdYlBu"))(256)
    
    # Create the heatmap annotation
    ha <- HeatmapAnnotation(Status = all.meta1[subj,]$disease,col = col1)
    
    h1<-Heatmap(as.matrix(database[posir,subj]), name = "Methylation",clustering_distance_rows = "pearson",show_row_names=FALSE,show_column_names=FALSE,show_column_dend = TRUE,row_dend_reorder = TRUE,column_title = "",col = col_fun,column_km = (length(col1$Status)+1), show_parent_dend_line = FALSE, top_annotation = ha) +
      Heatmap(annotation_chr[posir,]$chr1, name = "Chromosome", width = unit(5, "mm"),
              col =  rainbow(length(unique(annotation_chr[posir,]$chr1))))
    reset_par()
    png(pngname, width = 10, height =10, units = 'in', res = 300)
    print(h1)
    dev.off()
    
    # #ordered by chromosome
    # ha <- HeatmapAnnotation(Status = all.meta1[subj,]$disease,col = col1)
    # 
    # h2<-Heatmap(as.matrix(database[posir,subj]), name = "Methylation",row_order=posir,show_row_names=FALSE,show_column_names=FALSE,column_title = "",col = col_fun, top_annotation = ha) +
    #   Heatmap(state_case, name = "State Case", width = unit(5, "mm"),
    #           col =  rainbow(length(unique(state_case))))+
    #   Heatmap(state_control, name = "State Control", width = unit(5, "mm"),
    #           col =  rainbow(length(unique(state_control))))+
    #   
    #   reset_par()
    # png(pngname1, width = 10, height =10, units = 'in', res = 300)
    # print(h2)
    # dev.off()
    # 
    #ordered by sample
    ha <- HeatmapAnnotation(Status = all.meta1[subj,]$disease,col = col1)
    
    h3<-Heatmap(as.matrix(database[posir,subj]), name = "Methylation",column_order=subj,show_row_names=FALSE,show_column_names=FALSE,column_title = "",col = col_fun, top_annotation = ha) +
      Heatmap(annotation_chr[posir,]$chr1, name = "Chromosome", width = unit(5, "mm"),
              col =  rainbow(length(unique(annotation_chr[posir,]$chr1))))
    reset_par()
    png(pngname2, width = 10, height =10, units = 'in', res = 300)
    print(h3)
    dev.off()
    
    #correlation
    #h2<-Heatmap(cor_n31, name = "correlation",row_order=rownames(cor_n31),column_order=colnames(cor_n31),show_row_names=FALSE,show_column_names=FALSE) +
    #  Heatmap(annotation_chr[rownames(cor_n31),]$chr1, name = "chr", width = unit(5, "mm"),
    #          col =  rainbow(length(annotation_chr[rownames(cor_n31),]$chr1)))
    #reset_par()
    #png("heatmap_cor.png", width = 10, height =10, units = 'in', res = 300)
    #print(h2)
    #dev.off()
    
    #correlation importance
    #png("importance_cor.png", width = 10, height =10, units = 'in', res = 300)
    #corr_cross(cor_n3, # name of dataset
    #           max_pvalue = 0.05, # display only significant correlations (at 5% level)
    #           top = 50)
    #dev.off()
    
    #correlation test
    #library(Hmisc)
    #res2<-rcorr(as.matrix(cor_n31))
    #flattenCorrMatrix(res2$r, res2$P)
    
    #library("PerformanceAnalytics")
    #my_data <- cor_n31[1:100,1:100]
    #chart.Correlation(my_data, histogram=TRUE, pch=19)
    #}
    
    anotcor<-annotation[n3,]
    print("CpGs characteristics")
    print(length(n3))
    print(sort(table(anotcor$chr)/length(n3))*100)
    
    #plot genes frequency
    #reset_par()
    #par(mfrow=c(2,1))
    barcor<-bargen(n3,database,0)
    outgenes<-barcor[[1]]
    outcounts<-barcor[[2]]
    genescgs<-barcor[[3]]
    
    #characteristics of the CpGs
    islands<-anotcor$Relation_to_Island
    print(sort(table(islands)/length(n3))*100)
    
    regulatory<-anotcor$Regulatory_Feature_Group
    print(sort(table(regulatory)/length(n3))*100)
    
    #group<-anotcor$UCSC_RefGene_Group
    #print(sort(table(group)/length(n3))*100)
    
  }
})

#######################
#Plot colored functions

plotmap<-function(graphic2,data,text,cat,text1){
  
  coul <- rainbow(3)[1:nlevels(as.factor(data))]
  #brewer.pal(nlevels(as.factor(data)), "Set2")
  # Map the color to cylinders
  my_color <- coul[as.numeric((data))]
  names(my_color)<-names(data)
  
  edge.start <- ends(graphic2, es=E(graphic2), names=T)[,1]
  edge.col <- data[edge.start]
  
  #coordinates
  #coords <- layout.davidson.harel(graphic2)
  
  # plot
  reset_par()
  png(text1, width = 6, height =6, units = 'in', res = 600)
  par(bg="black")
  set.seed(4)
  
  data_aux<-as.numeric(data)
  names(data_aux)<-names(data)
  
  mst.plot.mod(graphic2,bg="black", v.size=2,e.size=.20,
               mst.e.size=1.2, sf=0,mst.edge.col="white",colors= "white", vertex.color = "green",
               expression=data_aux,v.lab=FALSE,v.lab.cex=0.7,v.lab.col="blue",edge.col.wt=as.numeric(edge.col))
  
  legend(x=0.1 , y=1, 
         legend=paste(cat, sep=""), 
         col = c("red" ,"deeppink", "blue") , 
         bty = "n", pch=20 , pt.cex = 2, cex =0.9,
         text.col="white" , horiz = F,title =text)
  dev.off()
  
  #plot(graphic2, 
  #            vertex.size=2,
  #            vertex.color=my_color, 
  #            edge.color=edge.col,
  #            vertex.label=NA,
  #            vertex.label.cex=0.7,
  #            vertex.label.color="white",
  #            vertex.frame.color="transparent",layout=layout.graphopt)
  # title and legend
  #legend(x=0.1 , y=1, 
  #       legend=paste(cat, sep=""), 
  #       col = rainbow(3) , 
  #       bty = "n", pch=20 , pt.cex = 2, cex =0.7,
  #       text.col="white" , horiz = F,title =text)
  
  #data_aux<-as.numeric(data)
  #names(data_aux)<-names(data)
  
  #mst.plot.mod(graphic2,bg="black", v.size=2,e.size=.25,
  #             mst.e.size=1.2, sf=0,mst.edge.col="white",colors= "white",
  #             expression=data_aux,vertex.color =coul)
  
}

####################################################
#Manhattan plot of methylation values per chromosome

manhattan_methyl<-function (x, chr = "CHR", bp = "BP", p = "P",
                            snp = "SNP", col = c("gray10", "gray60"), 
                            chrlabs = NULL, suggestiveline = -log10(1e-05), genomewideline = -log10(5e-08), 
                            highlight = NULL, logp = TRUE, annotatePvalmax = NULL, annotatePvalmin = NULL,annotateTop = TRUE,...)
{
  CHR = BP = P = index = NULL
  if (!(chr %in% names(x))) 
    stop(paste("Column", chr, "not found!"))
  if (!(bp %in% names(x))) 
    stop(paste("Column", bp, "not found!"))
  if (!(p %in% names(x))) 
    stop(paste("Column", p, "not found!"))
  if (!(snp %in% names(x))) 
    warning(paste("No SNP column found. OK unless you're trying to highlight."))
  if (!is.numeric(x[[chr]])) 
    stop(paste(chr, "column should be numeric. Do you have 'X', 'Y', 'MT', etc? If so change to numbers and try again."))
  if (!is.numeric(x[[bp]])) 
    stop(paste(bp, "column should be numeric."))
  if (!is.numeric(x[[p]])) 
    stop(paste(p, "column should be numeric."))
  d = data.frame(CHR = x[[chr]], BP = x[[bp]], P = x[[p]])
  if (!is.null(x[[snp]])) 
    d = transform(d, SNP = x[[snp]])
  d <- subset(d, (is.numeric(CHR) & is.numeric(BP) & is.numeric(P)))
  d <- d[order(d$CHR, d$BP), ]
  if (logp) {
    d$logp <- -log10(d$P)
  }
  else {
    d$logp <- d$P
  }
  
  d$pos = NA
  d$index = NA
  ind = 0
  for (i in unique(d$CHR)) {
    ind = ind + 1
    d[d$CHR == i, ]$index = ind
  }
  nchr = length(unique(d$CHR))
  if (nchr == 1) {
    d$pos = d$BP
    ticks = floor(length(d$pos))/2 + 1
    xlabel = paste("Chromosome", unique(d$CHR), "position")
    labs = ticks
    ylabel= "Methylation"
  }
  else {
    lastbase = 0
    ticks = NULL
    for (i in unique(d$index)) {
      if (i == 1) {
        d[d$index == i, ]$pos = d[d$index == i, ]$BP
      }
      else {
        lastbase = lastbase + tail(subset(d, index == 
                                            i - 1)$BP, 1)
        d[d$index == i, ]$pos = d[d$index == i, ]$BP + 
          lastbase
      }
      ticks = c(ticks, (min(d[d$index == i, ]$pos) + max(d[d$index == 
                                                             i, ]$pos))/2 + 1)
    }
    xlabel = "Chromosome"
    labs <- unique(d$CHR)
    ylabel= "Methylation"
  } 
  
  xmax = ceiling(max(d$pos) * 1.03)
  xmin = floor(max(d$pos) * -0.03)
  def_args <- list(xaxt = "n", bty = "n", xaxs = "i", 
                   yaxs = "i", las = 1, pch = 20, xlim = c(xmin, xmax), 
                   ylim = c(floor(min(d$logp)), ceiling(max(d$logp))), xlab = xlabel, ylab =ylabel)
  dotargs <- list(...)
  do.call("plot", c(NA, dotargs, def_args[!names(def_args) %in% 
                                            names(dotargs)]))
  if (!is.null(chrlabs)) {
    if (is.character(chrlabs)) {
      if (length(chrlabs) == length(labs)) {
        labs <- chrlabs
      }
      else {
        warning("You're trying to specify chromosome labels but the number of labels != number of chromosomes.")
      }
    }
    else {
      warning("If you're trying to specify chromosome labels, chrlabs must be a character vector")
    }
  }
  if (nchr == 1) {
    axis(1, ...)
  }
  else {
    axis(1, at = ticks, labels = labs, ...)
  }
  col = rep(col, max(d$CHR))
  if (nchr == 1) {
    with(d, points(pos, logp, pch = 20, col = col[1], ...))
  }
  else {
    icol = 1
    for (i in unique(d$index)) {
      with(d[d$index == unique(d$index)[i], ], points(pos, 
                                                      logp, col = col[icol], pch = 20, ...))
      icol = icol + 1
    }
  }
  if (suggestiveline) 
    abline(h = suggestiveline, col = "blue")
  if (genomewideline) 
    abline(h = genomewideline, col = "red")
  if (!is.null(highlight)) {
    if (any(!(highlight %in% d$SNP))) 
      warning("You're trying to highlight SNPs that don't exist in your results.")
    d.highlight = d[which(d$SNP %in% highlight), ]
    with(d.highlight, points(pos, logp, col = "green3", 
                             pch = 20, ...))
  }
  if (!is.null(annotatePvalmin)) {
    topHits = subset(d, P <= annotatePvalmin)
    par(xpd = TRUE)
    if (annotateTop == FALSE) {
      with(subset(d, P<= annotatePvalmin), textxy(pos, P, 
                                                  offset = 0.625, labs = topHits$SNP, cex = 0.45), 
           ...)
    }
    else { with(subset(d, P<= annotatePvalmin), textxy(pos, P, 
                                                       offset = 0.625, labs =topHits$SNP, cex = 0.45), 
                ...)
      #topHits <- topHits[order(topHits$P), ]
      #topSNPs <- NULL
      #for (i in unique(topHits$CHR)) {
      #  chrSNPs <- topHits[topHits$CHR == i, ]
      #  topSNPs <- rbind(topSNPs, chrSNPs[1, ])
      #}
      #textxy(topSNPs$pos, (topSNPs$P), offset = 0.625, 
      #       labs = topSNPs$SNP, cex = 0.5, ...)
    }
  }
  if (!is.null(annotatePvalmax)) {
    topHits = subset(d, P >= annotatePvalmax)
    par(xpd = TRUE)
    if (annotateTop == FALSE) {
      with(subset(d, P>= annotatePvalmax), textxy(pos, P, 
                                                  offset = 0.625, labs = topHits$SNP, cex = 0.45), 
           ...)
    }
    else { with(subset(d, P>= annotatePvalmax), textxy(pos, P, 
                                                       offset = 0.625, labs =topHits$SNP, cex = 0.45), 
                ...)
      #topHits <- topHits[order(topHits$P), ]
      #topSNPs <- NULL
      #for (i in unique(topHits$CHR)) {
      #  chrSNPs <- topHits[topHits$CHR == i, ]
      #  topSNPs <- rbind(topSNPs, chrSNPs[1, ])
      #}
      #textxy(topSNPs$pos, (topSNPs$P), offset = 0.625, 
      #       labs = topSNPs$SNP, cex = 0.5, ...)
    }
  }
  par(xpd = FALSE)
}


##################################
#Plot genes from the selected CpGs

bargen <- function (cg,data,num) 
{
  
  #genes associated to the site
  gensig<-annotation0[cg,]
  gensig
  gene1 <- gensig$UCSC_RefGene_Name
  
  gene1_d <- as.character(gene1)
  names(gene1_d)<-rownames(gensig)
  
  if(length(gene1_d)>=1){
    gene1_d10<-strsplit(gene1_d, c(";"))
    gene1_d1<-unlist(gene1_d10)
    gene1_d11<-data.frame(ID = rep(names(gene1_d10), sapply(gene1_d10, length)),Obs = unlist(gene1_d10))
    gene1_d11<-as.data.frame(gene1_d11)
    gene1_d12<-as.character(gene1_d11$Obs)
    names(gene1_d12)<-gene1_d11$ID
    #gene1_d1 <- as.data.frame(gsub(";.*$", "", gene1_d))
    #colnames(gene1_d1)='gene'
    counts_gen <- table((gene1_d12))
    counts_gen_name <- rownames(counts_gen)
    
    #relative percentage based on the total CGs related to that gene 
    gensigall<-annotation0[rownames(data),]
    gensigall
    gene1all <- gensigall$UCSC_RefGene_Name
    gene1_dall<- as.character(gene1all)
    
    gene1_d10all<-strsplit(gene1_dall, ";")
    gene1_d1all<-unlist(gene1_d10all)
    
    counts_genall <- table(gene1_d1all)
    counts_gen_nameall <- rownames(counts_genall)
    
    counts_gen_num<-as.data.frame(counts_gen)
    counts_gen_numall<-as.data.frame(counts_genall)
    colnames(counts_gen_numall)<-c("Var1", "Freqall")
    colnames(counts_gen_num)<-c("Var1", "Freq")
    counts_all<-merge(counts_gen_num, counts_gen_numall,by="Var1")
    head(counts_all)
    
    #percentage based on the genes associated to the total CGs selected
    counts_all$Freqall1<-length(cg)
    
    attach(counts_all)
    counts_all$counts_per<-(Freq/Freqall)*100
    detach(counts_all)
    
    graphics::barplot(counts_all$counts_per, main="Gene distribution",names.arg=as.character(counts_all$Var1),
                      col=c(hue_pal()(dim(counts_all)[1])),  cex.names = 0.65, cex.lab= 0.65,las=2)
    
    #genes more present
    counts_sig<- counts_all[counts_all$counts_per>max(counts_all$counts_per)/2,]
    graphics::barplot(counts_sig$counts_per, main="Gene distribution - Genes more present",names.arg=as.character(counts_sig$Var1),
                      col=c(hue_pal()(dim(counts_sig)[1])),  cex.names = 0.65, cex.lab= 0.65,las=2)
    
    gene1_d123<-gene1_d12[counts_all$Var1]
    out <- list(counts_all$Var1, counts_all$counts_per,gene1_d12)  
    return(out)
  }
}

########################################
#Plot site region from the selected CpGs

sitegen <- function (cg1,cg2,data,num) 
{
  
  #genes associated to the CG
  gensig<-annotation[cg1,]
  gensig
  gene1 <- gensig$UCSC_RefGene_Group
  
  gene1_d <- as.character(gene1)
  names(gene1_d)<-rownames(gensig)
  
  if(length(gene1_d)>=1){
    gene1_d10<-strsplit(gene1_d, c(";"))
    gene1_d1<-unlist(gene1_d10)
    gene1_d11<-data.frame(ID = rep(names(gene1_d10), sapply(gene1_d10, length)),Obs = unlist(gene1_d10))
    gene1_d11<-as.data.frame(gene1_d11)
    gene1_d12<-as.character(gene1_d11$Obs)
    names(gene1_d12)<-gene1_d11$ID
    #gene1_d1 <- as.data.frame(gsub(";.*$", "", gene1_d))
    #colnames(gene1_d1)='gene'
    counts_gen <- table((gene1_d12))
    counts_gen_name <- rownames(counts_gen)
    
    gensigall<-annotation[cg2,]
    gensigall
    gene1all <- gensigall$UCSC_RefGene_Group
    gene1_dall<- as.character(gene1all)
    
    gene1_d10all<-strsplit(gene1_dall, ";")
    gene1_d1all<-unlist(gene1_d10all)
    
    counts_genall <- table(gene1_d1all)
    counts_gen_nameall <- rownames(counts_genall)
    
    counts_gen_num<-as.data.frame(counts_gen)
    counts_gen_numall<-as.data.frame(counts_genall)
    colnames(counts_gen_numall)<-c("Var1", "Freq")
    counts_gen_numall$category="Hypo"
    counts_gen_numall$Freqall1<-sum(counts_gen_numall$Freq)
    
    colnames(counts_gen_num)<-c("Var1", "Freq")
    counts_gen_num$category="Hyper"
    counts_gen_num$Freqall1<-sum(counts_gen_num$Freq)
    
    counts_all<-rbind(counts_gen_num, counts_gen_numall)
    head(counts_all)
    
    attach(counts_all)
    counts_all$counts_per<-round((Freq/Freqall1)*100,1)
    detach(counts_all)
    
    col_fun <- colorRampPalette(brewer.pal(10, "RdYlBu"))(2)
    
    gplot<-ggplot(data=counts_all, aes(x=Var1, y=counts_per, fill=category)) +
      geom_bar(stat="identity", position=position_dodge())+
      geom_text(aes(label=counts_per), vjust=1.6, color="black",
                position = position_dodge(0.9), size=4)+
      scale_fill_manual(values=c("#74ADD1","#F46D43"))+xlab("Genomic region")+ylab("Percentage")+
      theme(text = element_text(size=15)) 
    print(gplot)

    gene1_d123<-gene1_d12[counts_all$Var1]
    out <- list(counts_all$Var1, counts_all$counts_per,gene1_d12)  
    return(out)
  }
}

##################################################
#Plot relationship to islands of the selected CpGs

islandgen <- function (cg1,cg2,data,num) 
{
  
  #genes associated to the CG
  gensig<-annotation[cg1,]
  gensig
  gene1 <- gensig$Relation_to_Island
  
  gene1_d <- as.character(gene1)
  names(gene1_d)<-rownames(gensig)
  
  if(length(gene1_d)>=1){
    gene1_d10<-strsplit(gene1_d, c(";"))
    gene1_d1<-unlist(gene1_d10)
    gene1_d11<-data.frame(ID = rep(names(gene1_d10), sapply(gene1_d10, length)),Obs = unlist(gene1_d10))
    gene1_d11<-as.data.frame(gene1_d11)
    gene1_d12<-as.character(gene1_d11$Obs)
    names(gene1_d12)<-gene1_d11$ID
    #gene1_d1 <- as.data.frame(gsub(";.*$", "", gene1_d))
    #colnames(gene1_d1)='gene'
    counts_gen <- table((gene1_d12))
    counts_gen_name <- rownames(counts_gen)
    
    gensigall<-annotation[cg2,]
    gensigall
    gene1all <- gensigall$Relation_to_Island
    gene1_dall<- as.character(gene1all)
    
    gene1_d10all<-strsplit(gene1_dall, ";")
    gene1_d1all<-unlist(gene1_d10all)
    
    counts_genall <- table(gene1_d1all)
    counts_gen_nameall <- rownames(counts_genall)
    
    counts_gen_num<-as.data.frame(counts_gen)
    counts_gen_numall<-as.data.frame(counts_genall)
    colnames(counts_gen_numall)<-c("Var1", "Freq")
    counts_gen_numall$category="Hypo"
    counts_gen_numall$Freqall1<-sum(counts_gen_numall$Freq)
    
    colnames(counts_gen_num)<-c("Var1", "Freq")
    counts_gen_num$category="Hyper"
    counts_gen_num$Freqall1<-sum(counts_gen_num$Freq)
    
    counts_all<-rbind(counts_gen_num, counts_gen_numall)
    head(counts_all)
    
    attach(counts_all)
    counts_all$counts_per<-round((Freq/Freqall1)*100,1)
    detach(counts_all)
    
    col_fun <- colorRampPalette(brewer.pal(10, "RdYlBu"))(2)
    
    gplot<-ggplot(data=counts_all, aes(x=Var1, y=counts_per, fill=category)) +
      geom_bar(stat="identity", position=position_dodge())+
      geom_text(aes(label=counts_per), vjust=1.6, color="black",
                position = position_dodge(0.9), size=4)+
      scale_fill_manual(values=c("#74ADD1","#F46D43"))+xlab("Island or Ocean")+ylab("Percentage")+
      theme(text = element_text(size=15)) 
    print(gplot)

    gene1_d123<-gene1_d12[counts_all$Var1]
    out <- list(counts_all$Var1, counts_all$counts_per,gene1_d12)  
    return(out)
  }
}

#######################
#Mapper-based functions

mapperVertices <- function(m, pt_labels) {
  
  # Hovering over vertices gives the point labels:
  # convert the list of vectors of point indices to a list of vectors of labels
  labels_in_vertex <- lapply( m$points_in_vertex, FUN=function(v){ pt_labels[v] } )
  nodename <- sapply( sapply(labels_in_vertex, as.character), paste0, collapse=", ")
  nodename <- paste0("V", 1:m$num_vertices, ": ", nodename )
  
  # Hovering over vertices gives the point indices:
  # list the points in each vertex
  # nodename <- sapply( sapply(m$points_in_vertex, as.character), paste0, collapse=", ")
  # concatenate the vertex number with the labels for the points in each vertex
  #nodename <- paste0("V", 1:m$num_vertices, ": ", nodename )
  
  nodegroup <- m$level_of_vertex
  nodesize <- sapply(m$points_in_vertex, length)
  
  return(data.frame( Nodename=nodename, 
                     Nodegroup=nodegroup, 
                     Nodesize=nodesize ))
  
}

mapperEdges_w<-function (m) {
  linksource <- c()
  linktarget <- c()
  linkvalue <- c()
  k <- 1
  for (i in 2:m$num_vertices) {
    for (j in 1:(i - 1)) {
      if (m$adjacency[i, j] > 0) {
        linksource[k] <- i
        linktarget[k] <- j
        linkvalue[k] <- m$adjacency[i, j]
        k <- k + 1
      }
    }
  }
  return(data.frame(Linksource = linksource, Linktarget = linktarget, 
                    Linkvalue = linkvalue))
}

# the inverse function
lsfi_from_lsmi <- function( lsmi, num_intervals ) {
  lsfi <- lsmi[1]
  if (length(num_intervals) > 1) {
    for (i in 2:length(num_intervals)) {
      lsfi <- lsfi + prod(num_intervals[1:(i-1)]) * (lsmi[i]-1)
    }
  }
  return(lsfi)
}
# function from the level set flat index (lsfi) to the level set multi-index (lsmi)
lsmi_from_lsfi <- function( lsfi, num_intervals ) {
  # inputs:
  # lsfi = an integer in the range 1:prod(v)
  # num_intervals = c(i1,i1,...) a vector of numbers of intervals
  # output:
  # f+1 = a vector of multiindices with length filter_output_dim
  j <- c(1,num_intervals) # put 1 in front to make indexing easier in the product prod(j[1:k])
  f <- c()
  for (k in 1:length(num_intervals)) {
    # use lsfi-1 to shift from 1-based indexing to 0-based indexing
    f[k] <- floor( (lsfi-1) / prod(j[1:k])) %% num_intervals[k]
  }
  #print(f+1)
  # lsmi = f+1 = level set multi index
  return(f+1) # shift from 0-based indexing back to 1-based indexing
}

cluster_cutoff_at_first_empty_bin <- function(heights, diam, num_bins_when_clustering) {
  
  # if there are only two points (one height value), then we have a single cluster
  if (length(heights) == 1) {
    if (heights == diam) {
      cutoff <- Inf
      return(cutoff)
    }
  }
  
  bin_breaks <- seq(from=min(heights), to=diam, 
                    by=(diam - min(heights))/num_bins_when_clustering)
  if (length(bin_breaks) == 1) { bin_breaks <- 1 }
  
  myhist <- hist(c(heights,diam), breaks=bin_breaks, plot=FALSE)
  z <- (myhist$counts == 0)
  if (sum(z) == 0) {
    cutoff <- Inf
    return(cutoff)
  } else {
    #  which returns the indices of the logical vector (z == TRUE), min gives the smallest index
    cutoff <- myhist$mids[ min(which(z == TRUE)) ]
    return(cutoff)
  }
  
}

mapper_corr_w <- function (object,dist_object, filter_values, num_intervals, percent_overlap, 
                           num_bins_when_clustering, meth,pero0,corm0) 
{
  filter_values <- data.frame(filter_values)
  num_points <- dim(filter_values)[1]
  filter_output_dim <- dim(filter_values)[2]
  num_intervals<-c(int,int)
  num_levelsets <- prod(num_intervals)
  filter_min <- as.vector(sapply(filter_values, min))
  filter_max <- as.vector(sapply(filter_values, max))
  interval_width <- (filter_max - filter_min)/num_intervals
  vertex_index <- 0
  level_of_vertex <- c()
  points_in_vertex <- list()
  points_in_level_set <- vector("list", num_levelsets)
  vertices_in_level_set <- vector("list", num_levelsets)
  method <- meth
  object1<-object
  d<-num_bins_when_clustering
  for (lsfi in 1:num_levelsets) {
    lsmi <- lsmi_from_lsfi(lsfi, num_intervals)
    lsfmin <- filter_min + (lsmi - 1) * interval_width - 
      0.5 * interval_width * percent_overlap/100
    lsfmax <- lsfmin + interval_width + interval_width * 
      percent_overlap/100
    for (point_index in 1:num_points) {
      if (all(lsfmin <= filter_values[point_index, ] & 
              filter_values[point_index, ] <= lsfmax)) {
        points_in_level_set[[lsfi]] <- c(points_in_level_set[[lsfi]], 
                                         point_index)
      }
    }
    points_in_this_level <- points_in_level_set[[lsfi]]
    num_points_in_this_level <- length(points_in_level_set[[lsfi]])
    if (num_points_in_this_level == 0) {
      num_vertices_in_this_level <- 0
    }
    if (num_points_in_this_level == 1) {
      num_vertices_in_this_level <- 1
      level_internal_indices <- c(1)
      level_external_indices <- points_in_level_set[[lsfi]]
    }
    if (num_points_in_this_level > 1) {
      
      if (method=="Forgy"){
        #level_dist_object <- 1-abs(dist_object[points_in_this_level,points_in_this_level])
        level_dist_object <- sqrt(1-(dist_object[points_in_this_level,points_in_this_level])^2)
        #diag(level_dist_object)<-0
        level_max_dist <- max(level_dist_object)
        
        if (dim(level_dist_object)[1] < d){
          set.seed(12345)
          level_kclust <- kmeans(level_dist_object,centers=1, algorithm="Forgy",iter.max=50)
        }
        if (dim(level_dist_object)[1] >= d){
          set.seed(12345)
          level_kclust <- kmeans(level_dist_object,centers=num_bins_when_clustering, algorithm="Forgy",iter.max=50,nstart=10)
        }
        level_internal_indices <- level_kclust$cluster
        level_external_indices <- points_in_this_level
        num_vertices_in_this_level <- max(level_internal_indices)
      }
      
      if (method=="Lloyd"){
        #level_dist_object <- 1-abs(dist_object[points_in_this_level,points_in_this_level])
        #diag(level_dist_object)<-0
        level_dist_object <- sqrt(1-(dist_object[points_in_this_level,points_in_this_level])^2)
        level_max_dist <- max(level_dist_object)
        
        if (dim(level_dist_object)[1] < d){
          set.seed(12345)
          level_kclust <- kmeans(level_dist_object,centers=1, algorithm="Lloyd",iter.max=50)
        }
        if (dim(level_dist_object)[1] >= d){
          set.seed(12345)
          level_kclust <- kmeans(level_dist_object,centers=num_bins_when_clustering, algorithm="Lloyd",iter.max=50,nstart=10)
        }
        level_internal_indices <- level_kclust$cluster
        level_external_indices <- points_in_this_level
        num_vertices_in_this_level <- max(level_internal_indices)
      }
      
      if (method=="Hartigan-Wong"){
        #level_dist_object <- 1-abs(dist_object[points_in_this_level,points_in_this_level])
        #diag(level_dist_object)<-0
        level_dist_object <- sqrt(1-(dist_object[points_in_this_level,points_in_this_level])^2)
        level_max_dist <- max(level_dist_object)
        
        if (dim(level_dist_object)[1] < d){
          set.seed(12345)
          level_kclust <- kmeans(level_dist_object,centers=1, algorithm="Hartigan-Wong",iter.max=50)
        }
        if (dim(level_dist_object)[1] >= d){
          set.seed(12345)
          level_kclust <- kmeans(level_dist_object,centers=num_bins_when_clustering, algorithm="Hartigan-Wong",iter.max=50,nstart=10)
        }
        level_internal_indices <- level_kclust$cluster
        level_external_indices <- points_in_this_level
        num_vertices_in_this_level <- max(level_internal_indices)
      }
      
      if (method=="MacQueen"){
        #level_dist_object <- 1-abs(dist_object[points_in_this_level,points_in_this_level])
        #diag(level_dist_object)<-0
        level_dist_object <- sqrt(1-(dist_object[points_in_this_level,points_in_this_level])^2)
        level_max_dist <- max(level_dist_object)
        
        if (dim(level_dist_object)[1] < d){
          set.seed(12345)
          level_kclust <- kmeans(level_dist_object,centers=1, algorithm="MacQueen",iter.max=50)
        }
        if (dim(level_dist_object)[1] >= d){
          set.seed(12345)
          level_kclust <- kmeans(level_dist_object,centers=num_bins_when_clustering, algorithm="MacQueen",iter.max=50,nstart=10)
        }
        level_internal_indices <- level_kclust$cluster
        level_external_indices <- points_in_this_level
        num_vertices_in_this_level <- max(level_internal_indices)
      }
      
      if (method=="hclust"){
        #level_dist_object <- as.dist(as.matrix(dist_object)[points_in_this_level, 
        #    points_in_this_level])
        #level_dist_object <- as.dist(1-abs(dist_object[points_in_this_level, points_in_this_level]))
        level_dist_object <-as.dist(sqrt(1-(dist_object[points_in_this_level,points_in_this_level])^2))
        level_max_dist <- max(level_dist_object)
        set.seed(12345)
        level_hclust <- hclust(level_dist_object, method = "single")
        level_heights <- level_hclust$height
        level_cutoff <- cluster_cutoff_at_first_empty_bin(level_heights, 
                                                          level_max_dist, num_bins_when_clustering)
        level_external_indices <- points_in_this_level[level_hclust$order]
        level_internal_indices <- as.vector(cutree(list(merge = level_hclust$merge, 
                                                        height = level_hclust$height, labels = level_external_indices), 
                                                   h = level_cutoff))
        num_vertices_in_this_level <- max(level_internal_indices)
      }
    }
    if (num_vertices_in_this_level > 0) {
      vertices_in_level_set[[lsfi]] <- vertex_index + (1:num_vertices_in_this_level)
      for (j in 1:num_vertices_in_this_level) {
        vertex_index <- vertex_index + 1
        level_of_vertex[vertex_index] <- lsfi
        points_in_vertex[[vertex_index]] <- level_external_indices[level_internal_indices == j]
      }
    }
  }
  adja <- mat.or.vec(vertex_index, vertex_index)
  for (lsfi in 1:num_levelsets) {
    lsmi <- lsmi_from_lsfi(lsfi, num_intervals)
    for (k in 1:filter_output_dim) {
      if (lsmi[k] < num_intervals[k]) {
        lsmi_adjacent <- lsmi + diag(filter_output_dim)[, 
                                                        k]
        lsfi_adjacent <- lsfi_from_lsmi(lsmi_adjacent, 
                                        num_intervals)
      }
      else {
        next
      }
      if (length(vertices_in_level_set[[lsfi]]) < 1 | length(vertices_in_level_set[[lsfi_adjacent]]) < 
          1) {
        next
      }
      for (v1 in vertices_in_level_set[[lsfi]]) {
        for (v2 in vertices_in_level_set[[lsfi_adjacent]]) {
          adja[v1, v2] <- (length(intersect(points_in_vertex[[v1]], 
                                            points_in_vertex[[v2]])) > 0)
          adja[v2, v1] <- adja[v1, v2]
        }
      }
      for (v1 in vertices_in_level_set[[lsfi]]) {
        for (v2 in vertices_in_level_set[[lsfi_adjacent]]) {
          ind0<-intersect(points_in_vertex[[v1]], points_in_vertex[[v2]])
          ind1<-unique(c(points_in_vertex[[v1]], points_in_vertex[[v2]]))
          indmax<-max(length(points_in_vertex[[v1]]), length(points_in_vertex[[v2]]))
          
          if (length(ind1)>=2){
            adja00<-abs(cor(object1[,ind1]))
            diag(adja00)<-0
            delta<-mean(apply(adja00,1,mean))
          }else {delta<-0}
          
          if (length(ind0)>=pero0*indmax | delta>=corm0) {
            adja[v1, v2] <- 1
          }
          
          if (adja[v1,v2]==1){adja[v1,v2]=delta}
          
          adja[v2, v1] <- adja[v1, v2]
        }	
      }
    }
  }
  mapperoutput <- list(adjacency = adja, num_vertices = vertex_index, 
                       level_of_vertex = level_of_vertex, points_in_vertex = points_in_vertex, 
                       points_in_level_set = points_in_level_set, vertices_in_level_set = vertices_in_level_set)
  class(mapperoutput) <- "TDAmapper"
  return(mapperoutput)
}

defmacro<-function (..., expr) 
{
  expr <- substitute(expr)
  a <- substitute(list(...))[-1]
  nn <- names(a)
  if (is.null(nn)) 
    nn <- rep("", length(a))
  for (i in 1:length(a)) {
    if (nn[i] == "") {
      nn[i] <- paste(a[[i]])
      msg <- paste(a[[i]], "not supplied")
      a[[i]] <- substitute(stop(foo), list(foo = msg))
    }
    if (nn[i] == "DOTS") {
      nn[i] <- "..."
      a[[i]] <- formals(function(...) {
      })[[1]]
    }
  }
  names(a) <- nn
  a <- as.list(a)
  ff <- eval(substitute(function() {
    tmp <- substitute(body)
    eval(tmp, parent.frame())
  }, list(body = expr)))
  formals(ff) <- a
  mm <- match.call()
  mm$expr <- NULL
  mm[[1]] <- as.name("macro")
  attr(ff, "source") <- c(deparse(mm), deparse(expr))
  ff
}
