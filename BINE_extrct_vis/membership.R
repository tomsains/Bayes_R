# This script 
source("../Bayes_R/BINE_extrct_vis/membership_functions.R")


prefix_list = list.files(pattern = pattern)

for (i in 1:length(prefix_list)){
  # folder containing the results of BPTS

  prefix = prefix_list [i]
  folder=paste(prefix,sep="")
  
  
  # load the prefix list and determine the index of the prefix under consideration in the list
  list.files(pattern = pattern)
  prefix_ID=which(prefix==prefix_list)
  
  # excluded ensembles because only active at the beginning of the recording
  exclusions<-vector("list",length(prefix_list))
  print(i)
  
 for (j in 1:length(exclusions)) {
   exclusions [[j]] = c(-1)
 }
  print(i)
  exc=exclusions[[which(prefix==prefix_list)]]
  
  # file containing infos about side and anterior-posterior location of the ensemble
  #side_AP<-read.table(paste(folder,"/side_AP.dat",sep=""))
  #names(side_AP)<-c("enID","RL","AP")
  
  #Thresholds used for the analysis
  pmu_thresh=0.005
  lambda1_thresh=0.05
  lambda0_thresh=0.05 
  gs_thresh=5
  
  # number of samples extracted from the MCMC
  nsamples=20
  
  # load input binary matrix S
  loadS=T
  if(loadS) s<-read.table(paste(data_location, prefix, spikes_suffix ,sep=""));
  # c<- read.table(paste(data_location, prefix, cal_suffix ,sep=""));
  # load results
  mem_traj<-read.table(paste(folder,"/membership_traj.dat",sep=""))
  centers<-read.table(paste(data_location ,prefix,centers_suffix,sep=""));
  
  if (reg == TRUE) centers[,1:2]<-read.table(paste(data_location, prefix,"_registered_cell_centers.dat",sep=""));
  
  print(i)
  selection<-read.table(paste(folder,"/selection.dat",sep=""))$V1
  centers<-centers[selection+1,]
  FreeEnergy<-read.table(paste(folder,"/F.dat",sep=""))
  tot_samples=nrow(mem_traj)
  mem<-apply(mem_traj[(tot_samples-nsamples):tot_samples,],2,function(x) which.max(tabulate(x+1)))
  tab=sort(table(mem),decreasing=T)
  enIDs=as.numeric(names(tab))
  #mem=as.numeric(mem_traj[3000,]+1)
  probs<-apply(mem_traj[(tot_samples-nsamples):tot_samples,],2,function(x) max(tabulate(x+1)/length(x)))
  
  # write centers and membership
  write.table(centers,file=paste(folder,"/centers.dat",sep=""),row.names=FALSE,col.names=FALSE)
  write.table(mem,file=paste(folder,"/mem.dat",sep=""),row.names=FALSE,col.names=FALSE)
  maxcl=max(mem);
  
  # define the list of parameters nad save it
  A=plot_par();
  save(A, file=paste(folder,"/pars.RData",sep=""))
  # define the selection to be applied to tab
  sel=as.vector(A$gs.mean>gs_thresh & 
                  A$lambda1.mean>lambda1_thresh & 
                  A$pmu.mean>pmu_thresh &
                  A$lambda0.mean<lambda0_thresh &
                  !as.numeric(gsub("V","",names(A$gs.mean)))%in% exc)
  sel [is.na(sel)] = FALSE  
  # ensemble selection
  
  ensel=enIDs[sel]
  print(i)
  # calculate side_AP_ordered
  #side_AP_ord=data.frame(enID=ensel,RL=NA,AP=NA)
  #counter=0;
  #for(i in ensel){
  #	counter=counter+1
  #	j=which(side_AP$enID==i)
  #	if(length(j)==0) cat("problems in finding ensemble ID!!!!\n")
  #	side_AP_ord[counter,]=side_AP[j,]
  #}

 
  omega=get_omega();
  #cordata=get_cormat_bayes(nav=200,th=0.05)
  #graph <- cordata$graph 
  #wc <- edge.betweenness.community(graph, weight=get.data.frame(graph)$weight, directed = FALSE, bridges=TRUE)
  omega_extended = matrix(0,nrow=nrow(omega),ncol=ncol(s))
  for(k in 1:nrow(omega)){ omega_extended[k,colSums(s[selection+1,])>15]=omega[k,]}
  
  #mem2=mem
  #mem2[ !mem2 %in% ensel] = -1
  #mem3=mem2; for(i in ensel) mem3[mem2==i]=(which(ensel[order(subnet_mem)]==i))
  #write.table(mem3,file=paste(folder,"/mem3.dat",sep=""),col.names=F,row.names=F)
  
#  list_of_images <- vector("list", length = length(ensel))
 # for (j in 1:length(ensel)) {
  #  col = brewer.pal(n = 8, name = "Dark2") [4]
   # list_of_images [[j]] <- ggplot(data.frame(centers), aes(x = V1,y= V2)) +geom_point(size = 8, alpha = 0.4, col ="#666666")  + geom_point(data = data.frame(centers [mem==ensel [j] & probs>0.99,]), aes(x = V1,y=V2), size = 10, col =col)+theme_classic()
  #}
  
#  for (a in 1:length(list_of_images)) {
 #   tiff(paste("../../../../figure_materials/Chapter_1/Figure_3/Assembly_map_Matches_matched_to_firing",folder,"_", a,".tif", sep =""), width = 1000,height= 1000)
  #  print(list_of_images [[a]])
   # dev.off()
  #}
  
  
 # for (f in 1:length(ensel)){
  #  tiff(paste("../../../../figure_materials/Chapter_1/Figure_4/Assembly_firing", folder,"_", f,".tif", sep =""), width = 1200,height= 300)
   # plot_assembly_w_omega(f)
    #dev.off()
  #}
 
  print(i)
  png(paste(folder, "/maps.png", sep =""), width = 2000, height = 1000)
  plot_assemblies(-1, .99)
  dev.off()
  
  
  png(paste(folder, "/omega.png", sep =""), width = 1500, height = 1000)
  plot_omega()
  dev.off()
  
  png(paste(folder, "/sorted_raster.png", sep =""), width = 3000, height = 2000)
  sort_by_correlation(omega)
  dev.off()
  
  png(paste(folder, "/raster.png", sep =""), width = 1500, height = 1000)
  raster_plot(s, paste(folder))
  dev.off()
  
  #new_co <- vector("numeric", length = length(ensel))
  #for (i in 1:length(ensel)) {new_co [i] <- unlist(new_coherence(i, 0.2))}
  print(i)
  
  beg_Fire = beg_fire_ten(omega_extended, bin_size = 50) 
  beg_firing = data.frame(Beg_index = beg_Fire, Beg_thresh = beg_Fire > 0.5)
  
  write.table(Assembly_centers(centers), file = paste(folder, '/assembly_centers.dat', sep = ""))
  write.table(beg_firing, file = paste(folder, '/begin_firing.dat', sep = ""))
  write.table(omega [!beg_firing$Beg_thresh,], paste(folder, '/omega_be.dat', sep = ""))
  write.table(omega_extended [!beg_firing$Beg_thresh,], paste(folder, '/omega_ext_be.dat', sep = ""))
  write.table(ensel [!beg_firing$Beg_thresh], paste(folder, '/ensel_be.dat', sep = ""))
  write.table(ensel [beg_firing$Beg_thresh ], paste(folder, '/ensel_excluded.dat', sep = ""))
  write.table(omega [beg_firing$Beg_thresh,], paste(folder, '/omega_excluded.dat', sep = ""))
  
  write.table(omega, file = paste(folder, "/omega.dat", sep =""))
  write.table(omega_extended, file = paste(folder, "/omega_extended.dat", sep =""))
  write.table(eigen_decomp_of_cov(), file = paste(folder, "/eigen_decomp_of_cov.dat", sep ="") )
  write.table(total_area(), file = paste(folder, "/ensem_area.dat", sep =""))
  write.table(within_cluster_corr(), file = paste(folder, "/wcc.dat", sep ="")) 
  write.table(apply(omega_extended, 1, function(x) {sum(x)/ncol(omega_extended)}), file = paste(folder, "/Assembly_Freq.dat", sep ="")) 
  write.table(x = LI(mid = midline [i,]), file = paste(folder, "/LI.dat", sep =""))
  
  print(i)
  png(paste("../../../../figure_materials/Chapter_1/Figure_4/", folder, "Assembly_plot.png", sep =""), width = 1500, height = 500)
  plot_assembly_w_omega(1)
  dev.off()
  
  
  #png(paste(folder, "/exlusion_plot.png", sep =""), width = 1500, height = 1000)
  #check_exclusions(paste(folder))
  #dev.off()
  ensel = ensel [!beg_firing$Beg_thresh]
  omega = omega [!beg_firing$Beg_thresh,]
  omega_extended = omega_extended [!beg_firing$Beg_thresh,]

  png(paste(folder, "/sorted_raster.png", sep =""), width = 3000, height = 2000)
  sort_by_correlation(omega)
  dev.off()
  
  png(paste(folder, "/free_neurons.png", sep =""), width = 3000, height = 2000)
  raster_plot(s [probs < 0.99,], paste(folder))
  dev.off()
  
  #write.table(cbind(reverberating_activity(omega_extended = omega_extended), new_co), paste(folder, '/assembly_duration_and_co.dat', sep = ""))
  
 
  
  
  #uncomment these to make changes to C1 F4-5
  # lateral_index_calc
 
  #Li <- LI(mid = midline [i,])
  #NULL_Li <- NULL_LI(mid = midline [i,], samples = 2000)
  #LI_data_frame = data.frame(c(rep("Li", length(Li)), rep("NULL_LI", length(as.vector(NULL_Li)))))
  #colnames(LI_data_frame) <- c("sample")
  #LI_data_frame$LI_value <- c(Li, NULL_Li)
  #write.table(LI_data_frame, paste("../../../../figure_materials/Chapter_1/Figure_4/", folder, "_LI_vs_NULL.dat", sep =""))
  
  
  # Within culster correlation
  #wcc <- within_cluster_corr()
  #NULL_wcc <- NULL_within_cluster_corr(samples = 2000)
  #wcc_data_frame  = data.frame(c(rep("wcc", length(wcc)), rep("NULL_wcc", length(as.vector(NULL_wcc)))))
  #colnames(wcc_data_frame) <- c("sample")
  #wcc_data_frame$wcc_value<- c(wcc, NULL_wcc)
  #write.table(wcc_data_frame, paste("../../../../figure_materials/Chapter_1/Figure_4/", folder, "_wcc_vs_NULL.dat", sep = ""))
  
  # compactness
  #compct <- eigen_decomp_of_cov()
  #NULL_compct <- NULL_eigen_decomp_of_cov(samples = 2000)
  
  #compct_data_frame  = data.frame(c(rep("comp", length(compct)), rep("NULL_comp", length(as.vector(NULL_compct)))))
  #colnames(compct_data_frame) <- c("sample")
  #compct_data_frame$comp_value<- c(compct, NULL_compct)
  #write.table(compct_data_frame, paste("../../../../figure_materials/Chapter_1/Figure_4/", folder, "_compct_vs_NULL.dat", sep =""))
  
}


