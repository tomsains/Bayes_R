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
  
 for (i in 1:length(exclusions)) {
   exclusions [[i]] = c(-1)
 }
  
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
  # ensemble selection
  ensel=enIDs[sel]
  
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
  for(i in 1:nrow(omega)) omega_extended[i,colSums(s[selection+1,])>15]=omega[i,]
  
  #mem2=mem
  #mem2[ !mem2 %in% ensel] = -1
  #mem3=mem2; for(i in ensel) mem3[mem2==i]=(which(ensel[order(subnet_mem)]==i))
  #write.table(mem3,file=paste(folder,"/mem3.dat",sep=""),col.names=F,row.names=F)
  
 
  
  png(paste(folder, "/maps.png", sep =""), width = 1500, height = 1000)
  plot_assemblies(-1, .99)
  dev.off()
  
  
  png(paste(folder, "/omega.png", sep =""), width = 1500, height = 1000)
  plot_omega()
  dev.off()
  
  png(paste(folder, "/sorted_raster.png", sep =""), width = 1500, height = 1000)
  sort_by_correlation(omega)
  dev.off()
  
  png(paste(folder, "/raster.png", sep =""), width = 1500, height = 1000)
  raster_plot(s, paste(folder))
  dev.off()
  
  new_co <- vector("numeric", length = length(ensel))
  for (i in 1:length(ensel)) {new_co [i] <- unlist(new_coherence(i, 0.2))}
  
  
  beg_Fire = beg_fire(omega_extended, bin_size = 50) 
  beg_firing = data.frame(Beg_index = beg_Fire, Beg_thresh = beg_Fire > 0.3)
  
  write.table(beg_firing, file = paste(folder, '/begin_firing.dat', sep = ""))
  write.table(omega [!beg_firing$Beg_thresh,], paste(folder, '/omega_be.dat', sep = ""))
  write.table(omega_extended [!beg_firing$Beg_thresh,], paste(folder, '/omega_ext_be.dat', sep = ""))
  write.table(ensel [!beg_firing$Beg_thresh], paste(folder, '/ensel_be.dat', sep = ""))
  write.table(ensel [beg_firing$Beg_thresh ], paste(folder, '/ensel_excluded.dat', sep = ""))
  write.table(omega [beg_firing$Beg_thresh,], paste(folder, '/omega_excluded.dat', sep = ""))
  
  write.table(cbind(reverberating_activity(omega_extended = omega_extended), new_co), paste(folder, '/assembly_duration_and_co.dat', sep = ""))

  
  
  png(paste(folder, "/exlusion_plot.png", sep =""), width = 1500, height = 1000)
  check_exclusions(paste(folder))
  dev.off()
  
  
  write.table(omega, file = paste(folder, "/omega.dat", sep =""))
  write.table(omega_extended, file = paste(folder, "/omega_extended.dat", sep =""))
  write.table(eigen_decomp_of_cov(), file = paste(folder, "/eigen_decomp_of_cov.dat", sep ="") )
  write.table(total_area(), file = paste(folder, "/ensem_area.dat", sep =""))
  write.table(within_cluster_corr(), file = paste(folder, "/wcc.dat", sep ="")) 
  write.table(apply(omega_extended, 1, function(x) {sum(x)/ncol(omega_extended)}), file = paste(folder, "/Assembly_Freq.dat", sep ="")) 
  write.table(x = Lateral_index(), file = paste(folder, "/LI.dat", sep =""))
}


