# This script 
library(data.table)
library(bigmemory)
source("../Bayes_R/BINE_extrct_vis/ensemble_features_func.R")


PMAX=200
ACTIVE_NEURONS_THRESHOLD=15
FEATURE_NAMES=c("PREFIX_ID",
                "GENOTYPE",
                "REARING_CONDITIONS",
                "LATERAL_INDEX",
                "NULL_LATERAL_INDEX",
                "EVENT_COHERENCE",
        				"ASSEMBLY_ACTIVITY",
        				"ASSEMBLY_SIZE",
        				"ASSEMBLY_COMPACTNESS",
        				"ASSEMBLY_NOISE",
        				"MEAN_ASSEMBLY_DURATION",
        				"NUMBER_OF_ASSEMBLY_EVENTS",
        				"FIRING_IN_THE_BEGINNING",
        				"WITHIN_ASSEMBLY_CORR");				

prefix_list = list.files(pattern = pattern)

for (i in 4:length(prefix_list)){
  # folder containing the results of BPTS
  prefix = prefix_list [2]
  folder=paste(prefix,sep="")
  
  # load the prefix list and determine the index of the prefix under consideration in the list
  list.files(pattern = pattern)
  prefix_ID=which(prefix==prefix_list)
  
  # number of samples extracted from the MCMC
  nsamples=200
  
  # load input binary matrix S
  s<-read.big.matrix(filename = paste(data_location, prefix, spikes_suffix ,sep=""), sep = " ", header =FALSE, type = "short");
  
  # load results
  tot_samples_mem=as.numeric(system(paste("cat ",folder,"/membership_traj.dat | wc -l",sep=""), intern = TRUE))
  mem_traj<-fread(paste(folder,"/membership_traj.dat",sep=""),skip=tot_samples_mem-nsamples)+1
  centers<-as.matrix(fread(paste(data_location ,prefix,centers_suffix,sep="")), ncol =3);
  
  if (reg) centers[,1:2]<-read.table(paste(data_location, prefix,"_registered_cell_centers.dat",sep=""));
  
  selection<-read.table(paste(folder,"/selection.dat",sep=""))$V1
  centers<-centers[selection+1,]

  # write centers and membership
  write.table(centers,file=paste(folder,"/centers.dat",sep=""),row.names=FALSE,col.names=FALSE)

  # Build the omega tensor (number of ensembles) X (time frames) X (samples)
  # note: PMAX is set to 200. 
  
  tot_samples_omega= as.numeric(system(paste("cat ",folder,"/omega_traj.dat | wc -l",sep=""), intern = TRUE))/PMAX
  omega_traj=as.matrix(fread(file=paste(folder,"/omega_traj.dat",sep=""),skip=(tot_samples_omega-nsamples)*PMAX))
  #omega_traj = r(filename =paste(folder,"/omega_traj.dat",sep=""), sep = " ", header = FALSE, skip=tot_samples-nsamples)
  times=ncol(omega_traj)

  omega_traj=lapply(split(omega_traj,
                          rep(1:nsamples,each=PMAX)
					),
					matrix,PMAX,times)
  
  
  
  # not all samples will have the same cluster number - here we try to select those samples that have the most likely cluster 
  cluster_num_per_sample  <- Mode(unlist(lapply(omega_traj, function(x) sum(!is.na(rowSums(x))))))
  samples_with_modecl <- lapply(omega_traj, function(x) sum(!is.na(rowSums(x)))) == cluster_num_per_sample
  omega_traj <- omega_traj [samples_with_modecl]
  mem_traj <- mem_traj [samples_with_modecl,]
  nsamples <- length(omega_traj)
  
  omega_traj <- lapply(omega_traj, function(x) {x [!is.na(rowSums(x)),]})
  omega_mem_traj = vector("list",nsamples)
  for(i in 1:nsamples) omega_mem_traj[[i]] = list(omega=omega_traj[[i]],mem=mem_traj[i,])
		  
  # create omega extend for each sample
  omega_extended_tensor = array(0,dim=c(nrow(omega_traj$`1`),ncol(s),nsamples))
  for (i in 1:nsamples)  omega_extended_tensor[,colSums(s[selection+1,])>ACTIVE_NEURONS_THRESHOLD,i] = omega_traj [[i]]
  
  smooth_omega_ext <- get_omega_smooth_tensor(omega_extended_tensor)
  smooth_omega_mem_traj <- omega_mem_traj
  for (i in 1:dim(smooth_omega_ext) [3]) {smooth_omega_mem_traj [[i]]$omega <- smooth_omega_ext [,,i]}

  feature_list_bysamples = vector("list",length(FEATURE_NAMES))
  names(feature_list_bysamples)=FEATURE_NAMES
  
  feature_list_bysamples$PREFIX_ID = folder
  feature_list_bysamples$GENOTYPE = genotype
  feature_list_bysamples$REARING_CONDITIONS = Rearing_conditions
  feature_list_bysamples$LATERAL_INDEX = as.matrix(do.call(cbind, lapply(omega_mem_traj ,FUN = LI, mid = midline [prefix_ID,])))
  feature_list_bysamples$NULL_LATERAL_INDEX = do.call(cbind, lapply(omega_mem_traj, NULL_LI, mid = midline [prefix_ID,]))
  feature_list_bysamples$FIRING_IN_THE_BEGINNING <- apply(smooth_omega_ext, 3, beg_fire_ten)
  feature_list_bysamples$NUMBER_OF_ASSEMBLY_EVENTS  <- apply(smooth_omega_ext, 3, number_bursting_Events)
  feature_list_bysamples$ASSEMBLY_ACTIVITY <- apply(omega_extended_tensor, 3, Assembly_Activity)
  feature_list_bysamples$ASSEMBLY_SIZE <- as.matrix(do.call(cbind, lapply(omega_mem_traj ,FUN = Assembly_size)))
  feature_list_bysamples$ASSEMBLY_COMPACTNESS <- as.matrix(do.call(cbind, lapply(omega_mem_traj ,FUN = Assembley_Compactness)))
  feature_list_bysamples$EVENT_COHERENCE <- as.matrix(do.call(cbind, lapply(smooth_omega_mem_traj ,FUN = new_coherence)))
  feature_list_bysamples$MEAN_ASSEMBLY_DURATION <- apply(smooth_omega_ext, 3, function(x) {unlist(lapply(apply(x, 1, response_duration), mean))})
  feature_list_bysamples$WITHIN_ASSEMBLY_CORR <- lapply(smooth_omega_mem_traj, within_cl_corr)
  save(feature_list_bysamples,file =  paste(folder, "/", folder, "_feature_list_bysamples.RData", sep = ""))
  
  rm(list = c("feature_list_bysamples", "s", "omega_extended_tensor", "smooth_omega_ext",   "smooth_omega_mem_traj", "omega_mem_traj", "omega_traj", "centers"))
  gc()
}


