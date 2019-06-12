# This script 
library(data.table)
source("../Bayes_R/BINE_extrct_vis/membership_functions.R")

tot_samples=8000
PMAX=200
ACTIVE_NEURONS_THRESHOLD=15
FEATURE_NAMES=c("LATERAL_INDEX",
                "EVENT_COHERENCE",
				"ASSEMBLY_ACTIVITY",
				"ASSEMBLY_SIZE",
				"ASSEMBLY_COMPACTNESS",
				"ASSEMBLY_NOISE",
				"MEAN_ASSEMBLY_DURATION",
				"MAX_ASSEMBLY_DURATION",
				"NUMBER_OF_ASSEMBLY_EVENTS",
				"FIRING_IN_THE_BEGINNING");				

prefix_list = list.files(pattern = pattern)

for (i in 1:length(prefix_list)){
  # folder containing the results of BPTS
  prefix = prefix_list [i]
  folder=paste(prefix,sep="")
  
  # load the prefix list and determine the index of the prefix under consideration in the list
  list.files(pattern = pattern)
  prefix_ID=which(prefix==prefix_list)
  
  # number of samples extracted from the MCMC
  nsamples=20
  
  # load input binary matrix S
  s<-read.table(paste(data_location, prefix, spikes_suffix ,sep=""));
  
  # load results
  mem_traj<-fread(paste(folder,"/membership_traj.dat",sep=""),skip=tot_samples-nsamples)
  centers<-fread(paste(data_location ,prefix,centers_suffix,sep=""));
  
  if (reg) centers[,1:2]<-read.table(paste(data_location, prefix,"_registered_cell_centers.dat",sep=""));
  
  selection<-read.table(paste(folder,"/selection.dat",sep=""))$V1
  centers<-centers[selection+1,]

  # write centers and membership
  write.table(centers,file=paste(folder,"/centers.dat",sep=""),row.names=FALSE,col.names=FALSE)

  # Build the omega tensor (number of ensembles) X (time frames) X (samples)
  # note: PMAX is set to 200. 

  omega_traj=fread(file=paste(folder,"omega_traj.dat",sep=""),skip=tot_samples-nsamples)
  times=ncol(omegatraj)

  omega_traj=lapply(split(omega_traj,
                          rep(1:nsamples,each=PMAX)
					),
					matrix,PMAX,times)

  omega_mem_traj = vector("list",nsamples)
  for(i in 1:nsamples) omega_mem_traj[[i]] = list(omega=omega_traj[[i]],mem=mem_traj[i,])
		  

  omega_extended_tensor = array(0,dim=c(PMAX,ncol(s),nsamples))
  omega_extended_tensor[,colSums(s[selection+1,])>ACTIVE_NEURONS_THRESHOLD,] = omega_traj

  feature_list_bysamples = vector("list",length(FEATURE_NAMES))
  names(feature_list_bysamples)=FEATURE_NAMES


}


