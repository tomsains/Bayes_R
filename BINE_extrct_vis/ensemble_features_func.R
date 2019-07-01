# functions for extracting ensemble features from samples

#librarys..........................................................................................................................
library(sm)
library(zoo)


# a function to calulate the mode of of the samples (Ie. the most likely sample value)
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

guassian_smooth <- function(x){
  if (sum(x) == 0) {
     sm_tr <- rep(0, length(x))
  }
  else{
    sm_tr = density(which(x == 1), width = 10,  from = 0, to = length(x), n = length(x))$y
  }
  return(sm_tr)
}

get_omega_smooth <- function(o) {
  smooth_omega <- t(apply(o, 1, guassian_smooth))
  binary_smooth <- smooth_omega*100 > 0.13
  print(1)
  return(binary_smooth)
}

get_omega_smooth_tensor <- function(o) {
  #tensor <- apply(o, 3, get_omega_smooth)
  pb <- txtProgressBar(min = 0, max = dim(omega_extended_tensor)[3], style = 3)
  tensor <- o
  for (i in 1:dim(omega_extended_tensor)[3]) {
   tensor [,,i] <- get_omega_smooth(o [,,i])
   setTxtProgressBar(pb, i)
  }
  return(tensor)
}

get_s_Smooth <- function(s) {
  s <- s [selection+1,]
  smooth_s <- t(apply(as.matrix(s), 1, guassian_smooth))
  cors <- cor(t(smooth_s*100))
  return(abs(mean(cors)))
}


####################################################################################################################################
# Paramters/metrics by samples
####################################################################################################################################

LI <- function(x, mid = midline [1,]) {
  # X the omega_traj_object
  LI <- vector("numeric", length = cluster_num_per_sample)
  for (i in 1:cluster_num_per_sample){
    assembly_mem <- centers [as.vector(as.numeric(x$mem) == i),]
    if (is.null(dim(assembly_mem))) {
      LI [i] <- NA
    }
    else {
      LH = sum(assembly_mem [,2] > (mid$slope *assembly_mem [,1] +mid$intercept))
      RH <- sum(assembly_mem [,2] < (mid$slope*assembly_mem [,1] +mid$intercept))
      ALL <- nrow(assembly_mem)
      LI [i] <- ((LH-RH)/ALL)
    }
  }
  return(LI)
}

NULL_LI <- function(x, mid = midline [1,]) {
  LI <- vector("numeric", length = cluster_num_per_sample)
  for (i in 1:as.numeric(cluster_num_per_sample)){
    assembly_mem <- centers [sample(as.vector(as.numeric(x$mem) == i)),]
    if (is.null(dim(assembly_mem))) {
      LI [i] <- NA
    }
    else {
      LH = sum(assembly_mem [,2] > (mid$slope *assembly_mem [,1] +mid$intercept))
      RH <- sum(assembly_mem [,2] < (mid$slope*assembly_mem [,1] +mid$intercept))
      ALL <- nrow(assembly_mem)
      LI [i] <- ((LH-RH)/ALL)
    }
  }
  return(LI)
}



beg_fire_ten <- function(omega, bin_size=50){
  beg_rate <- rowSums(omega [,1:bin_size])
  av_rate <- vector(length = nrow(omega))
  for (i in 1:nrow(omega)) {
    trace = omega [i,] [10061:ncol(omega)]
    reshaped_end = matrix(trace, ncol = bin_size)
    av_rate [i] <- max(rowSums(reshaped_end))
  }
  return((beg_rate-av_rate)/50)
}

response_duration <- function(t) {
  runs <- rle(t)
  return(runs$lengths[runs$values==1])
}

number_bursting_Events <- function(omega) {
  durations <- function(x) {apply(x, 1, response_duration)}
  durations <- durations(omega)
  number_events <-lapply(durations, function(x) {x <- x [x > 2]; length(x)})
  return(unlist(number_events))
}


Assembly_Activity <- function(omega) {
  assembly_freq <- apply(omega, 1, function(x) {sum(x)/ncol(omega)})
  return(assembly_freq)
}


Assembly_size <- function(x) {
  size <- vector("numeric", length = cluster_num_per_sample)
  for (i in 1:cluster_num_per_sample){
    size [i] <- sum(as.vector(as.numeric(x$mem) == i))
  }
  return(size)
}

Assembley_Compactness <- function(x) {
  Compactness <- vector("numeric", length = cluster_num_per_sample)
  for (i in 1:cluster_num_per_sample){
    assembly_mem <- centers [as.vector(as.numeric(x$mem) == i),]
    if (is.null(dim(assembly_mem))) {
      Compactness [i] <- NA
    }
    else {
      Compactness [i] <- (eigen(cov(assembly_mem))$values [1]*eigen(cov(assembly_mem))$values [2])
    }
  }
    return(Compactness)
}


new_coherence <- function(x, omega_ext_smooth, thresh = 0.01) {
  coherence <- vector("numeric", length = cluster_num_per_sample)
  for (i in 1:cluster_num_per_sample){
    assembly_mem <- centers [as.vector(as.numeric(x$mem) == i),]
    w_start <- which(diff(omega_ext_smooth > thresh)  == 1)
    w_end <- which(diff(omega_extended_smooth > thresh) == -1)
    awa <- matrix(NA,nrow = nrow(cells), ncol = length(w_start))
    if (length(w_start) == 0) {
      return(0)
    }
  else{
    for (i in 1:length(w_start)) {awa [,i] <- rowSums(cells [,w_start[i]:w_end[i]])}
    return(mean(colSums(awa >0)/nrow(awa)))
  }
}


#####################################################################################################################################
# Plots
#####################################################################################################################################

LI_vs_NULL <- function(){
  LI <- as.vector(feature_list_bysamples$LATERAL_INDEX)
  LI <- LI [!is.na(LI)]
  LI_NULL <- as.vector(feature_list_bysamples$NULL_LATERAL_INDEX)
  LI_NULL <- LI_NULL [!is.na(LI_NULL)]
  dims = dim(feature_list_bysamples$LATERAL_INDEX)
  sm.density.compare(abs(c(LI, LI_NULL)),  rep(1:2,each =  length(LI), xlim = c(0,1)), xlab="Lateral Index")
}