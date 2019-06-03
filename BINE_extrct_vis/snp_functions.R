mean_corr <- function(folder, save_name) {
  files <- list.files(path = folder, pattern = "all_cells_spikes.dat")
  mean_corrs <- vector("numeric", length(files))
  for (i in 1:length(files)) {
    mat <- cor(t(as.matrix(read.table(paste(folder, files[[i]], sep ="")))))
    mean_corrs [i] <- mean(abs(mat[lower.tri(mat)]), na.rm = TRUE)
  }
  write.table(mean_corrs, file = paste("mean_across_fish/mean_corr_", save_name, ".dat", sep =""))
}

plot_all_F <- function(pattern) {
  files <- list.files(pattern = pattern)
  par(mfrow = c(3,2))
  for (i in 1:length(files)) {
    plot(apply(X = read.table(file = paste(files [[i]], "/F.dat" ,sep ="")), 1, FUN = mean), type ="l", ylab = "likelyhood", xlab = paste(files [[i]]))
  }
}

Check_midlines <- function(midline, folder) {
  files <- list.files(folder, pattern = "cut_ordered")
  par(mfrow = c(2,4))
  for (i in 1:length(files)) {
    centers <- as.matrix(read.table(paste(folder, files [i], sep ="")))
    plot(centers [,1], centers [,2], cex = 0.2, pch = 19)
    abline(a = midline [i,1], b= midline [i,2], col = "red", lwd = 2 )
  }
}



###########################################################################################################################
# functions to caluculate response amplitude
seg1d <- function(trace,threshold=0){
  counter=1;
  out_trace=rep(0,length(trace));
  i=1;
  start=trace[1];
  
  while(i<length(trace)){
    if((i==1 & start==0) | i>1) while(trace[i]<threshold & i<length(trace)) i=i+1;
    while(trace[i]>threshold & i<=length(trace)){
      out_trace[i]=counter;
      i=i+1;
    }
    counter=counter+1;
  }
  return(out_trace)
}



response_amplitude <- function(trace) {
  seg_trace <- seg1d(trace , threshold = 0.1)
  amp <- by(trace, seg_trace, max)
  return(as.numeric(amp))
}
###########################################################################################################################
###########################################################################################################################
#function to calculate all single neurons properties



snp <- function(folder, save_name) {
  files <- list.files(path = folder, pattern = "all_cells_spikes.dat")
  mean_corrs <- vector("numeric", length(files))
  freq <- vector("numeric", length(files))
  cell_num <- vector("numeric", length(files))
  for (i in 1:length(files)) {
    s <- as.matrix(read.table(paste(folder, files[[i]], sep ="")))
    mat <- cor(t(s))
    mean_corrs [i] <- mean(abs(mat[lower.tri(mat)]), na.rm = TRUE)
    cell_num [i] <- nrow(s)
    freq [i]  <- mean(rowSums(s)/60) 
  }
  write.table(mean_corrs, file = paste("mean_across_fish/mean_corr_", save_name, ".dat", sep =""))
  write.table(cell_num,   file = paste("mean_across_fish/cell_num_", save_name, ".dat", sep =""))
  write.table(freq,  file = paste("mean_across_fish/sn_freq_", save_name, ".dat", sep =""))
}
