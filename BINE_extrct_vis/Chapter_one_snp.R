# Chapter one Single neuron properties and population correlations - full distributions
library(data.table)
library(sm)
library(ggplot2)
library(viridis)
library(DescTools)
library(fields)
library(matrixStats)

setwd("~/../../media/thomas_sainsbury/Samsung_T5/SeG/results/Baysian_network_inference/data_for_Gibbs/7_dpf/")


circular_permutation <- function(t) {
  permute_by <- round(runif(1,1, length(t)))
  pt <- c(tail(t, -permute_by), head(t, permute_by))
  return(pt)
}

n_cells <- function(s) {
  return(nrow(s))
}

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


response_duration <- function(t) {
  runs <- rle(t)
  return(runs$lengths[runs$values==1])
}

response_amplitude <- function(trace) {
  seg_trace <- seg1d(trace , threshold = 0.1)
  amp <- by(trace, seg_trace, max)
  return(mean(as.numeric(amp)))
}

single_neuron_properties_full <- function(folder, age, genotype, rearing_condition){
  s_names <- list.files(path = folder, pattern = "all_cells_spikes.dat")
  cal_names <- list.files(path = folder, pattern = "all_cells_cal.dat")
  data_set_names <- list(length(s_names))
  data_Set_number <- list(length(s_names))
  cell_numbers <- list(length(s_names))
  firing_frequency <- list(length(s_names))
  response_amp <- list(length(s_names))
  response_dur <- list(length(s_names))
  Age <- list(length(s_names))
  Genotype <- list(length(s_names))
  Rearing_conditions <- list(length(s_names))
  
  #pb = txtProgressBar(min = 0, max = length(s_names), initial = 0) 
  
  for (i in 1:length(s_names)) {
    print(i)
    s <- as.matrix(fread(paste(folder, s_names [[i]], sep ="")))
    c <- as.matrix(fread(paste(folder, cal_names [[i]], sep ="")))
    print("loaded")
    fil_vec <- (rowSums(c > 0.2) > 1)
    c <- c  [fil_vec,]
    s <- s [fil_vec,]
    print(nrow(c) == nrow(s))
    cell_numbers [[i]] <- nrow(s)
    firing_frequency [[i]] <- (rowSums(s)/ncol(s))*4.85
    response_amp [[i]] <- apply(c , 1, response_amplitude)
    response_dur [[i]] <- unlist(lapply(apply(c > 0.2, 1, response_duration), mean))/4.85
    print("calc")
    data_set_names [[i]] <- rep(s_names [[i]], nrow(s))
    data_Set_number [[i]] <- rep(as.numeric(i), nrow(s))
    Age [[i]] <- rep(age, nrow(s))
    Genotype [[i]] <-  rep(genotype, nrow(s))
    Rearing_conditions [[i]] <- rep(rearing_condition, nrow(s))
    
  
    #setTxtProgressBar(pb,i)
  }
  data_frame_SNP_f <- as.data.frame(cbind(unlist(data_set_names),unlist(data_Set_number), unlist(Genotype), unlist(Age), unlist(Rearing_conditions), as.numeric(as.vector(unlist(firing_frequency))), as.numeric(as.vector(unlist(response_amp))), as.numeric(as.vector(unlist(response_dur)))))
  colnames(data_frame_SNP_f) <- c("ID", "F_no.", "Genotype", "Age", "Rearing_conditions", "Firing_freq", "response_amp", "response_dur")
  print(ncol(data_frame_SNP_f) == 8)
  return(list(cell_numbers, data_frame_SNP_f))
}

snp <- single_neuron_properties_full(folder = "./", age = "7_dpf", genotype = "WT", rearing_condition = "")

####################################################################################################
# Density plots 
create_density_plot <- function(con = snp [[2]]$response_amp, xlab = "enter_xlab", ylim = c(0,5), xlim = c(0,10)) {
  ggplot(data = snp[[2]], aes(x = as.numeric(as.vector(con)))) +
    geom_line(aes(group=F_no.), color= "black", stat="density", size=1, alpha=0.3) +
    geom_line( color = "#762A83", size =2, stat="density") +  theme_classic(base_size = 30) + 
    ylim(ylim) + ylab("Density") + xlab(xlab) + xlim(xlim)
}

par(mar=c(5.1,8,4.1,2.1) )
pdf("../../../../../figure_materials/Chapter_1/Figure_2/Number_of_active_cells_WT_7_dpf.pdf")
barplot(unlist(snp[[1]]), names.arg = c(1,2,3,4,5,6,7,8), xlab ="Fish no.", ylab = "Number of Active Cells", cex.axis = 2, cex.names = 2, cex.lab = 2)
dev.off()

pdf("../../../../../figure_materials/Chapter_1/Figure_2/Response_duration_density_WT_7_dpf.pdf")
create_density_plot(snp [[2]]$response_dur, ylim = c(0,.7), xlab = "Response Duration (Seconds)", xlim = c(1,15))
dev.off()

pdf("../../../../../figure_materials/Chapter_1/Figure_2/Firing_frequency_density_WT_7_dpf.pdf")
create_density_plot(con = snp [[2]]$Firing_freq, ylim = c(0,55), xlim = c(0,.15), xlab = "Firing Frequency (Events/Second)")
dev.off()

pdf("../../../../../figure_materials/Chapter_1/Figure_2/Response_amplitude_density_WT_7_dpf.pdf")
create_density_plot(snp [[2]]$response_amp, ylim = c(0,5), xlab = "Response Amplitude (DF/F)", xlim = c(0.4,2))
dev.off()

# Example traces for diagram

c <- as.matrix(fread("180220_WT_h2b_gc6s_7dpf_f1_sa_aligned_all_cells_cal.dat"))
pdf("../../../../../figure_materials/Chapter_1/Figure_2/example_calcium trac.pdf")
plot(c [7,16500:16900], type ='l', ylim = c(0,0.9), lwd = 8)
dev.off()



###############################################################################################################
#                                          Figure 3
###############################################################################################################

raster_plot <- function(stmp1, m = "WT_NORM") {
  # stmp1 is the binary activity matrix
  circular_permutation <- function(t) {
    permute_by <- round(runif(1,1, length(t)))
    pt <- c(tail(t, -permute_by), head(t, permute_by))
    return(pt)
  }
  stmp1<- as.matrix(stmp1)
  stmp1_vec = as.numeric(unlist(stmp1))
  u1=rep(1:ncol(stmp1),1,each=nrow(stmp1))
  v1=rep(1:nrow(stmp1),ncol(stmp1))
  layout(matrix(c(1,1,1,1,1,1,2,2), nrow = 4, ncol = 2, byrow = TRUE))
  par(mar = c(0,6,5,2))
  plot(u1[stmp1_vec>0],v1[stmp1_vec>0],pch=19,cex=0.00005,col="black",cex.main=3, xaxt='n', xaxs = "i", yaxs ="i", ylim = c(0,nrow(stmp1)), ylab =  "Cell ID", cex.axis = 2, cex.lab = 3)
  title(main = list( m, cex = 2.5), adj = 0.5, line = 2)
  par(mar = c(6,6,0,2))
  plot((seq(1,ncol(stmp1))/4.85/60), colSums(stmp1), type ="h", axes=T, xaxs = "i", yaxs ="i", xlab = "Time (mins)", ylim = c(0,100), cex.lab=2.5, cex.axis=2, ylab = "Cell #" )
  abline(h = quantile(colSums(t(apply(stmp1, 1, circular_permutation))), .95), col ="red", lwd = 3)
}

save_rasters <- function(folder) {
  s_names <- list.files(path = folder, pattern = "all_cells_spikes.dat")
  for (i in 1:length(s_names)) {
    s <- as.matrix(fread(s_names [[i]]))
    png(paste("../../../../../figure_materials/Chapter_1/Figure_3/", s_names[[i]], "_Raster.png", sep = ""), width =3000, height = 1600)
    raster_plot(s, m = paste(s_names [[i]]))
    dev.off()
  }
}

save_rasters("./")
######################################################################################################################################


guassian_smooth <- function(x, width){
  if (sum(x) == 0) {
    sm_tr <- rep(0, length(x))
  }
  else{
    sm_tr = density(which(x == 1), width = width,  from = 0, to = length(x), n = length(x))$y
  }
  return(sm_tr)
}



correlations <- function(c) {
  c <- c[rowSums(c) > 0, ] 
  cor_c <- cor(t(c))
  cor_c <- as.vector(cor_c [lower.tri(cor_c)])
  shuffled_c <- t(apply(c, 1, circular_permutation))
  cor_c_sh <- cor(t(shuffled_c))
  cor_c_sh <- as.vector(cor_c_sh [lower.tri(cor_c_sh)])
  df <- data.frame(cbind(c(rep(0, length(cor_c)), rep(1, length(cor_c_sh))), c(cor_c, cor_c_sh)))
  return(df)
}



plot_correlations <- function(df, s_names){
  real_data <- hist(df$X2 [df$X1 == 0], breaks = seq(-1 ,1,0.02), plot=FALSE)
  real_data$counts=real_data$counts/sum(real_data$counts)
  plot(real_data)
  cp <- hist(df$X2 [df$X1 == 1], seq(-1 ,1,0.02), plot=FALSE)
  cp$counts=cp$counts/sum(cp$counts)
  plot(real_data$mids [35:100], real_data$counts [35:100], col = rgb(1,0,0, alpha = 1), log ="y", type ="l", lwd =4, xlab = "Correlation Coeff", ylab = "probablity")
  lines(cp$mids, cp$counts, col = rgb(0,0,1, alpha = 0.6), lwd =4)
  legend(x = 0.5, y = 1e-02, legend = c("data", "surrogate shuffled"), fill =  c("red", "blue"))
  title(s_names)
}


save_corr_plots <- function(folder) { 
  s_names <- list.files(path = folder, pattern = "all_cells_spikes.dat")
  for (i in 1:length(s_names)) {
    s <- as.matrix(fread(paste(folder, s_names [[i]], sep ="")))
    smooth_S <- t(apply(s, 1, guassian_smooth, 3))
    df <- correlations(smooth_S)
    real_data <- hist(df$X2 [df$X1 == 0], breaks = seq(-1 ,1,0.02), plot=FALSE)
    real_data$counts=real_data$counts/sum(real_data$counts)
    plot(real_data)
    cp <- hist(df$X2 [df$X1 == 1], seq(-1 ,1,0.02), plot=FALSE)
    cp$counts=cp$counts/sum(cp$counts)
    pdf(paste("../../../../../figure_materials/Chapter_1/Figure_3/", s_names [[i]], "_correlation_plot.pdf", sep =""))
  
    plot(real_data$mids [35:100], real_data$counts [35:100], col = rgb(1,0,0, alpha = 1), log ="y", type ="l", lwd =4, xlab = "Correlation Coeff", ylab = "probablity")
    lines(cp$mids, cp$counts, col = rgb(0,0,1, alpha = 0.6), lwd =4)
    legend(x = 0.5, y = 1e-02, legend = c("data", "surrogate shuffled"), fill =  c("red", "blue"))
    title(s_names [[i]])
    dev.off()
  }
}



get_MI_sim <- function(mat) {
  dice_coefficent <- function(a, b) {
    if (sum(a) == 0 || sum(b) == 0){
      return(0)
    }
    else {
      return(2*(a %*% b /(sum(a)+sum(b))))
    }
  }
  return(PairApply(t(mat), dice_coefficent, symmetric = TRUE))
}
  

sync_vs_non_sync_MI <- function(s) {
  sig_thresh <- quantile(colSums(t(apply(s, 1, circular_permutation))), .95)
  sync <- s [, colSums(s) > sig_thresh]
  non_syn <- s [,colSums(s) <  sig_thresh]
  MI_sync <- get_MI_sim(t(sync))
  MI_ns <- get_MI_sim(t(non_syn [, 2000:5000]))
  return(list(MI_sync, MI_ns))
}

all_sync_vs_non__sync_MI <- function(folder = "./"){
  s_names <- list.files(path = folder, pattern = "all_cells_spikes.dat")
  sync_MI_mean <- vector("numeric", length(s_names))
  non_sync_MI_mean <- vector("numeric", length(s_names))
  for (i in 1:length(s_names)){
    print(i)
    s <- as.matrix(fread(paste(folder, s_names [[i]], sep ="")))

    MI <- sync_vs_non_sync_MI(s)
    tiff(paste("../../../../../figure_materials/Chapter_1/Figure_3/", s_names [[i]], "MI_sim_mat_sych.tif", sep = ""), height = 2500, width =2705)
    image.plot(MI [[1]][200:ncol(MI[[1]]),200:ncol(MI[[1]])], col = viridis(200), zlim = c(0,0.6))
    dev.off()
  
    
    tiff(paste("../../../../../figure_materials/Chapter_1/Figure_3/", s_names [[i]], "_MI_sim_mat_non_sych.tif", sep = ""), height = 2505, width =2720)
    image.plot(MI [[2]][200:ncol(MI[[2]]),200:ncol(MI[[2]])], col = viridis(200), zlim = c(0,0.6))
    dev.off()
    sync_MI_mean [i] <- mean(MI [[1]][lower.tri(MI [[1]])])
    non_sync_MI_mean [i] <- mean(MI [[2]][lower.tri(MI [[2]])])
  }
 write.table(cbind(sync_MI_mean, non_sync_MI_mean), "MI_mean_sync_vs_non_sync.dataframe")
 return(cbind(sync_MI_mean, non_sync_MI_mean))

}

MI <- all_sync_vs_non__sync_MI("./")

############################## distance dependent correlations ############################################################################
midline <- data.frame(cbind(intercept = c(295, 312, 267, 345, 386, 327, 243, 362), slope = c(-1.23, -1.23, -1, -1.57, -1.828, -1.38, -1.012, -1.62)))





corr_by_dist <- function(cells, cal) {
  ecul <- as.vector(dist(cells,  method = "euclidean"))
  corrs <- as.vector(as.dist(cor(t(cal))))
  rand_corrs <- matrix(NA, length(corrs), ncol = 200)
  pb <- txtProgressBar(min = 0, max = ncol(rand_corrs), style = 3)
  for (i in 1:ncol(rand_corrs)) {
    rand_corrs [,i] <- as.vector(as.dist(cor(apply(cal, 1, circular_permutation))))
    setTxtProgressBar(pb, i)
  }
  sig_level <- apply(rand_corrs, 1, quantile, 0.95)
  rand_sd <- apply(rand_corrs, 1, sd)
  rand_mean <- apply(rand_corrs, 1, mean)
  plot(ecul, corrs ,pch=19, cex=.2,col=rgb(0,0,0,.1))
  #plot(lowess(ecul,corrs), type ="l", col = "red", lwd = 2)
  df <- data.frame(cbind(ecul, corrs, sig_level, rand_sd, rand_mean))
  return(df)
}



tect_hemisphere <- function(centers, cal, hemi = "l", mid = midline [6,]) {
  if (hemi == "l"){
    th <- centers [centers [,2] > (mid$slope *centers [,1] +mid$intercept),]
    resp <- cal [centers [,2] > (mid$slope *centers [,1] +mid$intercept),]
  }  
  if (hemi == "r") {
    th <- centers [centers [,2] < (mid$slope *centers [,1] +mid$intercept),]
    resp <- cal [centers [,2] < (mid$slope *centers [,1] +mid$intercept),]
  }
  return(list(th, resp))
} 

corr_by_dist_th <- function(centers, cal, midline_no =1 ) {
  aligned <- centers [rowSums(cal)> 0,]
  cal <- cal [rowSums(cal)> 0,]
  left_tect <- tect_hemisphere(aligned, cal, hemi = "l", midline [midline_no,])
  right_tect <- tect_hemisphere(aligned, cal, hemi = "r", midline [midline_no,])
  corr_dist_l <- corr_by_dist(left_tect [[1]], left_tect [[2]])
  corr_dist_r <- corr_by_dist(right_tect [[1]], right_tect [[2]])
  return(list(corr_dist_r, corr_dist_l))
}

corr_by_dist_all_folders <- function(folder = "./"){
  s_names <- list.files(path = folder, pattern = "all_cells_spikes.dat")
  c_names <- list.files(path = folder, pattern = "_all_cells_centers_cut_ordered.dat")
  midline <- data.frame(cbind(intercept = c(295, 312, 267, 345, 386, 327, 243, 362), slope = c(-1.23, -1.23, -1, -1.57, -1.828, -1.38, -1.012, -1.62)))
  for (i in 1:length(s_names)){
    centers <- as.matrix(fread(paste(folder, c_names [[i]], sep ="")))
    s <- as.matrix(fread(paste(folder, s_names [[i]], sep ="")))
    centers <- centers [rowSums(s) > 0,]
    s <- s [rowSums(s) > 0,]
    Smooth_S <- t(apply(s, 1, guassian_smooth, 3))
    
    a <- corr_by_dist_th(centers = centers, cal = Smooth_S, midline_no = i)
    a<-do.call(rbind,a)
    S_a <- a [abs(a$corrs) > abs(a$sig_level)*2,]
    S_a <- S_a [S_a$corrs > 0.01,]
    
    
    tiff(paste("../../../../../figure_materials/Chapter_1/Figure_3/", s_names [[i]], "corrs_by_dist.tif", sep = ""), height = 1000, width = 1000)
    #par(mar = c(5,5,5,5))
    ggplot(S_a, aes(x = ecul*0.7, y= abs(corrs))) + geom_point(size = 0.5, alpha =0.3) + theme_classic(base_size = 70) + ylab("Correlation Coeff") + xlab("Distance (microns)") + theme(panel.grid = element_blank(), panel.border = element_blank()) + geom_smooth(se = TRUE, size = 2)
    dev.off()
    
    write.table(a, file = paste(s_names [[i]], "corr_by_dist.dat"))
  }
}
