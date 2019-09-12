# These are the functions to calculate assembly paramemters
library(splancs)
library(igraph)
library(viridis)
library(DescTools)
library(manipulate)
library(fields)

# functions to calculate the size of the tectum

guassian_smooth <- function(x, width){
  if (sum(x) == 0) {
    sm_tr <- rep(0, length(x))
  }
  else{
    sm_tr = density(which(x == 1), width = width,  from = 0, to = length(x), n = length(x))$y
  }
  return(sm_tr)
}

cha<-function(x,y){
  chull(x,y)->i
  return(areapl(cbind(x[i],y[i])))
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


total_area <- function(minprob = 0.99){
  areas = vector("numeric", length = length(as.numeric(names(tab)[sel])))
  densities = vector("numeric", length = length(as.numeric(names(tab)[sel])))
  counter = 1
  for(i in as.numeric(ensel)){
    areas [counter] = cha(centers[mem==i & probs>minprob,1], centers[mem==i & probs>minprob,2])
    densities [counter] = length(centers[mem==i & probs>minprob,1])/(areas [counter])
    counter = counter +1
  }
  return(data.frame(areas, densities))
}


#finding and ploting midline of tectum and finding the anterior or posterior in order to calculate assembly spatial parameters
find_midline <- function(centers){
  plot_points <- function(data, a, b){
    smoothScatter(data[,1],data[,2], transformation = function(x) x^.25)
    points(data[,1],data[,2], cex= 0.2)
    abline(a = a, b = b, col ="red", lwd =5)
  }
  manipulate(plot_points(centers [,1:2], a, b), a =slider(-100,300), b = slider(-3,-0.5))
}

plot_midline <- function() {
  smoothScatter(centers[,1],centers[,2], transformation = function(x) x^.25)
  points(centers[,1],centers[,2])
  abline(a = midline [prefix_ID,1], b= midline [prefix_ID,2], col = "red", lwd = 2 )
}

plot_all_centers <- function(folder){
  files <- list.files(path = folder, pattern = "registered")
  plot(read.table(paste(folder, files [[1]], sep  = "")) [,1:2])
  for (i in 1:length(files)) {
    points(read.table(paste(folder, files [[i]], sep  = "")) [,1:2])
  }
  
}

Lateral_index <- function(){
  LI <- vector("numeric", length = length(ensel))
  counter = 1
  for (i in as.numeric(ensel)){
    assembly_mem <- centers[mem==i & probs>.99, 1:2]
    LH = sum(assembly_mem [,2] > (midline$slope [prefix_ID] *assembly_mem [,1] +midline$intercept [prefix_ID]))
    RH <- sum(assembly_mem [,2] < (midline$slope [prefix_ID]*assembly_mem [,1] +midline$intercept [prefix_ID]))
    ALL <- nrow(assembly_mem)
    LI [counter] <- ((LH-RH)/ALL)
    counter = counter + 1
  }
  return(LI)
}


LI <- function(mid = midline [1,]) {
  # X the omega_traj_object
  LI <- vector("numeric", length = length(ensel))
  counter = 1
  for (i in as.numeric(ensel)){
    assembly_mem <- centers [mem==i & probs>.99, 1:2]
    if (is.null(dim(assembly_mem))) {
      LI [i] <- NA
    }
    else {
      LH = sum(assembly_mem [,2] > (mid$slope *assembly_mem [,1] +mid$intercept))
      RH <- sum(assembly_mem [,2] < (mid$slope*assembly_mem [,1] +mid$intercept))
      ALL <- nrow(assembly_mem)
      LI [counter] <- ((LH-RH)/ALL)
      counter = counter +1
    }
  }
  return(LI)
}


NULL_LI <- function(mid = midline [1,], samples = 2000) {
  # X the omega_traj_object
  all_samples_LI <-list(length = samples)
  for (j in 1:samples){
    LI <- vector("numeric", length = length(ensel))
    counter = 1
    for (i in as.numeric(ensel)){
      assembly_mem <- centers [sample(mem==i & probs>.99),]
      if (is.null(dim(assembly_mem))) {
        LI [i] <- NA
      }
      else {
        LH = sum(assembly_mem [,2] > (mid$slope *assembly_mem [,1] +mid$intercept))
        RH <- sum(assembly_mem [,2] < (mid$slope*assembly_mem [,1] +mid$intercept))
        ALL <- nrow(assembly_mem)
        LI [counter] <- ((LH-RH)/ALL)
        counter = counter +1
      }
      }
      all_samples_LI [[j]] <- LI
    }
  return(as.matrix(do.call(cbind, all_samples_LI)))
}


rotate_tect_by_line <- function(centers = centers) {
  angle <- atan(midline$slope [prefix_ID]) +pi
  rot_angle <- angle - (pi/2)
  trans_mat <- matrix(c(cos(rot_angle), sin(rot_angle), -sin(rot_angle), cos(rot_angle)), 2, 2)
  rot_cent <- as.matrix(centers [,1:2]) %*% trans_mat
  a = seq(1,100,1)
  b = midline$slope[1]*seq(1,100,1)+midline$intercept[1]
  new_inter =  cbind(a, b) %*% trans_mat
  #shift_cent <- mid$intercept [prefix_ID]* cos(rot_angle)
  #print(shift_cent)
  rot_cent [,1] <- rot_cent [,1] - new_inter [1,1]
  return(rot_cent)
}


mean_by_bin <- function(centers) {
  new_cent <- rotate_tect_by_line(centers = centers )
  RL <- list(new_cent [new_cent [,1] > 0,], new_cent [new_cent [,1] < 0,])
  theta = c(pi/-3, pi/3)
  rotRL <- RL
  means <- list("numeric", length(theta))
  line <- list("numeric", length(theta))
  for (i in 1:length(theta)) {
    rotRL [[i]] <- cbind(Rotate(RL [[i]], theta = theta [i])$x, Rotate(RL [[i]], theta = theta [i])$y)
    cuts <- cut(x = rotRL [[i]] [,1], breaks = 10)
    means [[i]] <- cbind(by(data = rotRL [[i]] [,1], INDICES = cuts, FUN =function(x) {mean(x)}), by(data = rotRL [[i]] [,2], INDICES = cuts, FUN =function(x) {mean(x)}))
    line [[i]] <- cbind(Rotate(approx(means [[i]], n = 100), theta = theta [-i])$x,Rotate(approx(means [[i]], n = 100), theta = theta [-i])$y)
    
  }
  return(list(new_cent, line, RL))
} 




#plot raster plot of all neurons and the timepoints where there is sig sync activity - compeared to shuffled data
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
  plot(u1[stmp1_vec>0],v1[stmp1_vec>0],pch=19,cex=0.0001,col="black",cex.main=3, xaxt='n', xaxs = "i", yaxs ="i", ylim = c(0,nrow(stmp1)), ylab =  "Cell ID", cex.axis = 2, cex.lab = 3)
  title(main = list( m, cex = 2.5), adj = 0.5, line = 2)
  par(mar = c(6,6,0,2))
  plot((seq(1,ncol(s))/4.85/60), colSums(stmp1), type ="h", axes=T, xaxs = "i", yaxs ="i", xlab = "Time (mins)", ylim = c(0,100), cex.lab=2.5, cex.axis=2, ylab = "Cell #" )
  abline(h = quantile(colSums(t(apply(stmp1, 1, circular_permutation))), .95), col ="red", lwd = 3)
}

eigen_decomp_of_cov <- function() {
  whole_tectal_area <- ((eigen(cov(centers [,1:2]))$values [1]*eigen(cov(centers [,1:2]))$values [2])*pi)
  counter = 1
  norm_eigen_decomp <- vector("numeric", length = length(ensel))
  for (i in as.numeric(ensel)) {
    assembly_mem <- centers [mem==i & probs>.99,1:2]
    norm_eigen_decomp [counter] <- (eigen(cov(assembly_mem))$values [1])*(eigen(cov(assembly_mem))$values [2])*(pi)
    counter = counter +1
  }
  return(norm_eigen_decomp)
}

NULL_eigen_decomp_of_cov <- function(samples = 2000) {
  NULLS <- list(length = samples)
  for (j in 1:samples) {
    whole_tectal_area <- ((eigen(cov(centers [,1:2]))$values [1]*eigen(cov(centers [,1:2]))$values [2])*pi)
    counter = 1
    norm_eigen_decomp <- vector("numeric", length = length(ensel))
    for (i in as.numeric(ensel)) {
      a <- sum(mem==i & probs>.99)
      assembly_mem <- centers [sample(seq(1,nrow(centers)), a),1:2]
      norm_eigen_decomp [counter] <- (eigen(cov(assembly_mem))$values [1])*(eigen(cov(assembly_mem))$values [2])*(pi)
      counter = counter +1
      }
    NULLS [[j]] <- norm_eigen_decomp
  }
  return(as.matrix(do.call(cbind, NULLS)))
} 

within_cluster_corr <- function(){
  ensemble_mean <- vector("numeric", length = length(ensel))
  counter = 1
  spikes <- t(apply(s [selection +1,], 1, guassian_smooth, 3))
  for (i in as.numeric(ensel)){
    
    cor_mat <- cor(t(spikes[mem==i & probs>.99,]))
    lower <- cor_mat[lower.tri(cor_mat)]
    ensemble_mean [counter] <- mean(lower)
    counter = counter+1
  }
  return(ensemble_mean)
}

NULL_within_cluster_corr <- function(samples) {
  NULLS <- list(length = samples)
  spikes <- t(apply(s [selection +1,], 1, guassian_smooth, 3))
  for (j in 1:samples){
    ensemble_mean <- vector("numeric", length = length(ensel))
    counter = 1
    for (i in as.numeric(ensel)){
      a <- sum(mem==i & probs>.99)
      cor_mat <- cor(t(spikes[sample(seq(1,nrow(centers)), a),]))
      lower <- cor_mat[lower.tri(cor_mat)]
      ensemble_mean [counter] <- mean(lower)
      counter = counter+1
    }
  NULLS [[j]] <- ensemble_mean
  }
  return(as.matrix(do.call(cbind, NULLS)))
}



sort_by_correlation <-function(omega) {
  d <- dist(omega)
  
  ensel_sort <- ensel [hclust(d)$order]
  cell_list <- list(length(ensel_sort))
  for (i in 1:length(ensel_sort)) {
    cells <- s[selection+1, ] [mem == ensel_sort [i] &  probs > 0.99,]
    dis <- dist(cells)
    cells <- cells [hclust(dis)$order, ]
    cell_list [[i]] <- cells
  }
  raster_plot(do.call("rbind", cell_list))
}

# funciton to identify ensembles that are fireing at the begining

beg_fire <- function(omega_ext, bin_size){
  beg_rate <- rowSums(omega_ext [,1:bin_size]>0.2)
  av_rate <- vector(length = nrow(omega_ext))
  for (i in 1:nrow(omega_ext)) {
    trace = omega_ext [i,] [10061:ncol(omega_ext)] > 0.2
    reshaped_end = matrix(trace, ncol = bin_size)
    av_rate [i] <- max(rowSums(reshaped_end))
  }
  return((beg_rate-av_rate)/50)
}

########################################################################################################################

# Assembly params
plot_par <- function(P,plot=FALSE){
  
  assembly_centers<-matrix(NA,length(enIDs),3)
  pmu=as.matrix(read.table(paste(folder,"/pmu.dat",sep="")))[(tot_samples-nsamples+1):tot_samples,enIDs]
  lambda0=as.matrix(read.table(paste(folder,"/lambda0.dat",sep="")))[(tot_samples-nsamples+1):tot_samples,enIDs]
  lambda1=as.matrix(read.table(paste(folder,"/lambda1.dat",sep="")))[(tot_samples-nsamples+1):tot_samples,enIDs]
  gs=as.matrix(read.table(paste(folder,"/n.dat",sep="")))[(tot_samples-nsamples+1):tot_samples,enIDs]
  
  counter=0;
  for(i in enIDs){
    counter=counter+1;
    assembly_centers[counter,] <- apply(centers[mem==i,],2,mean);
  }
  
  if(plot){
    X11();layout(matrix(1:4,1,4))
    boxplot(pmu[,sel],main="pmu")
    boxplot(lambda0[,sel],main="lambda0")
    boxplot(lambda1[,sel],main="lambda1")
    boxplot(gs[,sel],main="group sizes")
  }
  
  return(list(pmu.mean=apply(pmu,2,mean),
              lambda0.mean=apply(lambda0,2,mean),
              lambda1.mean=apply(lambda1,2,mean),
              centers=assembly_centers,
              gs.mean=apply(gs,2,mean),
              membership=mem))
  
}

plot_F <- function(){
  plot(apply(FreeEnergy,1,mean),type='l')
}

plot_P <- function(){
  ptrace<-read.table(paste(folder,"/P.dat",sep=""))$V1
  ptrace<-ptrace[(length(ptrace)-nsamples):length(ptrace)]
  plot(ptrace,type='l')
}

plot_assemblies_simple <- function(minprob=0){
  layout(matrix(1:35,nrow=5))
  par(mar=c(0,0,0,0));
  for(i in 1:35){
    plot(centers[,1:2],pch=19,cex=1,xaxt='n',yaxt='n')
    title(paste(i),line=-1)
    points(centers[mem==i & probs>minprob,1:2],col="red",pch=19)
  }
  plot(apply(FreeEnergy,1,mean),type='l')
}


plot_assemblies_custom <- function(print=F,snc=FALSE){
  
  subnet_mem<-wc$membership
  minprob=0.99
  
  if(print) png(paste("figures2/assemblies_",prefix,".png",sep=""),width=6*2,height=7*2,units='in',res=300)
  layout(matrix(1:42,ncol=6,byrow=T))
  par(mar=rep(0.3,4));
  counter=0;
  
  for(i in ensel[order(subnet_mem)]){
    counter=counter+1;			
    plot(centers[,1:2],pch=18,cex=.3,xaxt='n',yaxt='n',col=rgb(0,0,0,1))
    text(x=-100,y=100,paste(i),cex=2.5)
    #text(x=-100,y=100,i,cex=1.5)
    points(centers[mem==i & probs>minprob,1:2],col=viridis(max(subnet_mem))[sort(subnet_mem)[counter]],pch=19)
  }
  if(print) dev.off();
}

plot_cormat<- function(print=F){
  if(print) png(paste("figures2/cormat_",prefix,".png",sep=""),width=5*2,height=5*2,units='in',res=300)
  
  subnet_mem=wc$membership
  layout(matrix(1:4,2,2),widths=c(.5,5),heights=c(5,.5))
  par(mar=c(0,1,1,0))
  image(matrix(1:sum(sel),1,sum(sel)),col=viridis(max(subnet_mem))[subnet_mem[order(subnet_mem)]],xlab='',ylab='',xaxt='n',yaxt='n')
  par(mar=c(1,1,0,0))
  plot.new()
  par(mar=c(0,0,1,1))
  cmat=cor(t(omega[order(wc$membership),]))
  image(cmat,col=grey.colors(100),frame=F,axes=F)
  par(mar=c(1,0,0,1))
  image(matrix(1:sum(sel),sum(sel),1),col=viridis(max(subnet_mem))[subnet_mem[order(subnet_mem)]],xlab='',ylab='',xaxt='n',yaxt='n')
}



plot_assemblies <- function(x,minprob=0){
  if(x<0){
    layout(matrix(1:55,nrow=5))
    par(mar=c(0,0,0,0));
    for(i in ensel){
      
      plot(centers[,1:2],pch=18,cex=.3,xaxt='n',yaxt='n',col=rgb(0,0,0,1))
      title(paste(i),line=-1)
      points(centers[mem==i & probs>minprob,1:2],col="red",pch=19)
    }
    plot(apply(FreeEnergy,1,mean),type='l')
  }
  
  if(x>0) {
    plot(centers[,1:2],col="black",pch=19)
    points(centers[mem==x & probs>minprob,1:2],col="red",pch=19)
  }
  
  
}

get_omega <- function(print=FALSE,P,PMAX=200,nav=nsamples){
  system(paste("tail -n ",PMAX*nav," ",folder,"/omega_traj.dat > ",folder,"/omega_traj_tmp.dat",sep=""))
  omegatraj<-as.matrix(read.table(paste(folder,"/omega_traj_tmp.dat",sep="")))
  
  times=ncol(omegatraj)
  samples=nrow(omegatraj)/PMAX
  
  omegaprob=matrix(0,ncol=times,nrow=PMAX)
  for(i in (samples-nav):(samples-1)){
    omegaprob=omegaprob+omegatraj[ i*PMAX + 1:PMAX,]
  }
  
  omegaprob=omegaprob/nav
  omegaprob=omegaprob[ensel,]
  return(omegaprob)
}

get_omega0 <- function(print=FALSE,P,PMAX=200,nav=nsamples){
  system(paste("tail -n ",PMAX*nav," ",folder,"/omega_traj.dat > ",folder,"/omega_traj_tmp.dat",sep=""))
  omegatraj<-as.matrix(read.table(paste(folder,"/omega_traj_tmp.dat",sep="")))
  
  times=ncol(omegatraj)
  samples=nrow(omegatraj)/PMAX
  
  omegaprob=matrix(0,ncol=times,nrow=PMAX)
  for(i in (samples-nav):(samples-1)){
    omegaprob=omegaprob+omegatraj[ i*PMAX + 1:PMAX,]
  }
  
  omegaprob=omegaprob/nav
  return(omegaprob)
}


plot_omega <- function(print=FALSE,PMAX=200,subset=NA,nav=nsamples){
  system(paste("tail -n ",PMAX*nav," ",folder,"/omega_traj.dat > ",folder,"/omega_traj_tmp.dat",sep=""))
  omegatraj<-as.matrix(read.table(paste(folder,"/omega_traj_tmp.dat",sep="")))
  
  times=ncol(omegatraj)
  samples=nrow(omegatraj)/PMAX
  
  omegaprob=matrix(0,ncol=times,nrow=PMAX)
  for(i in (samples-nav):(samples-1)){
    omegaprob=omegaprob+omegatraj[ i*PMAX + 1:PMAX,]
  }
  
  par(mar=c(1,4,1,1))
  omegaprob=omegaprob/nav
  if(is.na(subset)){
    omegaprob=omegaprob[ensel,]
    image(x=1:times,y=1:sum(sel),z=t(omegaprob),col=grey.colors(100),yaxt='n',xlab="",ylab="",xaxt='n',cex.axis=2,cex.lab=2,breaks=seq(0,1,length=101))
    axis(2,at=1:sum(sel),labels=ensel)
  } else {
    omegaprob=omegaprob[subset,]
    image(x=1:times,y=1:length(subset),z=t(omegaprob),col=grey.colors(100),yaxt='n',xlab="synchronous frames",ylab="assemblies",xaxt='n',cex.axis=2,cex.lab=2,breaks=seq(0,1,length=101))
    axis(2,at=1:length(subset),labels=subset)
  }
  
}

plot_omega_custom <- function(print=FALSE,P,PMAX=100,subset=NA){
  
  if(print) png(paste("figures2/omega_",prefix,".png",sep=""),width=5*2,height=7*2,units='in',res=300)
  
  subnet_mem=wc$membership
  times=ncol(omega)
  
  par(mar=c(1,4,1,0))
  omegaprob=as.matrix(omega[order(subnet_mem),])
  layout(matrix(1:2,1,2),widths=c(5,.3))
  
  image(x=1:times,y=1:sum(sel),z=1-t(omegaprob),col=grey.colors(100),yaxt='n',xlab="",ylab="",xaxt='n',cex.axis=2,cex.lab=2,breaks=seq(0,1,length=101))
  
  #axis(2,at=1:sum(sel),labels=rep("",sum(sel)))
  axis(2,at=1:sum(sel),labels=ensel[order(subnet_mem)],las=2,cex.axis=2)
  
  par(mar=c(1,0,1,1))
  image(matrix(1:sum(sel),1,sum(sel)),col=viridis(max(subnet_mem))[subnet_mem[order(subnet_mem)]],xlab='',ylab='',xaxt='n',yaxt='n')
  if(print) dev.off()
  
}

autocor<- function(s){
  ls=length(s)
  tmat<-matrix(0,2,2)
  for(i in 1:(ls-1)){
    tmat[s[i+1]+1,s[i]+1]=tmat[s[i+1]+1,s[i]+1]+1;
  }
  tmat = t(t(tmat)/colSums(tmat))
  return(tmat[2,2])
}

cordist<- function(use="L"){
  
  if(use=="L"){
    sel_side=(side_AP$enID %in% ensel) & side_AP$RL==1 & !is.na(side_AP$RL);
  } else {
    sel_side=(side_AP$enID %in% ensel) & side_AP$RL==0 & !is.na(side_AP$RL);
  }
  
  ng=sum(sel_side)
  
  assembly_centers<-matrix(NA,ng,2);
  counter=0;
  for(i in side_AP$enID[sel_side]){
    counter=counter+1;
    assembly_centers[counter,]<-apply(centers[mem==i,1:2],2,mean)
  }
  omegamat=get_omega0()[side_AP$enID[sel_side],];
  cormat=cor(t(omegamat))[lower.tri(matrix(NA,ng,ng))];
  dist=as.matrix(dist(assembly_centers,upper=T))[lower.tri(matrix(NA,ng,ng))]
  
  return(list(cormat=cormat,dist=dist))
  
}

get_corr_AP <- function(print=F){
  
  cormat=cor(t(get_omega()));
  diag(cormat)=0;
  sel_A=(side_AP_ord$AP==0 & !is.na(side_AP_ord$AP))
  sel_P=(side_AP_ord$AP==1 & !is.na(side_AP_ord$AP))
  
  # use Left side 
  sel_side=(side_AP_ord$RL==1 & !is.na(side_AP_ord$RL));
  cormat_APL=as.vector(cormat[sel_side & sel_A, sel_side & sel_P]);
  cormat_AAL=as.vector(cormat[sel_side & sel_A, sel_side & sel_A]);
  cormat_PPL=as.vector(cormat[sel_side & sel_P, sel_side & sel_P]);
  
  # use right side 
  sel_side=(side_AP_ord$RL==0 & !is.na(side_AP_ord$RL));
  cormat_APR=as.vector(cormat[sel_side & sel_A, sel_side & sel_P]);
  cormat_AAR=as.vector(cormat[sel_side & sel_A, sel_side & sel_A]);
  cormat_PPR=as.vector(cormat[sel_side & sel_P, sel_side & sel_P]);
  if(print){
    write.table(c(cormat_APR,cormat_APL),file=paste(folder,"/cormatAP.dat",sep=""),col.names=F,row.names=F) 
    write.table(c(cormat_AAR,cormat_AAL),file=paste(folder,"/cormatAA.dat",sep=""),col.names=F,row.names=F) 
    write.table(c(cormat_PPR,cormat_PPL),file=paste(folder,"/cormatPP.dat",sep=""),col.names=F,row.names=F) 
  }
  
  return(list(AP=c(cormat_APR,cormat_APL),AA=c(cormat_AAR,cormat_AAL), PP=c(cormat_PPR,cormat_PPL)))
  
}

get_corr_AP_across <- function(print=F){
  
  cormat=cor(t(get_omega()));
  diag(cormat)=0;
  sel_A=(side_AP_ord$AP==0 & !is.na(side_AP_ord$AP))
  sel_P=(side_AP_ord$AP==1 & !is.na(side_AP_ord$AP))
  left_side=(side_AP_ord$RL==1 & !is.na(side_AP_ord$RL));
  right_side=(side_AP_ord$RL==0 & !is.na(side_AP_ord$RL));
  
  # use Left/right side 
  cormat_APL=as.vector(cormat[right_side & sel_A, left_side & sel_P]);
  cormat_AAL=as.vector(cormat[right_side & sel_A, left_side & sel_A]);
  cormat_PPL=as.vector(cormat[right_side & sel_P, left_side & sel_P]);
  
  # use right side 
  cormat_APR=as.vector(cormat[left_side & sel_A, right_side & sel_P]);
  cormat_AAR=as.vector(cormat[left_side & sel_A, right_side & sel_A]);
  cormat_PPR=as.vector(cormat[left_side & sel_P, right_side & sel_P]);
  if(print){
    write.table(c(cormat_APR,cormat_APL),file=paste(folder,"/cormatAP_across.dat",sep=""),col.names=F,row.names=F) 
    write.table(c(cormat_AAR,cormat_AAL),file=paste(folder,"/cormatAA_across.dat",sep=""),col.names=F,row.names=F) 
    write.table(c(cormat_PPR,cormat_PPL),file=paste(folder,"/cormatPP_across.dat",sep=""),col.names=F,row.names=F) 
  }
  
  return(list(AP=c(cormat_APR,cormat_APL),AA=c(cormat_AAR,cormat_AAL), PP=c(cormat_PPR,cormat_PPL)))
  
}

plot_corr_AP <- function(print=FALSE){
  u<-get_corr_AP(print);
  ind=seq(0,max(c(u$AA,u$PP)),0.01)
  if(print) pdf(paste("figures2/AP_correlations_",prefix,".pdf",sep=""),5,5)
  plot(ind,1-ecdf(u$AA)(ind),type='l',col="red",xlab="time correlation",ylab="frequency",cex.lab=1.5,cex.axis=1.5)
  lines(ind,1-ecdf(u$PP)(ind),col="blue")
  lines(ind,1-ecdf(u$AP)(ind),col="magenta")
  legend("topright",legend=c("AA","PP","AP"),lty=1,col=c("red","blue","magenta"))
  if(print) dev.off();
}

plot_corr_AP_across <- function(print=FALSE){
  u<-get_corr_AP_across(print);
  ind=seq(0,max(c(u$AA,u$PP)),0.01)
  if(print) pdf(paste("figures2/AP_correlations_across",prefix,".pdf",sep=""),5,5)
  plot(ind,1-ecdf(u$AA)(ind),type='l',col="red",xlab="time correlation",ylab="frequency",cex.lab=1.5,cex.axis=1.5)
  lines(ind,1-ecdf(u$PP)(ind),col="blue")
  lines(ind,1-ecdf(u$AP)(ind),col="magenta")
  legend("topright",legend=c("AA","PP","AP"),lty=1,col=c("red","blue","magenta"))
  if(print) dev.off();
}

plot_corr <- function(print=FALSE){
  uL=cordist("L")
  uR=cordist("R")
  if(print){
    write.table(cbind(c(uL$dist,uR$dist),c(uL$cormat,uR$cormat)),file=paste("figures2/corr_vs_dist_",prefix,".dat",sep=""),col.names=F,row.names=F)
  }
  if(print) pdf(paste("figures2/corr_vs_dist_",prefix,".pdf",sep=""),height=5,width=10)
  layout(matrix(1:2,1,2))
  plot(uL$dist,uL$cormat,main="OTL",xlab="distance",ylab="correlation",cex.lab=1.5)
  plot(uR$dist,uR$cormat,main="OTR",xlab="distance",ylab="correlation",cex.lab=1.5)
  if(print) dev.off();
}

plot_par_means<- function(pmu_thresh=0,gs_thresh=1,xlim=c(0,1),ylim=c(0,1)) {
  plot(A$lambda1.mean[A$pmu.mean>pmu_thresh & A$gs.mean>gs_thresh],A$lambda0.mean[A$pmu.mean>pmu_thresh & A$gs.mean>gs_thresh],cex=2,cex.axis=3,cex.lab=2,xlim=xlim,ylim=ylim)
}


plot_cells <- function(i,thresh=0){
  layout(matrix(1:2,2,1),heights=c(7,1))
  id=which(ensel==i)
  sum_act=apply(s[selection+1,],2,sum)
  par(mar=c(0,2,2,1))
  image(t(as.matrix(s)[selection+1,sum_act>thresh][mem==i & probs>0.9,]),axes=F,frame=F)
  par(mar=c(2,2,0,1))
  image(as.matrix(omega[id,]),frame=F,axes=F,col=gray.colors(100),breaks=seq(0,1,length=101))
}

plot_s <- function(thresh=15){
  mem2=mem;
  cls=as.numeric(gsub("V","",names(A$gs.mean)))[sel]
  for(i in 1:length(cls)){
    mem2[mem==cls[i]]=i
  }
  
  mem2[!mem%in%cls]=-1
  sum_act=apply(s[selection+1,],2,sum)
  ord2=order(mem2[mem2!=-1])
  
  layout(matrix(1:2,2,1))
  par(mar=c(0.1,0.1,0.1,0.1));
  image(1-t(as.matrix(s)[selection[mem2!=-1][ord2]+1,sum_act>thresh]),frame=F,xaxt='n',yaxt='n',col=gray.colors(100))
  image(1-t(as.matrix(s)[selection[mem2==-1]+1,sum_act>thresh]),frame=F,xaxt='n',yaxt='n',col=gray.colors(100))
}

plot2_s <- function(print=F,width=200){
  if(print) png(paste("figures2/",prefix,"_smatrix.png",sep=""),15*width,5*200) 
  stmp1=as.matrix(s)
  stmp1_vec = as.numeric(stmp1)
  u1=rep(1:ncol(stmp1),1,each=nrow(stmp1))
  v1=rep(1:nrow(stmp1),ncol(stmp1))
  plot(u1[stmp1_vec>0],v1[stmp1_vec>0],pch=19,cex=0.1,col=viridis(10)[4],axes=F,xlab="",ylab="")
  if(print) dev.off();
}

plot2_cells_beh <- function(i,b="rs",print=F,width=200,len=NA,
                            col=viridis(10)[4]){
  if(print) png(paste("figures2/",prefix,"_enIDs_",i,".png",sep=""),15*width,5*200*11/19.)
  layout(matrix(1:3,3,1),heights=c(5,1,3))
  
  if(is.na(len)) len=ncol(s)
  btraj=colMeans(s[selection+1,][mem==i & probs>0.99,]);
  
  tabsel=as.numeric(names(sort(table(mem),decreasing=T))[sel])
  enID=which(tabsel==i)
  stmp1=as.matrix(s)[selection+1,][mem==i & probs>0.99,]
  stmp1_vec=as.numeric(stmp1)
  u1=rep(1:ncol(stmp1),1,each=nrow(stmp1))
  v1=rep(1:nrow(stmp1),ncol(stmp1))
  omega_row=omega_extended[enID,1:len]
  
  par(mar=c(0,5,2,1))
  image(xlim=c(0,len),ylim=c(1,nrow(stmp1)),z=matrix(NA,2,2),xlab="",ylab="",xaxt='n',yaxt='n')
  points(u1[stmp1_vec>0],v1[stmp1_vec>0],pch=19,cex=0.1,col=col)
  par(mar=c(0,5,0,1))
  image(xlim=c(0,len),ylim=c(0,1.2),z=matrix(NA,2,2),xlab="",ylab="",xaxt='n',yaxt='n')
  abline(v=which(omega_row>0.8)[1:len],col=gray.colors(100))
  par(mar=c(4,5,0,1))
  image(xlim=c(0,len),ylim=range(btraj[1:len]),z=matrix(NA,1,2),cex.axis=2,xlab="",ylab="",xaxt='n',las=2)
  lines(btraj[1:len])
  if(print) dev.off()
}

group_plots <- function(start){
  if(any(is.na(start))){
    X11(); plot_assemblies(-1,.99);
    X11(); plot_omega();
    X11(); P<-read.table(paste(folder,"/P.dat",sep=""))$V1; barplot(table(P)); 
  } else {
    dev.set(start[1]); plot_assemblies(-1,.99);
    dev.set(start[2]); plot_omega();
    dev.set(start[3]); P<-read.table(paste(folder,"/P.dat",sep=""))$V1; barplot(table(P));
  }
}

get_props <- function(x) return(c(A$pmu.mean[x],A$lambda0.mean[x],A$lambda1.mean[x]))

plot_subfig <- function(){
  plot(A$pmu.mean[A$gs.mean>gs_thresh],A$lambda1.mean[A$gs.mean>gs_thresh],
       pch=19,xlab="activity", ylab="coherence",cex.lab=1.2,cex.axis=1.2,cex.lab=1.5); 
  lines(x=rep(pmu_thresh,2), y=c(0,1),lty=2,col="green",lwd=2)
  lines(x=c(0,1), y=rep(lambda1_thresh,2),lty=2,col="green",lwd=2)
  
}

plot_graph <- function(print=F,th=0,snc=FALSE,nav=nsamples){
  corrmat=cor(t(get_omega(nav=nav)))
  corrmat[corrmat<th]=0
  diag(corrmat)=0
  if(snc){
    graph<-graph.adjacency(corrmat[order(subnet_col),order(subnet_col)],weighted=TRUE,mode="lower")
  } else {
    graph<-graph.adjacency(corrmat,weighted=TRUE,mode="lower")
  }
  
  # 1 is left, 2 is right
  shapes=c("circle","square")
  if(snc){
    V(graph)$color=viridis(3)[subnet_col[order(subnet_col)]]
    V(graph)$shape=shapes[side_AP_ord$RL[order(subnet_col)]+1]
  } else {
    V(graph)$color[!is.na(side_AP_ord$RL)]=c("#ff9500","#0059ff")[side_AP_ord$RL[!is.na(side_AP_ord$RL)]+1]
  }
  V(graph)$label=ensel
  
  if(print) pdf(paste(folder,"/graph.pdf",sep=""))
  plot(graph)
  if(print) dev.off();
}

plot_subnet <- function(print=F,nav=200,th=0.05){
  subnet_mem<-wc$membership
  
  # 1 is left, 2 is right
  shapes=c("circle","square")
  V(graph)$shape[!is.na(side_AP_ord$RL)]=shapes[side_AP_ord$RL[!is.na(side_AP_ord$RL)]+1]
  V(graph)$shape[is.na(side_AP_ord$RL)]="none"
  V(graph)$label.color="black"
  
  if(print) pdf(paste("figures2/subnets_",prefix,".pdf",sep=""))
  
  plot(wc,
       graph,
       mark.col=viridis(max(wc$membership),alpha=0.15), #shaded community colors
       mark.border=viridis(max(wc$membership)), # border colors of the communities
       col=(viridis(max(wc$membership),alpha=0.4)[wc$membership])) # vertex colors
  
  if(print) dev.off();
}

get_feature_df<- function(){
  return(data.frame(gs=A$gs.mean[sel],pmu=A$pmu.mean[sel],lambda0=A$lambda0.mean[sel],lambda1=A$lambda1.mean[sel]))
}

get_cormat_bayes<- function(nav=10,PMAX=200,th,disp=FALSE){
  system(paste("tail -n ",PMAX*nav," ",folder,"/omega_traj.dat > ",folder,"/omega_traj_tmp.dat",sep=""))
  
  corsel=matrix(0,length(ensel),length(ensel))
  cortab=array(0, dim=c(length(ensel),length(ensel),nav))
  omegatraj<-as.matrix(read.table(paste(folder,"/omega_traj_tmp.dat",sep="")))
  
  times=ncol(omegatraj)
  samples=nrow(omegatraj)/PMAX
  
  omegaprob=matrix(0,ncol=times,nrow=PMAX)
  counter=0;
  for(i in (samples-nav):(samples-1)){
    counter=counter+1
    omegasample=omegatraj[ i*PMAX + 1:PMAX,]
    omegasample=omegasample[ensel,]
    corsample=cor(t(omegasample))
    corsample[is.na(corsample)]=0
    cortab[,,counter]=corsample
    corsel[corsample>th]=corsel[corsample>th]+1
    if(disp) image(corsel,col=gray.colors(100))
  }
  
  cormean=apply(cortab,c(1,2),mean)
  ad=cormean;
  ad[corsel/nav<0.9]=0
  diag(corsel)=0
  diag(ad)=0
  
  graph <- graph.adjacency(ad,weighted=TRUE,mode="lower")
  V(graph)$color[!is.na(side_AP_ord$RL)]=c("#ff9500","#0059ff")[side_AP_ord$RL[!is.na(side_AP_ord$RL)]+1]
  V(graph)$label=ensel
  
  return(list(cortab=cortab,corsel=corsel/nav,ad=ad,graph=graph))
}


##########################################################################################
check_exclusions <- function(folder_names) {
  omega <- as.matrix(read.table(file = paste(folder_names, "/omega_be.dat", sep ='')))
  omega_ex <- as.matrix(read.table(file = paste(folder_names, "/omega_excluded.dat", sep ='')))
  layout(matrix(c(1,1,1,1,1,1,2,2), nrow = 4, ncol = 2, byrow = TRUE))
  par(mar = c(0,6,5,2))
  image(t(omega))
  par(mar = c(6,6,0,2))
  
  image(t(omega_ex))
}

plot_assembly_w_omega <- function(ensel_num, xlim = c(1,17460)){
  cells <- s[selection+1, ] [mem == ensel [ensel_num] &  probs > 0.99,]
  cor_with_cluster <- cor(omega_extended [which(ensel == ensel [ensel_num]),] > 0.8, t(cells))
  cells <- cells [order(cor_with_cluster,decreasing = TRUE), ]
  layout(matrix(c(1,1,1,1,1,1,2,2), nrow = 4, ncol = 2, byrow = TRUE))
  par( mar = c(0,6,6,2))
  raster_plot2(cells, xlim = xlim)
  par( mar = c(2,6,0,2))
  plot(omega_extended [which(ensel == ensel [ensel_num]),], type ='h', xlim = xlim, xaxs = "i", yaxs ="i",  cex.axis =2, ylab = "prob Aseembly Active")
}



raster_plot2 <- function(stmp1, m = "WT_GRAV",axes =F, xlim = c(0,17460)) {
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
  plot(u1[stmp1_vec>0],v1[stmp1_vec>0],pch=19,cex=0.001,col="black",cex.main=3, xlim = xlim, xaxt='n',  xaxs = "i", yaxs ="i", ylim = c(0,nrow(stmp1)), ylab =  "Cell ID", cex.axis = 2, cex.lab = 3)
  title(main = list( m, cex = 2.5), adj = 0.5, line = 2)
}


Assembly_centers <- function(centers) {
  x_coords <- vector("numeric", length(ensel))
  y_coords <- vector("numeric", length(ensel))
  z_coords <- vector("numeric", length(ensel))
  counter = 1
  for (i in as.numeric(ensel)) {
    assembly_mem <- centers[mem==i & probs>.99, ]
    x_coords [counter] <- mean(assembly_mem [,1])
    y_coords [counter] <- mean(assembly_mem [,2])
    z_coords [counter] <- mean(assembly_mem [,3])
    counter = counter +1
   }
  return(data.frame(x_coords, y_coords, z_coords))
}