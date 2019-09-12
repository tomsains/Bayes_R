# Chapter 1 - Assembly figures
library(data.table)
library(sm)
library(ggplot2)
library(viridis)
library(DescTools)
library(fields)
library(matrixStats)
library(ggridges)
library(data.table)
library(sm)
library(ggplot2)
library(viridis)
library(DescTools)
library(fields)
library(matrixStats)

setwd("~/../../media/thomas_sainsbury/Samsung_T5/SeG/results/Baysian_network_inference/data_for_Gibbs/7_dpf/")


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
  plot(u1[stmp1_vec>0],v1[stmp1_vec>0],pch=19,cex=0.001,col="black",cex.main=3, xlim = xlim, xaxt='n',  xaxs = "i", yaxs ="i", ylim = c(0,nrow(stmp1)), ylab =  "Cell ID", cex.axis = 3, cex.lab = 3)
  title(main = list( m, cex = 2.5), adj = 0.5, line = 2)
}


plot_raster_sorted <- function(s, fr = c(1,ncol(s))) {
  c <- s [,fr [1]:fr[2]]
  ff <- time_to_fire(c)
  raster_plot2(s [order(ff),])
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


wcc_plot_1 <- read.table(file = "../../../../figure_materials/Chapter_1/Figure_4/180510_WT_h2b_gc6s_7dpf_f3_wcc_vs_NULL.dat")



plot_correlations <- function(df, s_names){
  real_data <- hist(df$X2 [df$X1 == 0], breaks = seq(0 ,0.5,0.02), plot=FALSE)
  real_data$counts=real_data$counts/sum(real_data$counts)
  plot(real_data)
  cp <- hist(df$X2 [df$X1 == 1], seq(-1 ,1,0.02), plot=FALSE)
  cp$counts=cp$counts/sum(cp$counts)
  plot(real_data$mids [35:100], real_data$counts [35:100], col = rgb(1,0,0, alpha = 1), log ="y", type ="l", lwd =4, xlab = "Correlation Coeff", ylab = "probablity")
  lines(cp$mids, cp$counts, col = rgb(0,0,1, alpha = 0.6), lwd =4)
  legend(x = 0.5, y = 1e-02, legend = c("data", "surrogate shuffled"), fill =  c("red", "blue"))
  title(s_names)
}


  


p_0.01 <- quantile(wcc_plot_1$wcc_value [wcc_plot_1$sample == "NULL_wcc"], 0.99)
ggplot(wcc_plot_1 [wcc_plot_1$sample == "wcc", ], aes(x=wcc_value)) + geom_hline(yintercept = 0,lwd =1.2) + geom_line(color = "#762A83", alpha = 1, stat = "density", lwd= 2) +
  xlim(c(-0.05,0.5)) + 
  #geom_vline( = quantile(wcc_plot_1$wcc_value [wcc_plot_1$sample == "NULL_wcc"], 0.99), colour="red", linetype = "longdash") +
  geom_segment(aes(x =  p_0.01 , y = 0, xend =p_0.01, yend = Inf),col = "#D73C4C", linetype = "dashed", lwd = 1.2, alpha = 0.1)+
  theme_classic(base_size = 30) +
  geom_boxploth(data = wcc_plot_1 [wcc_plot_1$sample == "NULL_wcc",], aes(wcc_value, y = -0.5), size = 2,  outlier.size = 1, outlier.alpha = 0.01) + ylab("Probability Density") +
  xlab("Mean Ensemble Correlation Coeff.") 
ggsave("../../../../figure_materials/Chapter_1/Figure_4/Within_Cluster_Correlation.pdf")

(sum(wcc_plot_1$wcc_value [wcc_plot_1$sample == "wcc"] > p_0.01)/length(wcc_plot_1$wcc_value [wcc_plot_1$sample == "wcc"]))*100



centers <- as.matrix(read.table("../data_for_Gibbs/7_dpf/180510_WT_h2b_gc6s_7dpf_f1_sa_aligned_all_cells_centers_cut_ordered.dat"))
s<-as.matrix(read.table("../data_for_Gibbs/7_dpf/180510_WT_h2b_gc6s_7dpf_f1_sa_aligned_all_cells_spikes.dat"))
mem<-as.matrix(read.table("180510_WT_h2b_gc6s_7dpf_f1/mem.dat"))
probs<-as.matrix(read.table("180510_WT_h2b_gc6s_7dpf_f1/P.dat"))  

smooth_S <- t(apply(s, 1, guassian_smooth, 3))
corrs <- cor(t(smooth_S [rowSums(smooth_S) > 0,]))
which.max(rowSums(corrs))


library(RColorBrewer)
cols <- brewer.pal(9, "RdPu")
cols

tiff("../../../../figure_materials/Chapter_1/Figure_3/seed_map_1.tif", width = 1200,height= 1000)
diag(corrs) <- 0
c_cor <- corrs [4,]

ordered_cent <- data.frame(centers) [rowSums(s) > 0,] [order(c_cor),] 
c_cor <- sort(c_cor)
ggplot(ordered_cent, aes(V1, V2)) + geom_point(aes(color = c_cor), size =7) + scale_color_viridis(option = "D") +theme_classic(base_size  = 40) + geom_point(data = ordered_cent [which(c_cor == 0),], aes(V1, V2), col ='red', size = 10 )

dev.off()

tiff("../../../../figure_materials/Chapter_1/Figure_3/seed_map_2.tif", width = 1200,height= 1000)
diag(corrs) <- 0
c_cor <- corrs [7,]

ordered_cent <- data.frame(centers) [rowSums(s) > 0,] [order(c_cor),] 
c_cor <- sort(c_cor)
ggplot(ordered_cent, aes(V1, V2)) + geom_point(aes(color = c_cor), size =7) + scale_color_viridis(option = "D") +theme_classic(base_size  = 40) + geom_point(data = ordered_cent [which(c_cor == 0),], aes(V1, V2), col ='red', size = 10 )

dev.off()


omega <- as.matrix(read.table(file = "180510_WT_h2b_gc6s_7dpf_f1/omega.dat"))

tiff("../../../../figure_materials/Chapter_1/Figure_4/Assembly_activity_map.tif", width = 3000, height = 2000)
par(cex.axis=3.5)
image(t(omega), col = grey.colors(n = 100, start = 0.98,0.2,gamma = 0.5,alpha = 0.8), axes =F)
axis(2, at = seq(0,1,1/nrow(omega)), labels=seq(0,nrow(omega)),tick=TRUE, cex = 6)
dev.off()


ncol(omega)





vs_null_df <- function(folder ="../../../../figure_materials/Chapter_1/Figure_4/", return = "comp") {
  s_names <- list.files(path = "./", pattern = "_WT_h2b_gc6s_7dpf") 
  each_fish_list <- vector("list", length(s_names))
  for (i in 1:length(each_fish_list)) {
    print(i)
    LI <- read.table(paste("../../../../figure_materials/Chapter_1/Figure_4/", s_names [[i]], "_LI_vs_NULL.dat", sep =""))
    WCC <- read.table(paste("../../../../figure_materials/Chapter_1/Figure_4/", s_names [[i]], "_wcc_vs_NULL.dat", sep = ""))
    comp <- read.table(paste("../../../../figure_materials/Chapter_1/Figure_4/", s_names [[i]], "_compct_vs_NULL.dat", sep =""))
    if (return == "comp"){
      df <- as.data.frame(cbind(rep(s_names [[i]], nrow(comp)), rep(i, nrow(comp))))
      df <- cbind(df, comp)
      each_fish_list [[i]] <-df
    }
    if (return == "LI") {
      df <- as.data.frame(cbind(rep(s_names [[i]], nrow(LI)), rep(i, nrow(LI))))
      df <- cbind(df, LI)
      each_fish_list [[i]] <-df
    }
    if (return == "WCC") {
      df <- as.data.frame(cbind(rep(s_names [[i]], nrow(WCC)), rep(i, nrow(WCC))))
      df <- cbind(df, WCC)
      each_fish_list [[i]] <- df
    }
  }
  return(do.call(rbind, each_fish_list))
}


ggplot(df9, aes (x = comp_value)) +geom_density(aes(fill = sample), alpha = 0.5) +scale_fill_brewer(palette="Dark2") +theme_classic() + 
  xlim(c(0,7.0e+07)) + geom_vline(xintercept = quantile(df9$comp_value [df9$sample == "NULL_comp"], 0.05), color = "red")

folder=paste("180510_WT_h2b_gc6s_7dpf_f1")


# load the prefix list and determine the index of the prefix under consideration in the list
list.files(pattern = pattern)
prefix_ID=which(prefix==prefix_list)

# excluded ensembles because only active at the beginning of the recording
exclusions<-vector("list",length(prefix_list))
print(i)

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


list_of_images <- vector("list", length = length(ensel))
for (j in 1:length(ensel)) {
    col = brewer.pal(n = 8, name = "Dark2") [4]
    list_of_images [[j]] <- ggplot(data.frame(centers), aes(x = V1,y= V2)) +geom_point(size = 8, alpha = 0.4, col ="#666666")  + geom_point(data = data.frame(centers [mem==ensel [j] & probs>0.99,]), aes(x = V1,y=V2), size = 10, col =col)+theme_classic()
}

for (i in 1:length(list_of_images)) {
  tiff(paste("../../../../figure_materials/Chapter_1/Figure_3/Assembly_map_Matches_matched_to_firing",folder, i,".tif", sep =""), width = 1000,height= 1000)
  print(list_of_images [[i]])
  dev.off()
}

pdf("../../../../figure_materials/Chapter_1/Figure_5/Compactness_vs_null.pdf")
ggplot(comp, aes(x= comp_value, y = V2, fill = sample)) +  geom_density_ridges( 
     scale = 1.9, alpha =0.6) +  scale_fill_manual(values = c("#D95F02", "#1B9E77"), labels = c("Compactness", "NULL model")) +theme_classic(base_size = 20) + xlab("Compactness") +ylab("Fish Number") + xlim(c(-5e+06,1e+08))  + theme(legend.position = "none")
dev.off()


pdf("../../../../figure_materials/Chapter_1/Figure_5/wcc_vs_null.pdf")
ggplot(wcc, aes(x=wcc_value, y=V2,  fill=sample)) +
  geom_density_ridges( 
    alpha = 0.7, scale = 1.9) + scale_fill_manual(values = c(  "#1B9E77", "#D95F02"), labels = c("NULL model", "Data")) +theme_classic(base_size = 20) + xlab("Within Ensemble Correlation Coeff.") +ylab("Fish Number")  + theme(legend.position = "none")
dev.off()

pdf("../../../../figure_materials/Chapter_1/Figure_5/LI_value_vs_null.pdf")
ggplot(li, aes(x=abs(LI_value), y=V2,  fill=sample)) +
  geom_density_ridges( 
    alpha = 0.7, scale = 1.5) + scale_fill_manual(values = c("#D95F02",  "#1B9E77"), labels = c("NULL model", "Data")) +theme_classic(base_size = 20) + xlab("Within Ensemble Correlation Coeff.") +ylab("Fish Number") +xlim(c(0,1)) + theme(legend.position = "none")
dev.off()

pdf("../../../../figure_materials/Chapter_1/Figure_5/Ensemble_size_number_neurons.pdf")
ggplot(df, aes(x = gs))  + geom_line(aes(group = ID), stat = "density", col = "grey")  +
  geom_line(aes(group=ID), color= "black", stat="density", size=1, alpha=0.2) +
  geom_line( color = "#762A83", size =2, stat="density") +  theme_classic(base_size = 30)  + ylab("Density") + xlab("Ensemble Size (# Neurons)") 
dev.off()





pdf("../../../../figure_materials/Chapter_1/Figure_5/Ensemble_size_number_neurons.pdf")
ggplot(df, aes(x = Assembly))  + geom_line(aes(group = ID), stat = "density", col = "grey")  +
  geom_line(aes(group=ID), color= "black", stat="density", size=1, alpha=0.2) +
  geom_line( color = "#762A83", size =2, stat="density") +  theme_classic(base_size = 30)  + ylab("Density") + xlab("Ensemble Size (# Neurons)") 
dev.off()


pdf("../../../../figure_materials/Chapter_1/Figure_5/lambda_1.pdf")
ggplot(df, aes(x = lambda1))  + geom_line(aes(group = ID), stat = "density", col = "grey")  +
  geom_line(aes(group=ID), color= "black", stat="density", size=1, alpha=0.2) +
  geom_line( color = "#762A83", size =2, stat="density") +  theme_classic(base_size = 30)  + ylab("Density") + xlab("Ensemble Synchrony") 
dev.off()


pdf("../../../../figure_materials/Chapter_1/Figure_5/noise.pdf")
ggplot(df, aes(x = lambda0))  + geom_line(aes(group = ID), stat = "density", col = "grey")  +
  geom_line(aes(group=ID), color= "black", stat="density", size=1, alpha=0.2) +
  geom_line( color = "#762A83", size =2, stat="density") +  theme_classic(base_size = 30)  + ylab("Density") + xlab("Noise") 
dev.off()



pdf("../../../../figure_materials/Chapter_1/Figure_5/Ensemble_Activity.pdf")
ggplot(df, aes(x = ((Assembly_freq)*4.8)*60))  + geom_line(aes(group = ID), stat = "density", col = "grey")  +
  geom_line(aes(group=ID), color= "black", stat="density", size=1, alpha=0.2) +
  geom_line( color = "#762A83", size =2, stat="density") +  theme_classic(base_size = 30)  + ylab("Density") + xlab("Assembly Activity (Events/min)") 
dev.off()







for  (i in 1:length(ensel)){
  tiff(paste("../../../../figure_materials/Chapter_1/Figure_4/Assembly_firing",i,".tif", sep =""), width = 1200,height= 300)
  plot_assembly_w_omega(i)
  dev.off()
}

ggplot(wcc, aes(x=wcc_value, y=V2, fill = sample)) +
  geom_density_ridges_gradient(jittered_points = FALSE, quantile_lines = 
                                 FALSE, quantiles = 2, scale=0.9, color='white') +
  geom_segment(data = iris_lines, aes(x = x0, xend = x0, y = as.numeric(Species),
                                      yend = as.numeric(Species) + .9),
               color = "red") +
  scale_y_discrete(expand = c(0.01, 0)) +
  theme_ridges(grid = FALSE, center = TRUE)
  


shorted_mat <- cbind(df$gs, df$abs_LI, df$norm_eigen_decomp, df$lambda1, df$lambda0,df$Assembly_freq, df$wcc)
colnames(shorted_mat) <- c("size (# neurons)", "Lateral Index", "Ensemble Area", "Synchrony",  "Asynchrony", "Firing Freq", "Within CC")
cor_params <- cor(shorted_mat, method = "spearman")

par(mar=c(7,4,4,2)+0.1) 
my_palette <- colorRampPalette(c("magenta", "black", "green"))(n = 75)
pdf("../../../../figure_materials/Chapter_1/Figure_6/Correlation_between_assembly_parameters.pdf")
heatmap.2(as.matrix(cor_params), col=my_palette ,
          density.info="none", trace="none", dendrogram=c("none"), 
          symm=F,symkey=T,symbreaks=T, scale="none",margins=c(12,12)) 
dev.off()



plot_2dkde <- function(data_frame, x= "gs", y = "lambda1", xlab = "Size (# neurons)", ylab = "Coherence", ylim = c(0,300)) {
  par(mar =c(5,6,2,2))
  heatmap <- kde2d(unlist(data_frame[x]), unlist(data_frame [y]), n= 100)
  
  heatmap$z <- heatmap$z/sum(heatmap$z)
  image.plot(heatmap, col = viridis(80), xlab = xlab, ylab = ylab, cex.axis = 2, cex.lab = 2, ylim = ylim)
}

pdf("../../../../figure_materials/Chapter_1/Figure_6/Abs_Assembly_freq.pdf")
plot_2dkde(df, "abs_LI", y = "Assembly_freq", xlab = "Lateral Index", ylab = "Ensemble firing freq", ylim = c(0.0009,0.015))
points(df$abs_LI, df$Assembly_freq, cex =0.2, pch = 19)
dev.off()

pdf("../../../../figure_materials//Chapter_1/Figure_6/size_sync.pdf")
plot_2dkde(df, "gs", y = "lambda1", xlab = "Ensemble Size (# Neurons)", ylab = "Synchrony", ylim = c(0.05,0.4))
points(x = df$gs, y = df$lambda1, cex = 0.2)
dev.off()


pdf("../../../../figure_materials//Chapter_1/Figure_6/size_Async.pdf")
plot_2dkde(df, "gs", y = "lambda0", xlab = "Ensemble Size (# Neurons)", ylab = "Asynchrony", ylim = c(0.001,0.03))
points(x = df$gs, y = df$lambda0, cex = 0.2)
dev.off()



bootstrap <- function(number_samples, percentage, method = "pearson") {
  cor_dist <- vector("numeric", length = number_samples)
  df_percentile <- (nrow(df)/100)*percentage
  Mutinfor <- vector("numeric", length = number_samples)
  for (i in 1:number_samples) {
    boots <- df [sample(nrow(df), size = df_percentile),]
    cor_dist [i] <- cor(x = boots$gs, y = boots$lambda0, method = method)
    Mutinfor [i] <- MutInf(x = boots$abs_LI, y = boots$Assembly_freq)
  }
  return(data.frame(cor_dist))
}





pdf("../../../../figure_materials/Chapter_1/Figure_6/size_lambda0_bootstrapping2.pdf", height = 2, width = 8 )
ggplot(boots3, mapping = aes(cor_dist)) +geom_density(fill = "#762A83", alpha = 0.4) + theme_classic() + xlim(c(-1,1))
dev.off()



pdf("../../../../figure_materials/Chapter_1/Figure_6/size_lambda1_bootstrapping2.pdf", height = 2, width = 8 )

ggplot(boots2, mapping = aes(cor_dist)) +geom_density(fill = "#762A83", alpha = 0.4) + theme_classic() + xlim(c(-1,1))
dev.off()

mean(unlist(boots2))
sd(unlist(boots2))


boots <- bootstrap(20000,40)
pdf("../../../../figure_materials/Chapter_1/Figure_6/Abs_Assembly_freq_bootstrapping.pdf", height = 2, width = 4)

ggplot(boots, mapping = aes(cor_dist)) +geom_density(fill = "#762A83", alpha = 0.4) + xlim(c(-1,1)) + theme_classic()
dev.off()

mean(unlist(boots))


image.plot(heatmap, col = viridis(80), xlab = 'Ensemble Size (# Neurons)', ylab = "Synchrony", cex.axis = 2, cex.lab = 3)