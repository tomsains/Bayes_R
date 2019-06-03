library(ggplot2)
library(ggsignif)
library(MASS)
library(gplots)


combine_data_tables <- function() {
  file_list <- list.files("mean_across_fish/", pattern = "mean_datatable")
  data_list <- list(length = length(file_list))
  for (i in 1:length(file_list)) {
    data_list [[i]] <- read.table(paste("mean_across_fish/", file_list [[i]], sep = ""))
  }
  return(do.call("rbind", data_list))
  #return(data_list)
}


combine_all_assemblies <- function() {
  file_list <- list.files("mean_across_fish/", pattern = "_all_assemblies")
  data_list <- list(length = length(file_list))
  for (i in 1:length(file_list)) {
    data_list [[i]] <- read.table(paste("mean_across_fish/", file_list [[i]], sep = ""))
  }
  return(do.call("rbind",data_list))
}




dat <- combine_data_tables()

data_summary <- function(data, varname, groupnames){
   require(plyr)
   summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
    sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                       varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}



dat$age <- factor(dat$age, levels = c("3_dpf", "5_dpf", "7_dpf"), ordered =TRUE)

dat$perc_ensmebled <- dat$ensembled_neurons/dat$cell_num 



# NORMAL DEV vs GRAVEL DEV
dev <- dat [dat$Genotype == "WT",]

plot_line_plot <-function(y="gs", ylab = "Assembly size\n # Neurons", col = c("#028e96","#ea0037")){
  ggplot(dev, aes_string(x = "age", y= y, group = "Rearing_conditions", col = "Rearing_conditions"))  + ylab(ylab)  + stat_summary(fun.y=mean, geom="line", lwd = 1)  + 
  stat_summary(fun.data=mean_se, geom="errorbar", width=0.2,lwd = 1) +
  stat_summary(fun.y=mean, geom="point", size = 5) + xlab("Age (dpf)") + theme_gray(base_size = 40) +
    scale_color_manual(values = col) + theme(legend.position="none")
}

png("../../Lab_meeting_plots/AssemblySize_dev.png", width = 800, height = 600)
plot_line_plot()
dev.off()

png("../../Lab_meeting_plots/Coherence_dev.png", width = 800, height = 600)
plot_line_plot(y = "lambda1", ylab = "Coherence")
dev.off()

png("../../Lab_meeting_plots/noise_dev.png", width = 800, height = 600)
plot_line_plot(y = "lambda0", ylab = "noise")
dev.off()


png("../../Lab_meeting_plots/Area_dev.png", width = 800, height = 600)
plot_line_plot(y= "areas", ylab = "area of assembly")
dev.off()



png("../../Lab_meeting_plots/Active_Cell_no_dev.png", width = 800, height = 600)
plot_line_plot(y= "cell_num", ylab = "Number of Active Cells")
dev.off()

png("../../Lab_meeting_plots/Tect_total_area_dev.png", width = 800, height = 600)
plot_line_plot(y= "tectum_total_area", ylab = "Tectal Total Area")
dev.off()

png("../../Lab_meeting_plots/LI_dev.png", width = 800, height = 600)
plot_line_plot(y= "abs_LI", ylab = "Lateral Index")
dev.off()

png("../../Lab_meeting_plots/Assembly_area_dev.png", width = 800, height = 600)
plot_line_plot(y= "norm_eigen", ylab = "Assembly Area\n Normalised to tectal area")
dev.off()

png("../../Lab_meeting_plots/SN_Freq_dev.png", width = 800, height = 600)
plot_line_plot(y= "sn_freq", ylab = "Frequency of single neurons")
dev.off()

plot_line_plot(y= "norm_eigen_decomp", ylab = "Assembly Area")


png("../../Lab_meeting_plots/Mean_corr_dev.png", width = 800, height = 600)
plot_line_plot(y= "mean_corr_total", ylab = "mean_corr_total")
dev.off()

png("../../Lab_meeting_plots/Assembly_no_dev.png", width = 800, height = 600)
plot_line_plot(y= "numb_ensembles", ylab = "Assembly number")
dev.off()

png("../../Lab_meeting_plots/Assembly_freq_dev.png", width = 800, height = 600)
plot_line_plot(y= "Assembly_freq", ylab = "Assembley Firing Freq")
dev.off()

png("../../Lab_meeting_plots/Percentage_of_assembly_neurons_dev.png", width = 800, height = 600)
plot_line_plot(y= "perc_ensmebled", ylab = "Neurons in an Assembly\n (%)")
dev.off()



# Plots for 7dpf data
WT <- dat [dat$age == "7_dpf" & dat$Genotype == "WT", ]
WT$order <- factor(WT$Rearing_conditions, levels = c("NORM", "GRAV"), ordered =TRUE)

plot_WT_7 <- function(y ="gs", ylim = c(20,65), ylab ="Assembly Size\n # Neurons") {
  ggplot(WT, aes_string(x = "order", y = y)) + geom_boxplot(fill = c("#ea0037", "#028e96")) + geom_jitter(width = 0.1) + ylim(ylim) +
    geom_signif(comparisons = list(c("NORM", "GRAV")), textsize = 10, map_signif_level=TRUE, test = "t.test") +
    ylab(ylab) + xlab("Genotype") + theme_gray(base_size = 40) 
  
}


png("../../Lab_meeting_plots/IOP_talk/WT_AssemblySize_dev.png", width = 800, height = 600)
plot_WT_7()
dev.off()

png("../../Lab_meeting_plots/IOP_talk/WT_Coherence_dev.png", width = 800, height = 600)
plot_WT_7(y = lambda1, ylab = "Coherence", ylim = c(0.1,0.24))
dev.off()

png("../../Lab_meeting_plots/IOP_talk/WT_noise_dev.png", width = 800, height = 600)
plot_WT_7(y = "lambda0", ylab = "noise",ylim = c(0.005,0.01))
dev.off()


png("../../Lab_meeting_plots/IOP_talk/Area_dev.png", width = 800, height = 600)
plot_WT_7(y= "areas", ylab = "area of assembly", ylim = c(0,25000))
dev.off()

png("../../Lab_meeting_plots/IOP_talk/WT_Active_Cell_no_dev.png", width = 800, height = 600)
plot_WT_7(y= "cell_num", ylab = "Number of Active Cells", ylim = c(1700,3500))
dev.off()

png("../../Lab_meeting_plots/IOP_talk/WT_LI_dev.png", width = 800, height = 600)
plot_WT_7(y= "abs_LI", ylab = "Lateral Index", ylim = c(0.4,0.9))
dev.off()

png("../../Lab_meeting_plots/IOP_talk/WT_Assembly_area_dev.png", width = 800, height = 600)
plot_WT_7(y= "norm_eigen", ylab = "Assembly Area\n Normalised to tectal area", ylim = c(0.0005,0.003))
dev.off()


png("../../Lab_meeting_plots/IOP_talk/WT_Assembly_area_dev.png", width = 800, height = 600)
plot_WT_7(y= "norm_eigen_decomp", ylab = "Assembly Area\n Normalised to tectal area", ylim = c(0,200))
dev.off()

png("../../Lab_meeting_plots/IOP_talk/SN_Freq_dev.png", width = 800, height = 600)
plot_WT_7(y= "sn_freq", ylab = "Frequency (Events/min)", ylim = c(0,2))
dev.off()


png("../../Lab_meeting_plots/IOP_talk/Mean_corr_dev.png", width = 800, height = 600)
plot_WT_7(y= "mean_corr_total", ylab = "mean_corr_total", ylim = c(0.005,0.01))
dev.off()

png("../../Lab_meeting_plots/IOP_talk/WT_Assembly_no_dev.png", width = 800, height = 600)
plot_WT_7(y= "numb_ensembles", ylab = "Assembly number", ylim = c(20,65))
dev.off()

png("../../Lab_meeting_plots/IOP_talk/WT_Assembly_freq_dev.png", width = 800, height = 600)
plot_WT_7(y= "Assembly_freq", ylab = "Assembley Firing Freq\n (events/min)", ylim = c(0,0.03))
dev.off()

png("../../Lab_meeting_plots/IOP_talk/WT_Percentage_of_assembly_neurons_dev.png", width = 800, height = 600)
plot_WT_7(y= "perc_ensmebled", ylab = "Percentage of neurons in an Assembly", ylim = c(0.4,0.8))
dev.off()


# plot 7dpf across conditions.........
dpf7 <- dat [dat$age == '7_dpf',]
dpf7$Genotype <- factor(dpf7$Genotype, levels = c("WT", "GRIN"), ordered =TRUE)
dpf7$Rearing_conditions <- factor(dpf7$Rearing_conditions, levels = c("NORM", "GRAV"), ordered =TRUE)

plot_dpf_7 <- function(y = "gs", ylab = "Assembly Size \n # neurons", ylim = c(20,65)) {
  ggplot(dpf7, aes_string(x = "Genotype", y = y, fill =  "Rearing_conditions")) + ylab(ylab) +geom_boxplot(outlier.shape = NA) +scale_fill_manual(values=c("#ea0037", "#028e96")) +
    geom_point(position=position_jitterdodge(jitter.width = 0.2)) +
    theme_gray(base_size = 40) + theme(legend.position="none") +ylim(ylim) 
}



png("../../Lab_meeting_plots/IOP_talk/GG_AssemblySize_dev.png", width = 800, height = 600)
plot_dpf_7()
dev.off()

png("../../Lab_meeting_plots/IOP_talk/GG_Coherence_dev.png", width = 800, height = 600)
plot_dpf_7(y = lambda1, ylab = "Coherence", ylim = c(0.1,0.24))
dev.off()

png("../../Lab_meeting_plots/IOP_talk/GG_noise_dev.png", width = 800, height = 600)
plot_dpf_7(y = "lambda0", ylab = "noise",ylim = c(0.005,0.01))
dev.off()


png("../../Lab_meeting_plots/IOP_talk/GG_Area_dev.png", width = 800, height = 600)
plot_dpf_7(y= "areas", ylab = "area of assembly", ylim = c(0,25000))
dev.off()

png("../../Lab_meeting_plots/IOP_talk/GG_Active_Cell_no_dev.png", width = 800, height = 600)
plot_dpf_7(y= "cell_num", ylab = "Number of Active Cells", ylim = c(1500,3500))
dev.off()


png("../../Lab_meeting_plots/IOP_talk/GG_Mean_corr_dev.png", width = 800, height = 600)
plot_dpf_7(y= "mean_corr_total", ylab = "mean_corr_total", ylim = c(0.005,0.01))
dev.off()

png("../../Lab_meeting_plots/IOP_talk/GG_Assembly_no_dev.png", width = 800, height = 600)
plot_dpf_7(y= "numb_ensembles", ylab = "Assembly number", ylim = c(16,68))
dev.off()

png("../../Lab_meeting_plots/IOP_talk/GG_Assembly_freq_dev.png", width = 800, height = 600)
plot_dpf_7(y= "Assembly_freq", ylab = "Assembley Firing Freq\n (events/min)", ylim = c(0,0.03))
dev.off()

png("../../Lab_meeting_plots/IOP_talk/GG_Percentage_of_assembly_neurons_dev.png", width = 800, height = 600)
plot_dpf_7(y= "perc_ensmebled", ylab = "", ylim = c(0.4,0.8))
dev.off()


png("../../Lab_meeting_plots/IOP_talk/GG_LI_dev.png", width = 800, height = 600)
plot_dpf_7(y= "abs_LI", ylab = "Lateral Index", ylim = c(0.5,0.8))
dev.off()

png("../../Lab_meeting_plots/IOP_talk/GG_Assembly_area_dev.png", width = 800, height = 600)
plot_dpf_7(y= "norm_eigen", ylab = "Assembly Area\n Normalised to tectal area", ylim = c(0.0005,0.003))
dev.off()


png("../../Lab_meeting_plots/IOP_talk/GG_Assembly_area_dev.png", width = 800, height = 600)
plot_dpf_7(y= "norm_eigen_decomp", ylab = "Assembly Area\n Normalised to tectal area", ylim = c(0,200))
dev.off()

png("../../Lab_meeting_plots/IOP_talk/GG_SN_Freq_dev.png", width = 800, height = 600)
plot_dpf_7(y= "sn_freq", ylab = "Frequency (Events/min)", ylim = c(0.8,2))
dev.off()







# Chapter_1 figures


WT_N_7 <- combine_all_assemblies()
WT_N_7 <- WT_N_7 [WT_N_7$age == "7_dpf" & WT_N_7$Genotype == "WT" & WT_N_7$Rearing_conditions == "NORM", ]

scatter_with_histograms <- function(df, x , y, xlab, ylab) {
  p <- ggplot(df, aes_string(x, y)) + geom_point() + geom_density_2d()  + theme_classic() + scale_fill_viridis(begin = 0, end= 1, option ="A") + 
    stat_density_2d(geom = "raster", aes(fill = stat(density)), n =50, contour = FALSE) + theme(legend.position="none") +
    ylab(ylab) + xlab(xlab)
  
  #+ 
    #scale_y_continuous(expand = c(0,0)) +  scale_x_continuous(expand = c(0,0)) 
  ggExtra::ggMarginal(p, type = "histogram", fill = "#595861", size = 2.5)
}

plot_density_plots <- function(x = "gs", xlab = "Size Assembly (# Neurons)" ) {
  ggplot(WT_N_7, aes_string(x=x)) + 
    geom_density(color = "#474747", fill =  "#028e96", alpha =.7) + theme_classic(base_size = 40) + xlab(xlab) +scale_y_continuous(expand = c(0,0)) +  scale_x_continuous(expand = c(0,0.0))
}


png("../../Lab_meeting_plots/Norm_7dpf_WT_plots/ch1_abs_LI_density.png", width = 800, height = 600)
plot_density_plots(x = "abs_LI", "Laterality Index")
dev.off()

png("../../Lab_meeting_plots/Norm_7dpf_WT_plots/ch1_LI_density.png", width = 800, height = 600)
plot_density_plots(x = "LI", "Laterality Index")
dev.off()

png("../../Lab_meeting_plots/Norm_7dpf_WT_plots/ch1_abs_Assembly_size.png", width = 800, height = 600)
plot_density_plots(x = "gs", "Size of Assembly (# Neurons)")
dev.off()

png("../../Lab_meeting_plots/Norm_7dpf_WT_plots/ch1_Areas_ch.png", width = 800, height = 600)
plot_density_plots(x = "areas", "Area (Convex hull)")
dev.off()


png("../../Lab_meeting_plots/Norm_7dpf_WT_plots/ch1_TTA.png", width = 800, height = 600)
plot_density_plots(x = "tectum_total_area", "Tectum Total Area (pixels)")
dev.off()


png("../../Lab_meeting_plots/Norm_7dpf_WT_plots/ch1_lambda1.png", width = 800, height = 600)
plot_density_plots(x = "lambda1", "Coherence")
dev.off()


png("../../Lab_meeting_plots/Norm_7dpf_WT_plots/ch1_lambda0.png", width = 800, height = 600)
plot_density_plots(x = "lambda0", "Noise")
dev.off()


png("../../Lab_meeting_plots/Norm_7dpf_WT_plots/ch1_Assembly_ff.png", width = 800, height = 600)
plot_density_plots(x = "Assembly_freq", "Assembly Firing Freq\n (Events/min)")
dev.off()


png("../../Lab_meeting_plots/Norm_7dpf_WT_plots/ch1_norm_eigen.png", width = 800, height = 600)
plot_density_plots(x = "norm_eigen", "Assembly area\n (Normalised)")
dev.off()

# Correlation between parameters
shorted_mat <- WT_N_7 [c("gs", "abs_LI", "lambda1", "lambda0","Assembly_freq", "wcc", "norm_eigen")]
colnames(shorted_mat) <- c("size (# neurons)", "Lateral Index",  "Coherence",  "noise", "Firing Freq", "Within CC", "Area")
cor_params <- cor(shorted_mat, method = "spearman")

par(mar=c(7,4,4,2)+0.1) 
my_palette <- colorRampPalette(c("magenta", "black", "green"))(n = 75)
png(filename="../../Lab_meeting_plots/Norm_7dpf_WT_plots/ch1_heatmap.png", width=800, height=750)
heatmap.2(as.matrix(cor_params), col=my_palette ,
          density.info="none", trace="none", dendrogram=c("none"), 
          symm=F,symkey=T,symbreaks=T, scale="none",margins=c(12,12)) 
dev.off()

############################# heatmaps
plot_2dkde <- function(data_frame, x= "gs", y = "lambda1", xlab = "Size (# neurons)", ylab = "Coherence", ylim = c(0,300)) {
  par(mar =c(5,6,2,2))
  heatmap <- kde2d(unlist(data_frame[x]), unlist(data_frame [y]), n= 100)
  library(RColorBrewer)
  rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
  r <- rf(32)
  heatmap$z <- heatmap$z/sum(heatmap$z)
  image(heatmap, col = r, xlab = xlab, ylab = ylab, cex.lab =3, cex.axis = 2, ylim = ylim)
}



png("../../Lab_meeting_plots/Norm_7dpf_WT_plots/ch1_heatmap_gs_lambda1.png", width = 800, height = 600)
plot_2dkde(data_frame = WT_N_7, x = "gs", y = "lambda1", xlab = "Size (# Neurons)", ylim = c(0.06,.53))
dev.off()


png("../../Lab_meeting_plots/Norm_7dpf_WT_plots/ch1_heatmap_LI_Area.png", width = 800, height = 600)
plot_2dkde(data_frame = WT_N_7, x = "LI", y = "norm_eigen", xlab = "Lateral Index", "Assembly Area", ylim = c(0,300))
dev.off()

png("../../Lab_meeting_plots/Norm_7dpf_WT_plots/ch1_heatmap_gs_AFF.png", width = 800, height = 600)
plot_2dkde(data_frame = WT_N_7, x = "gs", y = "Assembly_freq", xlab = "Size (# Neurons)", ylab = "Assembly Firing Freq", ylim = c(0.001,0.03))
dev.off()

png("../../Lab_meeting_plots/Norm_7dpf_WT_plots/ch1_heatmap_LI_AFF.png", width = 800, height = 600)
plot_2dkde(data_frame = WT_N_7, x = "abs_LI", y = "Assembly_freq", xlab = "lateral index", ylab = "Assembly Firing Freq", ylim = range(WT_N_7$Assembly_freq))
dev.off()

png("../../Lab_meeting_plots/Norm_7dpf_WT_plots/ch1_heatmap_Area_AFF.png", width = 800, height = 600)
plot_2dkde(data_frame = WT_N_7, x = "norm_eigen", y = "Assembly_freq", xlab = "Assembly Area (normalised)", ylab = "Assembly Firing Freq", ylim = range(WT_N_7$Assembly_freq))
dev.off()
