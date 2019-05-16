library(ggplot2)



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


dat$perc_ensmebled <- dat$ensembled_neurons/dat$cell_num 


plot_line_plot <-function(y="gs", ylab = "Assembly size\n # Neurons", col = c("#028e96","#ea0037")){
  ggplot(dat, aes_string(x = "age", y= y, group = "Rearing_conditions", col = "Rearing_conditions"))  + ylab(ylab)  + stat_summary(fun.y=mean, geom="line", lwd = 1)  + 
  stat_summary(fun.data=mean_se, geom="errorbar", width=0.2,lwd = 1) +
  stat_summary(fun.y=mean, geom="point", size = 5) + xlab("Age (dpf)") + theme_gray(base_size = 30) +
    scale_color_manual(values = col) + theme(legend.title=element_blank())
}

png("../../Lab_meeting_plots/AssemblySize_dev.png", width = 1000, height = 800)
plot_line_plot()
dev.off()

png("../../Lab_meeting_plots/Coherence_dev.png", width = 1000, height = 800)
plot_line_plot(y = "lambda1", ylab = "Coherence")
dev.off()

png("../../Lab_meeting_plots/noise_dev.png", width = 1000, height = 800)
plot_line_plot(y = "lambda0", ylab = "noise")
dev.off()


png("../../Lab_meeting_plots/Area_dev.png", width = 1000, height = 800)
plot_line_plot(y= "areas", ylab = "area of assembly")
dev.off()



png("../../Lab_meeting_plots/Active_Cell_no_dev.png", width = 1000, height = 800)
plot_line_plot(y= "cell_num", ylab = "Number of Active Cells")
dev.off()

png("../../Lab_meeting_plots/Tect_total_area_dev.png", width = 1000, height = 800)
plot_line_plot(y= "tectum_total_area", ylab = "Tectal Total Area")
dev.off()

png("../../Lab_meeting_plots/LI_dev.png", width = 1000, height = 800)
plot_line_plot(y= "abs_LI", ylab = "Lateral Index")
dev.off()

png("../../Lab_meeting_plots/Assembly_area_dev.png", width = 1000, height = 800)
plot_line_plot(y= "norm_eigen", ylab = "Assembly Area\n Normalised to tectal area")
dev.off()

png("../../Lab_meeting_plots/SN_Freq_dev.png", width = 1000, height = 800)
plot_line_plot(y= "sn_freq", ylab = "Frequency of single neurons")
dev.off()


png("../../Lab_meeting_plots/Mean_corr_dev.png", width = 1000, height = 800)
plot_line_plot(y= "mean_corr_total", ylab = "mean_corr_total")
dev.off()

png("../../Lab_meeting_plots/Assembly_no_dev.png", width = 1000, height = 800)
plot_line_plot(y= "numb_ensembles", ylab = "Assembly number")
dev.off()

png("../../Lab_meeting_plots/Assembly_freq_dev.png", width = 1000, height = 800)
plot_line_plot(y= "Assembly_freq", ylab = "Assembley Firing Freq")
dev.off()

png("../../Lab_meeting_plots/Percentage_of_assembly_neurons_dev.png", width = 1000, height = 800)
plot_line_plot(y= "perc_ensmebled", ylab = "Percentage of neurons in an Assembly")
dev.off()
