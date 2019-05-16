# make raster plots for each condition


s <- read.table("../data_for_Gibbs/7_dpf/180220_WT_h2b_gc6s_7dpf_f1_sa_aligned_all_cells_spikes.dat")
png("Plots/WT_NORM_RASTER.png", width = 1000, height = 600)
raster_plot(s)
dev.off()


s <- read.table("../data_for_Gibbs/Grav_7_dpf/180530_WT_grav_h2b_gc6s_7dpf_f1_sa__00001_scaled_aligned_all_cells_spikes.dat")
png("Plots/WT_GRAV_RASTER.png", width =1000,height =600)
raster_plot(s, "WT_GRAV")
dev.off()


s <- read.table("../data_for_Gibbs/Grin2a/180705_Grin2aKO_h2b_gc6s_7dpf_f2_sa__00001_scaled_aligned_all_cells_spikes.dat")
png("Plots/GRIN_NORM_RASTER.png",width = 1000, height =600)
raster_plot(s, "GRIN_NORM")
dev.off()


s <- read.table("../data_for_Gibbs/Grin_grav_7dpf/180927_grin_grav_h2b_gc6s_7dpf_f1_sa__00001_scaled_aligned_all_cells_spikes.dat")
png("Plots/GRIN_GRAV.png", width = 1000, height = 600)
raster_plot(s, "GRIN_GRAV")
dev.off()

s <- read.table("../data_for_Gibbs/Grav_3_dpf/190304_grav_3dpf_h2b_GCAMP6_sa_F1_00001_scaled_aligned_all_cells_spikes.dat")
raster_plot(s, "GRAV_3")

s <- read.table("../data_for_Gibbs/Grav_3_dpf/190304_grav_3dpf_h2b_GCAMP6_sa_F3_00002_scaled_aligned_all_cells_spikes.dat")
raster_plot(s, "GRAV_3")