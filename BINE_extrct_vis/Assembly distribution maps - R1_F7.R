# plot map of assembly centers - parametes


plot_all_centers <- function(folder){
  files <- list.files(path = folder, pattern = "registered")
  df <- as.data.frame(read.table(paste(folder, files [[1]], sep  = ""))) [,1:2]
  for (i in 1:length(files)) {
    df <- rbind(df, as.data.frame(read.table(paste(folder, files [[i]], sep  = ""))[,1:2]))
    
  }
  plot(df [,1:2])
  return(df)
}



all_centers  <- plot_all_centers(folder = "../data_for_Gibbs/7_dpf/")




df_aa <- read.table(paste("mean_across_fish/", genotype,"_", age,"_", Rearing_conditions, "_all_assemblies_datatable_registered.dat", sep =""))

high_LI <- df_aa [ df_aa$abs_LI > 0.5, ]


tiff("../../../../figure_materials/Chapter_1/Figure_7/assembly_density.tif", height = 1000, width = 1003, res = 100)
image(kde2d(high_LI$centers.1, high_LI$centers.2, n = 2000, lims = c(-180, 120, -120, 180),h = 25), col = viridis(200))
points(all_centers [,1:2], col = rgb(red = 1, green = 1, blue = 1, alpha = 0.1), pch =19, cex= 0.3)
dev.off()


tiff("../../../../figure_materials/Chapter_1/Figure_7/Assemblys_over_tectum.tif", height = 1000, width = 1003, res = 100)
ggplot(all_centers, aes(V1, V2)) + geom_point(alpha = 0.2, cex = 0.5) + geom_point(data = high_LI, aes(x = centers.1, y = centers.2), cex = 4) +theme_classic()
dev.off()


tiff("../../../../figure_materials/Chapter_1/Figure_7/synchrony_assembly_map.tif", height = 1000, width = 1153, res = 100)
ggplot(all_centers, aes(V1, V2)) + geom_point(alpha = 0.2, cex = 0.5) + geom_point(data = high_LI, aes(x = centers.1, y = centers.2, col = lambda1), alpha = 0.7, size = 6) +theme_classic()  + 
  scale_colour_viridis() 
dev.off()

tiff("../../../../figure_materials/Chapter_1/Figure_7/Asynchrony_assembly_map.tif", height = 1000, width = 1153, res = 100)
ggplot(all_centers, aes(V1, V2)) + geom_point(alpha = 0.2, cex = 0.5) + geom_point(data = high_LI, aes(x = centers.1, y = centers.2, col = lambda0), alpha = 0.7, size = 6) +theme_classic()  + 
  scale_colour_viridis()
dev.off()

tiff("../../../../figure_materials/Chapter_1/Figure_7/LI_assembly_map.tif", height = 1000, width = 1153, res = 100)
ggplot(all_centers, aes(V1, V2)) + geom_point(alpha = 0.2, cex = 0.5) + geom_point(data = high_LI, aes(x = centers.1, y = centers.2, col = LI), alpha = 0.7, size = 6) +theme_classic()  + 
  scale_colour_viridis()
dev.off()

tiff("../../../../figure_materials/Chapter_1/Figure_7/Assembly_freq_assembly_map.tif", height = 1000, width = 1153, res = 100)
ggplot(all_centers, aes(V1, V2)) + geom_point(alpha = 0.2, cex = 0.5) + geom_point(data = high_LI, aes(x = centers.1, y = centers.2, col = Assembly_freq), alpha = 0.7, size = 6) +theme_classic()  + 
  scale_colour_viridis()
dev.off()

tiff("../../../../figure_materials/Chapter_1/Figure_7/Size_assembly_map.tif", height = 1000, width = 1153, res = 100)
ggplot(all_centers, aes(V1, V2)) + geom_point(alpha = 0.2, cex = 0.5) + 
  geom_point(data = high_LI, aes(x = centers.1, y = centers.2, col = gs), alpha = 0.7, size = 6) +theme_classic()  + 
  scale_colour_viridis()
dev.off()

tiff("../../../../figure_materials/Chapter_1/Figure_7/wcc_assembly_map.tif", height = 1000, width = 1153, res = 100)
ggplot(all_centers, aes(V1, V2)) + geom_point(alpha = 0.2, cex = 0.5) + geom_point(data = high_LI, aes(x = centers.1, y = centers.2, col = wcc), alpha = 0.7, size = 6) +theme_classic()  + 
  scale_colour_viridis()
dev.off()


tiff("../../../../figure_materials/Chapter_1/Figure_7/area_assembly_map.tif", height = 1000, width = 1153, res = 100)
ggplot(all_centers, aes(V1, V2)) + geom_point(alpha = 0.2, cex = 0.5) + geom_point(data = high_LI, aes(x = centers.1, y = centers.2, col = norm_eigen_decomp), alpha = 0.7, size = 6) +theme_classic()  + 
  scale_colour_viridis()
dev.off()