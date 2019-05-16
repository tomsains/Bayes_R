# parameter_maps_across_the_tectum


plot_all_centers <- function(folder, pattern = "_sa_aligned_all_cells_centers_cut_ordered.dat") {
  files <- list.files(path = folder, pattern = pattern)
  list_centers <- list(length = length(files))
  plot(as.matrix(read.table(file = paste(folder,files[[1]], sep =""))) [,1:2], cex = 0.5, pch = 19, col = rgb(0,0,0,0.5))
  list_centers [[1]] <-  as.matrix(read.table(file = paste(folder,files[[1]], sep ="")))
  for (i in 2:length(files)) {
    points(as.matrix(read.table(file = paste(folder, "/", files[[i]],sep =""))) [,1:2], cex = 0.5, pch = 19, col = rgb(0,0,0,0.5))
    list_centers [[i]] <-  as.matrix(read.table(file = paste(folder, "/", files[[i]],sep =""))) 
  }
  return(do.call(rbind,list_centers))
}


# make maps
pointsize = 5
basesize = 30


m1 <- ggplot(data = df2 , aes(x = centers.1, y = centers.2))  +
  geom_point(aes(colour =  lambda1, size= pointsize)) + scale_colour_gradient2(low = "black", mid = "red",high = "yellow", midpoint = 0.3) +
  theme_classic(base_size = basesize) 

m2 <- ggplot(data = df2 , aes(x = centers.1, y = centers.2))  +
  geom_point(aes(colour =  lambda0, size= pointsize)) + scale_colour_gradient2(low = "black", mid = "red",high = "yellow", midpoint = 0.018) +
  theme_classic(base_size = basesize) 

m3 <- ggplot(data = df2 , aes(x = centers.1, y = centers.2))  +
  geom_point(aes(colour =  gs, size= pointsize)) + scale_colour_gradient2(low = "black", mid = "red",high = "yellow", midpoint = 60) +
  theme_classic(base_size = basesize) 

m4 <-  ggplot(data = df2 , aes(x = centers.1, y = centers.2))  +
  geom_point(aes(colour =  LI, size= pointsize)) + scale_colour_gradient2(low = "blue", mid = "black",high = "red", midpoint = 0) +
  theme_classic(base_size = basesize) 

m5 <-  ggplot(data = df2 , aes(x = centers.1, y = centers.2))  +
  geom_point(aes(colour =  Assembly_freq, size= pointsize)) + scale_colour_gradient2(low = "black", mid = "red",high = "yellow", midpoint = 0.015) +
  theme_classic(base_size = basesize) 

m6 <-  ggplot(data = df2 , aes(x = centers.1, y = centers.2))  +
  geom_point(aes(colour =  wcc, size= pointsize)) + scale_colour_gradient2(low = "black", mid = "red",high = "yellow", midpoint = 0.2) +
  theme_classic(base_size = basesize) 

m7 <-  ggplot(data = df2 , aes(x = centers.1, y = centers.2))  +
  geom_point(aes(colour =  areas, size= pointsize)) + scale_colour_gradient2(low = "black", mid = "red",high = "yellow", midpoint = 10000) +
  theme_classic(base_size = basesize) 


m8<-  ggplot(data = df2 , aes(x = centers.1, y = centers.2))  +
  geom_point(aes(colour =  norm_eigen, size= pointsize)) + scale_colour_gradient2(low = "black", mid = "red",high = "yellow", midpoint = 0.004) +
  theme_classic(base_size = basesize) 

m9 <-  ggplot(data = df2, aes(x = centers.1, y = centers.2))  +
  geom_point(aes(colour =  norm_eigen_decomp, size= pointsize)) + scale_colour_gradient2(low = "black", mid = "red",high = "yellow", midpoint = 0.004) +
  theme_classic(base_size = basesize) 

m10 <-  ggplot(data = df2, aes(x = centers.1, y = centers.2))  +
     geom_point(aes(colour =  pmu, size= pointsize)) + scale_colour_gradient2(low = "black", mid = "red",high = "yellow", midpoint = 0.1) +
     theme_classic(base_size = basesize) 

pdf("Plots/WT_norm_tectal_heatmaps_of_params_1.pdf", width = 50, height = 25)
multiplot(m1,m2,m3,m4,m5,m6, cols = 3)
dev.off()

pdf("Plots/WT_norm_tectal_heatmaps_of_params_2.pdf", width = 50, height = 25)
multiplot(m7,m8,m9,m10, cols = 2)
dev.off()

