time_to_fire <- function(c) {
  first_fires <- apply(c, 1, function(x) {which(x == !0) [1]})
  first_fires[is.na(first_fires)] <- 0
  return(first_fires)
}



raster_plot2 <- function(stmp1, m = "WT_GRAV",axes =F, xlim = xlim) {
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


plot_raster_sorted <- function(s, fr = c(1,ncol(s))) {
  c <- s [,fr [1]:fr[2]]
  ff <- time_to_fire(c)
  raster_plot(s [order(ff),])
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

plot_assembly_w_omega <- function(ensel_num, xlim = c(1,17460)){
  cells <- s[selection+1, ] [mem == ensel [ensel_num] &  probs > 0.99,]
  dis <- dist(cells)
  cells <- cells [hclust(dis)$order, ]
  layout(matrix(c(1,1,1,1,1,1,2,2), nrow = 4, ncol = 2, byrow = TRUE))
  par( mar = c(0,6,6,2))
  raster_plot(cells, xlim = xlim)
  par( mar = c(2,6,0,2))
  plot(omega_extended [which(ensel == ensel [ensel_num]),], type ='h', xlim = xlim, xaxs = "i", yaxs ="i",  cex.axis =2, ylab = "prob Aseembly Active")
}


plot_assembly_w_omega_ttf <- function(ensel_num, time_of_intrest =c(1, 17460), xlim = c(1,17460)){
  cells <- s[selection+1, ] [mem == ensel [ensel_num] &  probs > 0.99,]
  ff <- time_to_fire(cells [,time_of_intrest])
  dis <- dist(cells)
  cells <- cells [order(ff), ]
  layout(matrix(c(1,1,1,1,1,1,2,2), nrow = 4, ncol = 2, byrow = TRUE))
  par( mar = c(0,6,6,2))
  raster_plot2(cells, xlim = xlim)
  par( mar = c(2,6,0,2))
  plot(omega_extended [which(ensel == ensel [ensel_num]),], type ='h', xlim = xlim, xaxs = "i", yaxs ="i",  cex.axis =2, ylab = "A.A", cex.lab=1.5)
}



plot_assembly_w_omega_ttf_cal <- function(ensel_num, time_of_intrest =c(10, 30), xlim = c(1,17460)){
  cells <- c[selection+1, ] [mem == ensel [ensel_num] &  probs > 0.99,]
  ff <- time_to_fire(cells [,time_of_intrest] >.01 )
  dis <- dist(cells)
  cells <- cells [order(ff), ]
  layout(matrix(c(1,1,1,1,1,1,2,2), nrow = 4, ncol = 2, byrow = TRUE))
  par( mar = c(0,6,6,2))
  image(t(cells) [xlim[1]:xlim[2],])
  par( mar = c(2,6,0,2))
  plot(omega_extended [which(ensel == ensel [ensel_num]),], type ='h', xlim = xlim, xaxs = "i", yaxs ="i",  cex.axis =2, ylab = "A.A", cex.lab=1.5)
}


#png("../../Lab_meeting_plots/Assembly_w_Activity.png", width = 1000, height = 300)
#plot_assembly_w_omega_ttf(10, time_of_intrest = 2500:5000)
#dev.off()
