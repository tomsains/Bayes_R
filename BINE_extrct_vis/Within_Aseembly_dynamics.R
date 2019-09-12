time_to_fire <- function(c) {
  first_fires <- apply(c, 1, function(x) {which(x == !0) [1]})
  first_fires[is.na(first_fires)] <- 0
  return(first_fires)
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
  raster_plot2(do.call("rbind", cell_list))
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


new_coherence <- function(omega_extended) {
  omega_extended_smooth <- t(apply(omega_extended, 1, rollmean, k= 5, fill = c(0), align ="center"))
  start
}



reverberating_activity <- function(omega_extended) {
  response_duration <- function(t) {
    runs <- rle(t)
    return(runs$lengths[runs$values==1])
  }
  omega_extended <- t(apply(omega_extended, 1, rollmean, k= 5, fill = c(0), align ="center")) 
  durations <- apply(X = omega_extended, 1, response_duration)
  thresh_durations <- lapply(durations, function(x) {x [x > 3]})
  mean_Ass_duration <- lapply(thresh_durations, mean, na.rm =T)
  max_Ass_duration <- lapply(thresh_durations, max, na.rm = T)
  bursting_events <- lapply(thresh_durations, length)
  mean_Ass_duration [is.nan(unlist(mean_Ass_duration))] <- 0
  max_Ass_duration [max_Ass_duration == -Inf] <- 0
  return(data.frame(mean_Ass_duration = unlist(mean_Ass_duration), max_Ass_dur =unlist(max_Ass_duration), bursting_events = unlist(bursting_events)))
}



new_coherence <- function(ensel_num, thresh =0.2) {
  cells <- s[selection+1, ] [mem == ensel [ensel_num] &  probs > 0.99,]
  omega_extended_smooth <- rollmean(omega_extended [which(ensel == ensel [ensel_num]),] > 0.8, k= 5, fill = c(0), align ="center")
 
  w_start <- which(diff(omega_extended_smooth > thresh)  == 1)
  w_end <- which(diff(omega_extended_smooth > thresh) == -1)
  awa <- matrix(NA,nrow = nrow(cells), ncol = length(w_start))
  if (length(w_start) == 0) {
    return(0)
  }
    else{
    for (i in 1:length(w_start)) {awa [,i] <- rowSums(cells [,w_start[i]:w_end[i]])}
    return(mean(colSums(awa >0)/nrow(awa)))
  }
}



plot_assembly_w_omega_and_map <- function(ensel_num, xlim = c(1,17460)){
  cells <- s[selection+1, ] [mem == ensel [ensel_num] &  probs > 0.99,]
  cent <- centers [mem == ensel [ensel_num] &  probs > 0.99,]
  cor_with_cluster <- cor(omega_extended [which(ensel == ensel [ensel_num]),] > 0.8, t(cells))
  cells <- cells [order(cor_with_cluster,decreasing = TRUE), ]
  layout(matrix(c(1,2,2,1,2,2,0,3,3), nrow = 3, ncol = 3, byrow = TRUE))
  par(mar = c(2,2,2,2))
  plot(centers [,1:2], pch = 19 , cex = 0.2)
  points(cent [,1:2], col ='red', pch = 19)
  par( mar = c(0,2,2,2))
  raster_plot2(cells, xlim = xlim)
  par( mar = c(2,2,0,2))
  plot(omega_extended [which(ensel == ensel [ensel_num]),], type ='h', xlim = xlim, xaxs = "i", yaxs ="i",  cex.axis =2, ylab = "prob Aseembly Active")
}

for (i in 1:length(ensel)) {
  png(paste(folder, "/Assembly_", ensel [i], ".png", sep = ''), width = 1000, height = 300)
  plot_assembly_w_omega_and_map(i)
  dev.off()
}
