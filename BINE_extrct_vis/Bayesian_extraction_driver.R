# driver for extracting Bayesian data

LI_filter = FALSE

# WT data
#########################################################################################################################################################
#########################################################################################################################################################

# 7_dpf 

#########################################################################################################################################################

setwd("/media/thomas_sainsbury/Samsung_T5/Bayesian_network_inference/BPTS_Tom")

data_location = "../data_for_Gibbs/7_dpf/"
centers_suffix = "_sa_aligned_all_cells_centers_cut_ordered.dat"
spikes_suffix = "_sa_aligned_all_cells_spikes.dat"
midline <- data.frame(cbind(intercept = c(295, 312, 267, 345, 386, 327, 243, 362), slope = c(-1.23, -1.23, -1, -1.57, -1.828, -1.38, -1.012, -1.62)))
#midline <- data.frame(cbind(intercept = c(10, 10, 10, 10, 10, 10, 10, 10), slope = c(-0.97, -0.97, -0.97, -0.97, -0.97, -0.97, -0.97, -0.97)))
pattern <- "*WT_h2b_gc6s_7dpf_*"
reg = FALSE

genotype = "WT"
age = "7_dpf"
Rearing_conditions = "NORM"
corr_data = paste(genotype,"_", age,"_", Rearing_conditions, sep = "")


#snp(folder = data_location, save_name = corr_data)
#source("Generalised_membership_v4.R")
#source("Generalise_make_figures.R")
#Check_midlines(midline = midline, folder = data_location)
source("../Bayes_R/BINE_extrct_vis/ensemble_features_bySample.R")


# 5_dpf 

#########################################################################################################################################################



data_location = "../data_for_Gibbs/WT_5_dpf/"
centers_suffix = "_sa_aligned_all_cells_centers_cut_ordered.dat"
spikes_suffix = "_sa_aligned_all_cells_spikes.dat"


midline <- data.frame(cbind(intercept = c(400, 365, 396, 308, 393, 304, 200, 394, 398), slope = c(-2, -1.74, -1.87, -1.256, -1.788, -1.372, -0.64, -1.88, -1.91)))

pattern <- "*_WT_h2b_gc6s_5dpf_*"
reg = FALSE

genotype = "WT"
age = "5_dpf"
Rearing_conditions = "NORM"
corr_data = paste(genotype,"_", age,"_", Rearing_conditions, sep = "")
#Check_midlines(midline = midline, folder = data_location)

#snp(folder = data_location, save_name = corr_data)
#source("Generalised_membership_v4.R")
#source("Generalise_make_figures.R")
source("../Bayes_R/BINE_extrct_vis/ensemble_features_bySample.R")




# 3_dpf 

#########################################################################################################################################################



data_location = "../data_for_Gibbs/WT_3_dpf/"
centers_suffix = "_sa_aligned_all_cells_centers_cut_ordered.dat"
spikes_suffix = "_sa_aligned_all_cells_spikes.dat"


midline <- data.frame(cbind(intercept = c(400, 384.8, 382, 262, 277, 242), slope = c(-2,-1.85,-1.73, -1.12, -1.3, -0.98)))

pattern <- "*_WT_h2b_gc6s_3dpf_*"
reg = FALSE

genotype = "WT"
age = "3_dpf"
Rearing_conditions = "NORM"
corr_data = paste(genotype,"_", age,"_", Rearing_conditions, sep = "")
Check_midlines(midline = midline, folder = data_location)

#snp(folder = data_location, save_name = corr_data)


#source("Generalised_membership_v4.R")
#source("Generalise_make_figures.R")
source("../Bayes_R/BINE_extrct_vis/ensemble_features_bySample.R")




# WT GRAV 
########################################################################################################################################################

########################################################################################################################################################

#WT_7_dpf
########################################################################################################################################################
data_location = "../data_for_Gibbs/Grav_7_dpf/"
centers_suffix = "_scaled_aligned_all_cells_centers_cut_ordered.dat"
spikes_suffix = "_scaled_aligned_all_cells_spikes.dat"
cal_suffix ="_scaled_aligned_all_cells_cal.dat"


midline <- data.frame(cbind(intercept = c(321, 317, 280, 235,484,382, 477), slope = c(-1.468, -1.33, -1.15, -0.78, -2.57, -1.77, -2.36)))

pattern <- "*_WT_grav_h2b_gc6s_7dpf_*"
reg = FALSE

genotype = "WT"
age = "7_dpf"
Rearing_conditions = "GRAV"
corr_data = paste(genotype,"_", age,"_", Rearing_conditions, sep = "")

#Check_midlines(midline = midline, folder = data_location)
#snp(folder = data_location, save_name = corr_data)

#source("Generalised_membership_v4.R")
#source("Generalise_make_figures.R")

source("../Bayes_R/BINE_extrct_vis/ensemble_features_bySample.R")


#WT_5_dpf
########################################################################################################################################################
data_location = "../data_for_Gibbs/Grav_5_dpf/"
centers_suffix = "_scaled_aligned_all_cells_centers_cut_ordered.dat"
spikes_suffix = "_scaled_aligned_all_cells_spikes.dat"

midline <- data.frame(cbind(intercept = c(185, 160, 175, 360, 205), slope = c(-0.7, -0.4,-0.5,-1.8, -0.7)))

pattern <- "*_grav_5dpf_h2b_GCAMP6_sa_*"
reg = FALSE

genotype = "WT"
age = "5_dpf"
Rearing_conditions = "GRAV"
corr_data = paste(genotype,"_", age,"_", Rearing_conditions, sep = "")
Check_midlines(midline = midline, folder = data_location)

#snp(folder = data_location, save_name = corr_data)

#source("Generalised_membership_v4.R")
#source("Generalise_make_figures.R")
source("../Bayes_R/BINE_extrct_vis/ensemble_features_bySample.R")

#WT_3_dpf
########################################################################################################################################################
data_location = "../data_for_Gibbs/Grav_3_dpf/"
centers_suffix = "_scaled_aligned_all_cells_centers_cut_ordered.dat"
spikes_suffix = "_scaled_aligned_all_cells_spikes.dat"


midline <- data.frame(cbind(intercept = c(230,165), slope = c(-0.8, -0.3)))

pattern <- "*_grav_3dpf_h2b_GCAMP6_sa_*"
reg = FALSE

genotype = "WT"
age = "3_dpf"
Rearing_conditions = "GRAV"
corr_data = paste(genotype,"_", age,"_", Rearing_conditions, sep = "")
Check_midlines(midline = midline, folder = data_location)

#snp(folder = data_location, save_name = corr_data)

#source("Generalised_membership_v4.R")
#source("Generalise_make_figures.R")
source("../Bayes_R/BINE_extrct_vis/ensemble_features_bySample.R")












########################################################################################################################################################
########################################################################################################################################################

# GRIN2a

#Grin 7_dpf_Grav

data_location = "../data_for_Gibbs/Grin_grav_7dpf/"
centers_suffix = "_scaled_aligned_all_cells_centers_cut_ordered.dat"
spikes_suffix = "_scaled_aligned_all_cells_spikes.dat"


midline <- data.frame(cbind(intercept = c(335, 331, 449, 273), slope = c(-1.62, -1.41, -2.1, -1.12)))

pattern <- "*_grin_grav_h2b_gc6s_*"
reg = FALSE

genotype = "GRIN"
age = "7_dpf"
Rearing_conditions = "GRAV"
corr_data = paste(genotype,"_", age,"_", Rearing_conditions, sep = "")
Check_midlines(midline = midline, folder = data_location)

#snp(folder = data_location, save_name = corr_data)

#source("Generalised_membership_v4.R")
#source("Generalise_make_figures.R")
source("../Bayes_R/BINE_extrct_vis/ensemble_features_bySample.R")


#Grin 7_dpf_NORM

data_location = "../data_for_Gibbs/Grin2a/"
centers_suffix = "_scaled_aligned_all_cells_centers_cut_ordered.dat"
spikes_suffix = "_scaled_aligned_all_cells_spikes.dat"



midline <- data.frame(cbind(intercept = c(364, 298, 362, 357), slope = c(-1.61, -1.31, -1.7, -1.61)))


pattern <- "*Grin2a*"
reg = FALSE

genotype = "GRIN"
age = "7_dpf"
Rearing_conditions = "NORM"
corr_data = paste(genotype,"_", age,"_", Rearing_conditions, sep = "")
Check_midlines(midline = midline, folder = data_location)

#snp(folder = data_location, save_name = corr_data)

#source("Generalised_membership_v4.R")
#source("Generalise_make_figures.R")
source("../Bayes_R/BINE_extrct_vis/ensemble_features_bySample.R")

