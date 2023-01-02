source("setup.R")
out_dir <- "figures-main/" # output directory for figures
dir.create(out_dir, showWarnings=FALSE)

# set some default parameters
NHOST <- 10000 #per strain default ---> bump to 10k ?
ALPHA = 0.03
P_TAU = 0.4
COST = 0.15
TAU=0.35
T_MAX = 15*365
T_VAX=10*365
NITER <- no_cores*3
BETA <- 0.1
MU <- 0.0
dir.create("tmp", showWarnings=FALSE) #for storing temp files
#--------------------------#
#--------------------------#
# Main Experiments & Results
#--------------------------#
#--------------------------#

#############################
# RUN NEUTRALITY ANALYSIS
############################
setwd("tmp")
source("../neutrality.R")
setwd("../")
ggsave(paste0(out_dir, "neutrality_analysis.pdf"), gg_neutral, height=9, width=11)

#############################
# RUN WITHINHOST ILLUSTRATIONS
##############################
setwd("tmp")
source("../withinhost.R")
setwd("../")
ggsave(paste0(out_dir, "sens-res-main-fig.pdf"), gg_res, height=10, width=11)

#########################
# RUN COEXISTENCE
#########################
setwd("tmp")
source("../coexistence.R")
setwd("../")
ggsave(paste0(out_dir, "sens_res_coex.pdf"),gg_res, height=8, width=12)
ggsave(paste0(out_dir, "sens_res_coex_alpha.pdf"),gg_res_alpha, height=8, width=12)
#########################
# RUN VACCINATION ANALYSIS
#########################


setwd("tmp")
source("../2strain-vax.R")
setwd("../")
ggsave(paste0(out_dir, "2strain_vax.pdf"), gg_2strainvax, width=8, height=8)
ggsave(paste0(out_dir, "2strain_vax_traj.pdf"), gg_traj, width=8, height=8)
ggsave(paste0(out_dir, "2strain_vax_freqs.pdf"), gg_freqs, width=8, height=10)
###
setwd("tmp")
source("../3strain-vax.R")
setwd("../")
ggsave(paste0(out_dir, "3strain_vax.pdf"), gg_3strainvax, width=8, height=8)
ggsave(paste0(out_dir, "3strain_vax_freqs.pdf"), gg_freqs, width=8, height=10)

### Some trajectories too
setwd("tmp")
source("../vims_trajectories.R")
setwd("../")
ggsave(paste0(out_dir, "vims-toy.pdf"), gg_traj, height=10, width=10)
unlink("tmp", recursive=TRUE)



