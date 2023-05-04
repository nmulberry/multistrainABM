source("setup.R") 


###############################
# -----OPTIONS--------
###############################

# Path to model binary
multiabm <- "~/multistrainABM/target/release/multiabm_samplestrains"

# Whether or not to load generated data (if exists), or run all simulations 
# Note: some simulations may take a while to run
load_data <- TRUE

# If running, save or not
save_data <- TRUE

# If saving, specify directory
data_dir <- "./generated-data"
dir.create(data_dir, showWarnings = FALSE)

# Output directory (figures)
out_dir <- "./figures"
dir.create(out_dir, showWarnings = FALSE)

#############################
#-- MAIN SIMULATIONS-------
#############################

## Figure 3
source("withinhost-illustration.R") # pretty quick to run, no need to save data

# Rest of main figs
NITER <- no_cores*3
NHOST <- 7000
# these ones take a while, may not want to run locally
#source("figure-4.R") 

#source("figure-5.R") 

source("figure-6.R")
