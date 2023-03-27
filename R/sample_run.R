#############################
# SIMPLE EXAMPLE ON USING MODEL
#############################
source("setup.R")

sc_file <- "sc_pars.csv"
sero_file <- "sero_pars.csv"
ag_file <- "ag_pars.csv" 


## cascade
setup_sim_dir(
  nstrain=9, 
  nhost=5000, 
  strain_id=c(1,2,3,4,5,6,7,8,9), 
  A=c(1,1,1,2,2,2,3,3,3), 
  G=c(1,2,3,1,2,3,1,2,3), 
  O=c(1,1,1,1,1,1,1,1,1),
  V=c(1,1,1,0,0,0,0,0,0), 
  init = c(1,0,0,0,1,0,0,0,1),
  t_max=20*365,
  t_vax=10*365,
  beta=0.09,
  t_g=0.00005,
  t_l=0.00012)

write_sero_pars(alpha=c(0.02,0.024, 0.03),fname=sero_file)
write_ag_pars(kappa=c(0.6),fname=ag_file)
write_sc_pars(kappa=c(0.5,0.48, 0.46),fname=sc_file)

## Run
model_res <- run_simulations(sc_file, sero_file, ag_file, niter=20, label="Sample run")
## Plot relative strain frequencies
p1 <- plot_trajectories(model_res,plot_rel_freq=FALSE)



## no cascade?
setup_sim_dir(
  nstrain=9, 
  nhost=5000, 
  strain_id=c(1,2,3,4,5,6,7,8,9), 
  A=c(1,1,1,2,2,2,3,3,3), 
  G=c(1,2,3,1,2,3,1,2,3), 
  O=c(1,1,1,1,1,1,1,1,1),
  V=c(1,1,1,0,0,0,0,0,0), 
  init = c(1,0,0,0,1,0,0,0,1),
  t_max=20*365,
  t_vax=10*365,
  beta=0.09,
  t_g=0.00005,
  t_l=0.00012)

write_sero_pars(alpha=c(0.02,0.029, 0.035),fname=sero_file)
write_ag_pars(kappa=c(0.6),fname=ag_file)
write_sc_pars(kappa=c(0.5,0.44, 0.47),fname=sc_file)

## Run
model_res2 <- run_simulations(sc_file, sero_file, ag_file, niter=20, label="Sample run")
## Plot relative strain frequencies
p2 <- plot_trajectories(model_res2,plot_rel_freq=FALSE)

