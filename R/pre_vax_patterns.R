source("setup.R")

sc_file <- "sc_pars.csv"
sero_file <- "sero_pars.csv"
ag_file <- "ag_pars.csv" 


## understand 3 sero, 3 gpsc dynamics

setup_sim_dir(
  nstrain=9, 
  nhost=5000, 
  strain_id=c(1,2,3,4,5,6,7,8,9), 
  A=c(1,1,1,2,2,2,3,3,3), 
  G=c(1,2,3,1,2,3,1,2,3), 
  O=c(1,1,1,1,1,1,1,1,1),
  V=c(0,0,0,0,0,0,0,0,0), 
  init = c(1,0,0,0,1,0,0,0,1),
  t_max=30*365,
  t_vax=30*365,
  beta=0.09,
  t_g=0.00005,
  t_l=0.00012)

write_sero_pars(alpha=c(0.02,0.025, 0.03),fname=sero_file)
write_ag_pars(kappa=c(0.6),fname=ag_file)
write_sc_pars(kappa=c(0.4, 0.5, 0.6),fname=sc_file)

## Run
model_res <- run_simulations(sc_file, sero_file, ag_file, niter=20, label="Sample run")
## Plot relative strain frequencies
p1 <- plot_trajectories(model_res,plot_rel_freq=FALSE)

setup_sim_dir(
  nstrain=9, 
  nhost=5000, 
  strain_id=c(1,2,3,4,5,6,7,8,9), 
  A=c(1,1,1,2,2,2,3,3,3), 
  G=c(1,2,3,1,2,3,1,2,3), 
  O=c(1,1,1,1,1,1,1,1,1),
  V=c(0,0,0,0,0,0,0,0,0), 
  init = c(0,1,0,0,0,1,1,0,0),
  t_max=30*365,
  t_vax=30*365,
  beta=0.09,
  t_g=0.00005,
  t_l=0.00012)

write_sero_pars(alpha=c(0.02,0.025, 0.03),fname=sero_file)
write_ag_pars(kappa=c(0.6),fname=ag_file)
write_sc_pars(kappa=c(0.4, 0.5, 0.6),fname=sc_file)

## Run
model_res <- run_simulations(sc_file, sero_file, ag_file, niter=20, label="Sample run")
## Plot relative strain frequencies
p11 <- plot_trajectories(model_res,plot_rel_freq=FALSE)


setup_sim_dir(
  nstrain=9, 
  nhost=5000, 
  strain_id=c(1,2,3,4,5,6,7,8,9), 
  A=c(1,1,1,2,2,2,3,3,3), 
  G=c(1,2,3,1,2,3,1,2,3), 
  O=c(1,1,1,1,1,1,1,1,1),
  V=c(0,0,0,0,0,0,0,0,0), 
  init = c(1,0,0,0,1,0,0,0,1),
  t_max=30*365,
  t_vax=30*365,
  beta=0.09,
  t_g=0.00005,
  t_l=0.00012)

write_sero_pars(alpha=c(0.02,0.025, 0.03),fname=sero_file)
write_ag_pars(kappa=c(0.6),fname=ag_file)
write_sc_pars(kappa=c(0.5, 0.5, 0.5),fname=sc_file)

## Run
model_res <- run_simulations(sc_file, sero_file, ag_file, niter=20, label="Sample run")
## Plot relative strain frequencies
p2 <- plot_trajectories(model_res,plot_rel_freq=FALSE)



setup_sim_dir(
  nstrain=9, 
  nhost=5000, 
  strain_id=c(1,2,3,4,5,6,7,8,9), 
  A=c(1,1,1,2,2,2,3,3,3), 
  G=c(1,2,3,1,2,3,1,2,3), 
  O=c(1,1,1,1,1,1,1,1,1),
  V=c(0,0,0,0,0,0,0,0,0), 
  init = c(1,0,0,0,1,0,0,0,1),
  t_max=30*365,
  t_vax=30*365,
  beta=0.09,
  t_g=0.00005,
  t_l=0.00012)

write_sero_pars(alpha=c(0.03,0.03, 0.03),fname=sero_file)
write_ag_pars(kappa=c(0.6),fname=ag_file)
write_sc_pars(kappa=c(0.4, 0.5, 0.6),fname=sc_file)

## Run
model_res <- run_simulations(sc_file, sero_file, ag_file, niter=20, label="Sample run")
## Plot relative strain frequencies
p3 <- plot_trajectories(model_res,plot_rel_freq=FALSE)

## when get total exclusion of meta types?

