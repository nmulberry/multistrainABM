#############################
# SIMPLE EXAMPLE ON USING MODEL
#############################
source("setup.R")

res_file <- "res_pars.csv"
sero_file <- "sero_pars.csv"

setup_sim_dir(
  nstrain=4, 
  nhost=4000, 
  strain_id=c(1,2,3,4), 
  st=c(1,1,2,2), 
  rt=c(0,0,0,0), 
  mt=c(2,1,1,3),
  vt=c(1,1,0,0), 
  gt=c(1,1,1,1),
  t_max=10*365,
  t_vax=5*365,
  beta=0.1,
  mu=0.0,
  k_g=0.05)

write_sero_pars(alpha=c(0.03,0.03),fname_sero=sero_file) #strain 3 is ref. strain
write_res_pars(p_tau=0.4,cost_res=0.15,tau=0.35,fname_res=res_file)
## Run
model_res <- run_simulations(res_file, sero_file, niter=8, label="kg=0.05")
p1 <- plot_trajectories(model_res,plot_rel_freq=FALSE)

setup_sim_dir(
  nstrain=3, 
  nhost=4000, 
  strain_id=c(1,2,3), 
  st=c(1,2,2), 
  rt=c(0,0,0), 
  mt=c(1,1,3),
  vt=c(1,0,0), 
  gt=c(1,1,1),
  t_max=5*365,
  t_vax=5*365,
  beta=0.1,
  mu=0.0,
  k_g=0.25)
write_sero_pars(alpha=c(0.03,0.031),fname_sero=sero_file) #strain 3 is ref. strain
write_res_pars(p_tau=0.4,cost_res=0.15,tau=0.35,fname_res=res_file)
## Run
model_res <- run_simulations(res_file, sero_file, niter=8, label="kg=0.3")
## Plot relative strain frequencies
p2 <- plot_trajectories(model_res,plot_rel_freq=FALSE)


setup_sim_dir(
  nstrain=3, 
  nhost=4000, 
  strain_id=c(1,2,3), 
  st=c(1,2,2), 
  rt=c(0,0,0), 
  mt=c(1,1,3),
  vt=c(1,0,0), 
  gt=c(1,1,1),
  t_max=5*365,
  t_vax=0*365,
  beta=0.1,
  mu=0.0,
  k_g=0.25)
write_sero_pars(alpha=c(0.03,0.031),fname_sero=sero_file) #strain 3 is ref. strain
write_res_pars(p_tau=0.4,cost_res=0.15,tau=0.35,fname_res=res_file)
## Run
model_res <- run_simulations(res_file, sero_file, niter=8, label="kg=0.3")
## Plot relative strain frequencies
p3 <- plot_trajectories(model_res,plot_rel_freq=FALSE)

ggarrange(plotlist=list(p2,p3), common.legend=TRUE)
