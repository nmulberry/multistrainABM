#############################
# SIMPLE EXAMPLE ON USING MODEL
#############################
source("setup.R")

res_file <- "res_pars.csv"
sero_file <- "sero_pars.csv"

setup_sim_dir(
  nstrain=2, 
  nhost=10000, 
  strain_id=c(1,2), 
  st=c(1,1), 
  rt=c(0,1), 
  mt=c(1,3),
  vt=c(1,0), 
  t_max=10*365,
  t_vax=10*365,
  beta=0.1,
  mu=0.0)

write_sero_pars(alpha=c(0.025),fname_sero=sero_file)
write_res_pars(p_tau=0.4,cost_res=0.15,tau=0.35,fname_res=res_file)
## Run
model_res <- run_simulations(res_file, sero_file, niter=20, label="Sample run")
## Plot relative strain frequencies
p1 <- plot_trajectories(model_res,plot_rel_freq=TRUE)+ylim(c(0,1))

##############################################

res_file <- "res_pars.csv"
sero_file <- "sero_pars.csv"

setup_sim_dir(
  nstrain=2, 
  nhost=10000, 
  strain_id=c(1,2), 
  st=c(1,2), 
  rt=c(0,1), 
  mt=c(1,1),
  vt=c(1,0), 
  t_max=10*365,
  t_vax=10*365,
  beta=0.1,
  mu=0.00)

write_sero_pars(alpha=c(0.025, 0.03),fname_sero=sero_file)
write_res_pars(p_tau=0.4,cost_res=0.15,tau=0.35,fname_res=res_file)
## Run
model_res <- run_simulations(res_file, sero_file, niter=20, label="Sample run")
## Plot relative strain frequencies
p2 <- plot_trajectories(model_res,plot_rel_freq=TRUE)+ylim(c(0,1))




setup_sim_dir(
  nstrain=2, 
  nhost=10000, 
  strain_id=c(1,2), 
  st=c(1,2), 
  rt=c(0,1), 
  mt=c(1,3),
  vt=c(1,0), 
  t_max=10*365,
  t_vax=10*365,
  beta=0.1,
  mu=0.00)

write_sero_pars(alpha=c(0.025, 0.03),fname_sero=sero_file)
write_res_pars(p_tau=0.4,cost_res=0.15,tau=0.35,fname_res=res_file)
## Run
model_res <- run_simulations(res_file, sero_file, niter=20, label="Sample run")
## Plot relative strain frequencies
p3 <- plot_trajectories(model_res,plot_rel_freq=TRUE)+ylim(c(0,1))


ggarrange(plotlist=list(p1,p2,p3),nrow=3, common.legend=TRUE)

