#############################
# SIMPLE EXAMPLE ON USING MODEL
#############################
source("setup.R")

sc_file <- "sc_pars.csv"
sero_file <- "sero_pars.csv"
res_file <- "res_pars.csv" 


# init = c(0,1,0,0,0,0,0,1)
setup_sim_dir(
  nstrain=4, 
  nhost=5000, 
  strain_id=c(1,2,3,4), 
  A=c(1,1,2,2), 
  G=c(1,2,1,2), 
  R=c(0,0,0,0),
  V=c(1,1,0,0), 
  init = c(0,1,1,0),
  t_max=30*365,
  t_vax=20*365,
  beta=0.09,
  t_g=0.00005,
  t_l=0.00005)

write_sero_pars(alpha=c(0.02,0.025),fname=sero_file)
write_sc_pars(kappa=c(1.1,1.11),fname=sc_file)
write_res_pars(p_tau=0.2,cost_res=0.12,tau=0.8,fname_res=res_file)
## Run
model_res <- run_simulations(sc_file, sero_file, res_file, niter=20, label="sim")

#################3
## Plot relative strain frequencies
p1 <- plot_trajectories(model_res,plot_rel_freq=FALSE)

## 

##summarize
strain_pars <- read.csv("strain_pars.csv")
strain_pars$id <- paste0("strain ", seq(1,nrow(strain_pars)))

dat <- model_res %>%
      reshape2::melt(id=c("time", "iter","label"), 
        variable.name="id", 
        value.name="freq") %>%
      inner_join(strain_pars,by="id") %>%
      filter(time == max(time)) %>%
      group_by(label, id_m, id_s,id_r) %>%
      summarize(freq=median(freq))

dat$id_r <- as.factor(dat$id_r)
dat$id_m <- as.factor(dat$id_m)
dat$id_s <- as.factor(dat$id_s)

g1 <- ggplot(dat, aes(x=id_m, y=id_s, fill=freq))+
  geom_tile()+
  facet_grid(label~id_r)













# init = c(0,1,0,0,0,0,0,1)
setup_sim_dir(
  nstrain=9, 
  nhost=5000, 
  strain_id=c(1,2,3,4,5,6,7,8,9), 
  A=c(1,1,1,2,2,2,3,3,3), 
  G=c(1,2,3,1,2,3,1,2,3), 
  R=c(0,0,0,0,0,0,0,0,0),
  V=c(1,1,1,0,0,0,0,0,0), 
  init = c(0,0,1,0,1,0,1,0,0),
  t_max=30*365,
  t_vax=20*365,
  beta=0.09,
  t_g=0.00005,
  t_l=0.00005)

write_sero_pars(alpha=c(0.02,0.025,0.035),fname=sero_file)
write_sc_pars(kappa=c(1.1,1.12,1.14),fname=sc_file)
write_res_pars(p_tau=0.2,cost_res=0.12,tau=0.8,fname_res=res_file)
## Run
model_res <- run_simulations(sc_file, sero_file, res_file, niter=20, label="sim")

#################3
## Plot relative strain frequencies
p2 <- plot_trajectories(model_res,plot_rel_freq=FALSE)


