 #############################
# SAMPLE TRAJECTORIES
#############################


# SAME SERO
res_file <- "res_pars.csv"
sero_file <- "sero_pars.csv"

setup_sim_dir(
  nstrain=2, 
  nhost=5000, 
  strain_id=c(1,2), 
  st=c(1,1), 
  rt=c(0,1), 
  mt=c(1,3),
  vt=c(1,0), 
  t_max=10*365,
  t_vax=10*365,
  beta=BETA,
  mu=0.0)

write_sero_pars(alpha=c(0.025),fname_sero=sero_file)
write_res_pars(p_tau=0.4,cost_res=0.15,tau=0.35,fname_res=res_file)
## Run
model_res1 <- run_simulations(res_file, sero_file, niter=20, label="Antigenic")


## Run trajectories

traj_res1 <- run_trajectories(res_file, sero_file) %>%
  reshape2::melt(id="time") %>%
  mutate(label="Antigenic") 


##################################
##################################
setup_sim_dir(
  nstrain=2, 
  nhost=5000, 
  strain_id=c(1,2), 
  st=c(1,2), 
  rt=c(0,1), 
  mt=c(1,1),
  vt=c(1,0), 
  t_max=10*365,
  t_vax=10*365,
  beta=BETA,
  mu=0.0)

write_sero_pars(alpha=c(0.025, 0.03),fname_sero=sero_file)
write_res_pars(p_tau=0.4,cost_res=0.15,tau=0.35,fname_res=res_file)
## Run
model_res2 <- run_simulations(res_file, sero_file, niter=20, label="Metabolic")

## Run trajectories

traj_res2 <- run_trajectories(res_file, sero_file) %>%
  reshape2::melt(id="time") %>%
  mutate(label="Metabolic")



##################################
##################################
setup_sim_dir(
  nstrain=2, 
  nhost=5000, 
  strain_id=c(1,2), 
  st=c(1,2), 
  rt=c(0,1), 
  mt=c(1,3),
  vt=c(1,0), 
  t_max=10*365,
  t_vax=10*365,
  beta=BETA,
  mu=0.0)

write_sero_pars(alpha=c(0.025,0.03),fname_sero=sero_file)
write_res_pars(p_tau=0.4,cost_res=0.15,tau=0.35,fname_res=res_file)
## Run
model_res3 <- run_simulations(res_file, sero_file, niter=20, label="None")
## Run trajectories

traj_res3 <- run_trajectories(res_file, sero_file) %>%
  reshape2::melt(id="time") %>%
  mutate(label="None")



##################


setup_sim_dir(
  nstrain=2, 
  nhost=5000, 
  strain_id=c(1,2), 
  st=c(1,1), 
  rt=c(0,1), 
  mt=c(1,1),
  vt=c(1,0), 
  t_max=10*365,
  t_vax=10*365,
  beta=BETA,
  mu=0.0)

write_sero_pars(alpha=c(0.025),fname_sero=sero_file)
write_res_pars(p_tau=0.4,cost_res=0.15,tau=0.35,fname_res=res_file)
## Run
model_res4 <- run_simulations(res_file, sero_file, niter=20, label="Both")


## Run trajectories

traj_res4 <- run_trajectories(res_file, sero_file) %>%
  reshape2::melt(id="time") %>%
  mutate(label="Both")

traj_res <- rbind(traj_res1, traj_res2, traj_res3, traj_res4)
model_res <- rbind(model_res1, model_res2, model_res3, model_res4)
## PLOTS
traj_res$label <- factor(traj_res$label, levels=c("None", "Metabolic", "Antigenic", "Both"))
model_res$label <- factor(model_res$label, levels=c("None", "Metabolic", "Antigenic", "Both"))


p1 <- plot_trajectories(model_res, plot_rel_freq=FALSE)+
  scale_fill_manual(values=sens_res_pal)+
  scale_color_manual(values=sens_res_pal) +
  facet_grid(label~.)+
  labs(x="Time (yrs)", y="Prevalence", col="")

p2 <- ggplot(traj_res, aes(x=time, y=value, col=variable))+
  geom_line(size=1.25)+
  facet_grid(~label)+
  scale_color_manual(values=sens_res_pal) +
  facet_grid(label~.)+
  theme(strip.text.y = element_blank())+
  labs(x="Time (yrs)", y="Density", col="")
##################
ggarrange(plotlist=list(p2, p1), ncol=2, common.legend=TRUE, legend="bottom", labels=c("A", "B"),
    font.label=list(size=16))












