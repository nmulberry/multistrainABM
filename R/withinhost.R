

T_INV <- 10 # set invasion time
rho <- 0.01
# Look at the four different competition types: (with and without treatment)
# Resistant strain is always invading-type
# For each run, also look at the resulting pop-level trajectories for "high treament" and "low treatment" cases
res_file <- "res_pars.csv"
sero_file <- "sero_pars.csv"

#---NO COMPETITION
res1_notreat <- run_one_host(
  nstrain=2,
  strain_id=c(1,2),
  st=c(1,2), 
  rt=c(0,1), 
  mt=c(1,3),
  vt=c(0,0), 
  alpha=c(0.03), 
  p_tau=0, 
  cost_res=COST, 
  tau=TAU, 
  t_max=300,
  t_vax=300,
  t_inv=10,
  init_strain=c(1)
) %>%
  dplyr::select(-starts_with("I")) %>%
  reshape2::melt(id=c('time')) %>%
  mutate(value2=pmax(0,value-rho)) %>%
  group_by(time) %>%
  mutate(tot=sum(value2)) %>%
  ungroup() %>%
  mutate(rel_freq = value2/tot) %>%
  mutate(treat="Untreated") %>%
  mutate(comp="(a)")

# pop-level simulations
setup_sim_dir(
  nstrain=2, 
  nhost=2000, 
  strain_id=c(1,2), 
  st=c(1,2), 
  rt=c(0,1), 
  mt=c(1,5),
  vt=c(0,0), 
  t_max=5*365,
  t_vax=5*365,
  beta=BETA)

write_sero_pars(alpha=c(0.03,0.03), fname_sero=sero_file)
write_res_pars(p_tau=0.7,cost_res=COST,tau=TAU,fname_res=res_file)
## Run
res1_hightreat <- run_simulations(res_file, sero_file, niter=no_cores*3, label="(a)")%>%
  mutate(treat="High Treatment")
write_res_pars(p_tau=0.2,cost_res=COST,tau=TAU,fname_res=res_file)
## Run
res1_lowtreat <- run_simulations(res_file, sero_file, niter=no_cores*3,label="(a)")%>%
  mutate(treat="Low Treatment")




res1_treat <- run_one_host(
  nstrain=2,
  strain_id=c(1,2),
  st=c(1,2), 
  rt=c(0,1), 
  mt=c(1,3),
  vt=c(0,0), 
  alpha=c(0.03), 
  p_tau=1, 
  cost_res=COST, 
  tau=TAU, 
  t_max=300,
  t_vax=300,
  t_inv=10,
  init_strain=c(1)
) %>%
  dplyr::select(-starts_with("I")) %>%
  reshape2::melt(id=c('time')) %>%
  mutate(value2=pmax(0,value-rho)) %>%
  group_by(time) %>%
  mutate(tot=sum(value2)) %>%
  ungroup() %>%
  mutate(rel_freq = value2/tot)%>%
  mutate(treat="Treated")%>%
  mutate(comp="(a)")


#----WITHIN HOST COMPETITION
res2_notreat <- run_one_host(
  nstrain=2,
  strain_id=c(1,2),
  st=c(1,2), 
  rt=c(0,1), 
  mt=c(1,1),
  vt=c(0,0), 
  alpha=c(0.03), 
  p_tau=0, 
  cost_res=COST, 
  tau=TAU, 
  t_max=300,
  t_vax=300,
  t_inv=T_INV,
  init_strain=c(1)
)  %>%
  dplyr::select(-starts_with("I")) %>%
  reshape2::melt(id=c('time')) %>%
  mutate(value2=pmax(0,value-rho)) %>%
  group_by(time) %>%
  mutate(tot=sum(value2)) %>%
  ungroup() %>%
  mutate(rel_freq = value2/tot)%>%
  mutate(treat="Untreated")%>%
  mutate(comp="(b)")

res2_treat <- run_one_host(
  nstrain=2,
  strain_id=c(1,2),
  st=c(1,2), 
  rt=c(0,1), 
  mt=c(1,1),
  vt=c(0,0), 
  alpha=c(0.03), 
  p_tau=1, 
  cost_res=COST, 
  tau=TAU, 
  t_max=300,
  t_vax=300,
  t_inv=T_INV,
  init_strain=c(1)
) %>%
  dplyr::select(-starts_with("I")) %>%
  reshape2::melt(id=c('time')) %>%
  mutate(value2=pmax(0,value-rho)) %>%
  group_by(time) %>%
  mutate(tot=sum(value2)) %>%
  ungroup() %>%
  mutate(rel_freq = value2/tot)%>%
  mutate(treat="Treated")%>%
  mutate(comp="(b)")

# pop-level simulations
setup_sim_dir(
  nstrain=2, 
  nhost=2000, 
  strain_id=c(1,2), 
  st=c(1,2), 
  rt=c(0,1), 
  mt=c(1,1),
  vt=c(0,0), 
  t_max=5*365,
  t_vax=5*365,
  beta=BETA)

write_sero_pars(alpha=c(0.03,0.03), fname_sero=sero_file)
write_res_pars(p_tau=0.7,cost_res=COST,tau=TAU,fname_res=res_file)
## Run
res2_hightreat <- run_simulations(res_file, sero_file, niter=no_cores*3, label="(b)")%>%
  mutate(treat="High Treatment")

write_res_pars(p_tau=0.2,cost_res=COST,tau=TAU,fname_res=res_file)
## Run
res2_lowtreat <- run_simulations(res_file, sero_file, niter=no_cores*3,label="(b)")%>%
  mutate(treat="Low Treatment")


#----BETWEEN HOST COMPETITION 
res3_notreat <- run_one_host(
  nstrain=2,
  strain_id=c(1,2),
  st=c(1,1), 
  rt=c(0,1), 
  mt=c(1,3),
  vt=c(0,0), 
  alpha=c(0.03), 
  p_tau=0, 
  cost_res=COST, 
  tau=TAU, 
  t_max=300,
  t_vax=300,
  t_inv=T_INV,
  init_strain=c(1)
)  %>%
  dplyr::select(-starts_with("I")) %>%
  reshape2::melt(id=c('time')) %>%
  mutate(value2=pmax(0,value-rho)) %>%
  group_by(time) %>%
  mutate(tot=sum(value2)) %>%
  ungroup() %>%
  mutate(rel_freq = value2/tot)%>%
  mutate(treat="Untreated")%>%
  mutate(comp="(c)")

res3_treat <- run_one_host(
  nstrain=2,
  strain_id=c(1,2),
  st=c(1,1), 
  rt=c(0,1), 
  mt=c(1,3),
  vt=c(0,0), 
  alpha=c(0.03), 
  p_tau=1, 
  cost_res=COST, 
  tau=TAU, 
  t_max=300,
  t_vax=300,
  t_inv=T_INV,
  init_strain=c(1)
) %>%
  dplyr::select(-starts_with("I")) %>%
  reshape2::melt(id=c('time')) %>%
  mutate(value2=pmax(0,value-rho)) %>%
  group_by(time) %>%
  mutate(tot=sum(value2)) %>%
  ungroup() %>%
  mutate(rel_freq = value2/tot)%>%
  mutate(treat="Treated")%>%
  mutate(comp="(c)")

# pop-level simulations
setup_sim_dir(
  nstrain=2, 
  nhost=2000, 
  strain_id=c(1,2), 
  st=c(1,1), 
  rt=c(0,1), 
  mt=c(1,2),
  vt=c(0,0), 
  t_max=5*365,
  t_vax=5*365,
  beta=BETA)

write_sero_pars(alpha=c(0.03), fname_sero=sero_file)
write_res_pars(p_tau=0.7,cost_res=COST,tau=TAU,fname_res=res_file)
## Run
res3_hightreat <- run_simulations(res_file, sero_file, niter=no_cores*3, label="(c)")%>%
  mutate(treat="High Treatment")

write_res_pars(p_tau=0.2,cost_res=COST,tau=TAU,fname_res=res_file)
## Run
res3_lowtreat <- run_simulations(res_file, sero_file, niter=no_cores*3,label="(c)")%>%
  mutate(treat="Low Treatment")

#-----BOTH
res4_notreat <- run_one_host(
  nstrain=2,
  strain_id=c(1,2),
  st=c(1,1), 
  rt=c(0,1), 
  mt=c(1,1),
  vt=c(0,0), 
  alpha=c(0.03), 
  p_tau=0, 
  cost_res=COST, 
  tau=TAU, 
  t_max=300,
  t_vax=300,
  t_inv=T_INV,
  init_strain=c(1)
) %>%
  dplyr::select(-starts_with("I")) %>%
  reshape2::melt(id=c('time')) %>%
  mutate(value2=pmax(0,value-rho)) %>%
  group_by(time) %>%
  mutate(tot=sum(value2)) %>%
  ungroup() %>%
  mutate(rel_freq = value2/tot)%>%
  mutate(treat="Untreated")%>%
  mutate(comp="(d)")

res4_treat <- run_one_host(
  nstrain=2,
  strain_id=c(1,2),
  st=c(1,1), 
  rt=c(0,1), 
  mt=c(1,1),
  vt=c(0,0), 
  alpha=c(0.03), 
  p_tau=1, 
  cost_res=COST, 
  tau=TAU, 
  t_max=300,
  t_vax=300,
  t_inv=T_INV,
  init_strain=c(1)
) %>%
  dplyr::select(-starts_with("I")) %>%
  reshape2::melt(id=c('time')) %>%
  mutate(value2=pmax(0,value-rho)) %>%
  group_by(time) %>%
  mutate(tot=sum(value2)) %>%
  ungroup() %>%
  mutate(rel_freq = value2/tot)%>%
  mutate(treat="Treated")%>%
  mutate(comp="(d)")



# pop-level simulations
setup_sim_dir(
  nstrain=2, 
  nhost=2000, 
  strain_id=c(1,2), 
  st=c(1,1), 
  rt=c(0,1), 
  mt=c(1,1),
  vt=c(0,0), 
  t_max=5*365,
  t_vax=5*365,
  beta=BETA)

write_sero_pars(alpha=c(0.03), fname_sero=sero_file)
write_res_pars(p_tau=0.7,cost_res=COST,tau=TAU,fname_res=res_file)
## Run
res4_hightreat <- run_simulations(res_file, sero_file, niter=no_cores*3, label="(d)")%>%
  mutate(treat="High Treatment")

write_res_pars(p_tau=0.2,cost_res=COST,tau=TAU,fname_res=res_file)
## Run
res4_lowtreat <- run_simulations(res_file, sero_file, niter=no_cores*3,label="(d)") %>%
  mutate(treat="Low Treatment")




N <- 60


##### PLOT
res <- rbind(res1_notreat, res1_treat, res2_notreat, res2_treat, res3_notreat,
  res3_treat, res4_notreat, res4_treat) %>%
  mutate(variable=case_when(variable=="X1" ~ "strain 1", TRUE ~ "strain 2"))

#res$comp <- factor(res$comp, levels=c("No Competition", "Metabolic",
#  "Antigenic", "Both"))

gg_traj_treat <- ggplot(filter(res, treat=="Treated"),
   aes(x=time, y=value, col=variable))+
    scale_y_continuous(breaks=c(0,0.5,1), limits=c(0,1),expand = c(0, 0))+
    scale_x_continuous(breaks=c(0,N-10), limits=c(0,N),expand = c(0, 0))+
  geom_line(size=1.5)+
  facet_wrap(~comp, nrow=1)+
  scale_color_manual(values=sens_res_pal)+
  labs(x="Time (days)", y="Density", col="") +
#  theme(strip.text.x = element_blank())+
  ggtitle("Treated")+
  theme(legend.position="none")

gg_traj_notreat <- ggplot(filter(res, treat=="Untreated"),
 aes(x=time, y=value, col=variable))+
  scale_y_continuous(breaks=c(0,0.5,1), limits=c(0,1),expand = c(0, 0))+
  scale_x_continuous(breaks=c(0,N-10), limits=c(0,N),expand = c(0, 0))+
  geom_line(size=1.5)+
  facet_wrap(~comp,nrow=1)+
  scale_color_manual(values=sens_res_pal)+
  labs(x="Time (days)", y="Density", col="") +
 # theme(strip.text.x = element_blank())+
  ggtitle("Untreated")+
  theme(legend.position="none")


pop_res_hightreat <- rbind(res1_hightreat, res2_hightreat,
  res3_hightreat,res4_hightreat) %>%
  dplyr::select(-c(treat)) 

#pop_res_hightreat$label <- factor(pop_res_hightreat$label , 
#  levels=c("No Competition", "Metabolic",
#  "Antigenic", "Both"))


pop_res_lowtreat <- rbind(res1_lowtreat, res2_lowtreat,
  res3_lowtreat,res4_lowtreat)%>%
  dplyr::select(-c(treat))

#pop_res_lowtreat$label <- factor(pop_res_lowtreat$label , 
#  levels=c("No Competition", "Metabolic",
#  "Antigenic", "Both"))


gg_hightreat <- plot_trajectories(filter(pop_res_hightreat, time > 2*265), 
  plot_rel_freq=FALSE)+
  facet_wrap(~label, nrow=1)+
  scale_color_manual(values=sens_res_pal)+
  scale_fill_manual(values=sens_res_pal)+
  labs(x="Time (yrs)", y="Prevalence", col="")+
  theme(legend.position="none")+
  ggtitle(parse(text = paste0('"High treatment"', ' ~ (p[tau] == ', 0.7, ')')))+
 # theme(strip.text.x = element_blank())+
  scale_y_continuous(breaks=c(0,0.15,0.3),limits=c(0,0.3))+
  scale_x_continuous(breaks=c(2,5))

gg_lowtreat <- plot_trajectories(filter(pop_res_lowtreat, time > 2*265), plot_rel_freq=FALSE)+
  facet_wrap(~label, nrow=1)+
  scale_color_manual(values=sens_res_pal)+
  scale_fill_manual(values=sens_res_pal)+
  labs(x="Time (yrs)", y="Prevalence", col="")+
  theme(legend.position="none")+
  ggtitle(parse(text = paste0('"Low treatment"', ' ~ (p[tau] == ', 0.2, ')')))+
#  theme(strip.text.x = element_blank())+
  scale_y_continuous(breaks=c(0,0.15,0.3),limits=c(0,0.3))+
  scale_x_continuous(breaks=c(2,5))

gg_poplevel <- ggarrange(plotlist=list(gg_lowtreat, gg_hightreat),nrow=1)




# PUT FIGS TOGETHER
res2 <- res%>%
  mutate(comp=case_when(comp=="(a)" ~ "(a) No Competition",
    comp=="(b)"~"(b) Metabolic",
    comp=="(c)" ~ "(c) Antigenic",
    comp=="(d)"~"(d) Both"))
gg_system <- ggplot(res2) +
  facet_wrap(~ comp, nrow=1, strip.position="bottom")

## And plot transmission probabilities...
gg_probs_treat <- ggplot(filter(res, treat=="Treated"),
  aes(x=time, y=rel_freq, fill=variable))+
  geom_area()+
    scale_x_continuous(limits = c(0,N), expand = c(0, 0),breaks=c(0,N-10))+
    scale_y_continuous(limits = c(0,1), expand = c(0, 0),breaks=c(0,0.5,1))+
  facet_wrap(~comp,nrow=1)+
  scale_fill_manual(values=sens_res_pal)+
  labs(y='P(transmit)', x='Time (days)', fill="") +
  #theme(strip.text.x = element_blank())+
  theme(legend.position="none")  + ggtitle("Treated")

gg_probs_notreat <- ggplot(filter(res, treat=="Untreated"),
  aes(x=time, y=rel_freq, fill=variable))+
  geom_area()+
    scale_x_continuous(limits = c(0,N), expand = c(0, 0),breaks=c(0,N-10))+
    scale_y_continuous(limits = c(0,1), expand = c(0, 0),breaks=c(0,0.5,1))+
  facet_wrap(~comp,nrow=1)+
  scale_fill_manual(values=sens_res_pal)+
  labs(y='P(transmit)', x='Time (days)', fill="") +
 # theme(strip.text.x = element_blank())+
  theme(legend.position="none")  + ggtitle("Untreated")



gg_res <- ggarrange(plotlist=list(
  ggarrange(plotlist=list(ggplot()+theme_minimal(), gg_system, ggplot()+theme_minimal()), nrow=1,  widths=c(0.2,1,0.2)), 
  ggarrange(plotlist=list(
      gg_traj_notreat,
      gg_traj_treat)),
 # ggarrange(plotlist=list(
 #     gg_probs_notreat,
 #     gg_probs_treat)),
  gg_poplevel), nrow=3, labels=c("A","B","C"), heights=c(0.6,1,1),
  font.label=list(size=16))

