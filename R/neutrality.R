
four_pal <- c("strain 1"="#FAD510", "strain total"="#273046")

############################
# RUN TWO IDENTICAL STRAINS
############################
res_file <- "res_pars.csv"
sero_file <- "sero_pars.csv"

# ---- First look at two identical strains (same serotype, metatype and resistance type)
setup_sim_dir(
  nstrain=2, 
  nhost=10000, 
  strain_id=c(1,2), 
  st=c(1,1), 
  rt=c(0,0), 
  mt=c(1,1),
  vt=c(0,0), 
  t_max=T_MAX,
  t_vax=T_MAX,
  beta=BETA)
write_sero_pars(alpha=c(0.03),fname_sero=sero_file)
write_res_pars(p_tau=P_TAU,cost_res=COST,tau=TAU,fname_res=res_file)
# Run
model_res1 <- run_simulations(res_file, sero_file, niter=NITER, label="Two Identical Strains")
#gg_identical <- plot_trajectories(model_res1,plot_rel_freq=FALSE, with_total=TRUE)+
#  labs(x="Time (yrs)", y="Prevalence", col="")+scale_color_manual(values=four_pal)+
#  scale_fill_manual(values=four_pal)+ggtitle("Identical Strains")


model_res1_tot <- model_res1 %>% 
  group_by(iter, label, time) %>%
  summarize(`strain total` = `strain 1`+`strain 2`) 


# --- Now look at one strain
setup_sim_dir(
  nstrain=1, 
  nhost=10000, 
  strain_id=1, 
  st=1, 
  rt=0, 
  mt=1,
  vt=0, 
  t_max=T_MAX,
  t_vax=T_MAX,
  beta=BETA)
write_sero_pars(alpha=c(0.03),fname_sero=sero_file)
write_res_pars(p_tau=P_TAU,cost_res=COST,tau=TAU,fname_res=res_file)
# Run
model_res2 <- run_simulations(res_file, sero_file, niter=NITER, label="One Strain") %>%
  rename(`strain total`=`strain 1`)

model_res <- rbind(model_res1_tot, model_res2)

gg_poplevel <- plot_trajectories(model_res,plot_rel_freq=FALSE)+
  labs(x="Time (yrs)", y="Prevalence", col="")+scale_color_manual(values=four_pal)+
  scale_fill_manual(values=four_pal)+theme(legend.position="none")+
  facet_wrap(~label, nrow=2)




## initial relative frequencies preserved?

run_over_f <- function(f){
      X0 <- rep(0, 2) 
      I0 <- rep(0, 2)
      tot <- 0.1
      X0[1] <- (1-f)*tot
      X0[2] <- f*tot
      y0 <- c(X=X0, I=I0)
      kappa=1.1 # growth rate of sensitive strain
      eps_x=0.03 # wh death rate
      eps_a=0.02 # immune decay rate
      rho=0.01 # minimum wh pathogen density
      theta=70
      parameters <- c(nstrain=2,eps_a=eps_a,
        eps_x=eps_x,rho=rho, theta=theta,
        kappa=kappa, alpha=ALPHA, 
        tau=TAU, st=c(1,1), mt=c(1,1),rt=c(0,0),
        p_tau=P_TAU,cost_res=COST)

      times <- seq(0, 100, by = 1)
      out1 <- data.frame(ode(y=y0, times=times, func=withinhost_model, 
        parms=parameters) )   %>%
      dplyr::select(-starts_with("I")) %>%
      reshape2::melt(id='time') %>%
      mutate(frac=f)

      return(out1)
}


res <- seq(0.0, 0.5, by=0.1) %>% map_dfr(run_over_f)

g1 <- ggplot(res, aes(x=time, y=value, col=variable))+
  geom_line(size=1.5)+
  facet_wrap(~frac, nrow=1)+
  labs(x="Time (days)", y="Density", col="")+
  scale_color_manual(values=c("#E69F00","#56B4E9"))+
  theme(strip.text.x=element_blank(), legend.position="top")+
  scale_x_continuous(breaks=c(100))



res2 <- res %>%
  group_by(frac,time) %>%
  summarize(total=sum(value))

g2 <- ggplot(res2, aes(x=time, y=total))+
  geom_line(size=1.5)+
  facet_wrap(~frac, nrow=1)+
  labs(x="Time (days)", y="Density")+
  theme(strip.text.x=element_blank())+
  scale_x_continuous(breaks=c(100))


res3 <- res %>%
  group_by(frac,time) %>%
  mutate(total=sum(value)) %>%
  filter(variable=="X2") %>%
  mutate(frac2 = value/total) 

  

g3 <- ggplot(res3, aes(x=time, y=frac2))+
  geom_line(size=1.5)+
  facet_wrap(~frac, nrow=1)+
 labs(x="Time (days)", y="Relative Freq") +
   theme(strip.text.x=element_blank())+
  scale_x_continuous(breaks=c(100))

gg_wh <- ggarrange(plotlist=list(g1,g2,g3), nrow=3, labels=c("B","C","D"), font.label=c(size=16))

gg_neutral <- ggarrange(plotlist=list(gg_poplevel, gg_wh), labels=c("A",""), font.label=c(size=16),
  widths=c(1,2))
