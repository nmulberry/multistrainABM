


################################
# LOOK AT SENS/RES COEXISTENCE
#################################
ptau_vec <- seq(0.3,0.7,by=0.2)
tau_vec <- seq(0.0,1.0, by=0.1)
cost_vec <- seq(0.0,0.5,by=0.05)

res_file <- "res_pars.csv"
sero_file <- "sero_pars.csv"

pars <- expand.grid(alpha1=ALPHA, 
  alpha2=ALPHA, 
  comp_type=c("Antigenic", "Both"),
  rt1=1, 
  rt2=0, 
  p_tau=ptau_vec, 
  cost=cost_vec,
  tau=tau_vec,
  t_max=T_MAX, 
  t_vax=T_MAX)


#res_df <- pars %>% pmap(get_rel_freq_2strain)


#res_df <- do.call("rbind", res_df)
#res <- cbind(pars, res_df)

res <- readRDS("../generated/res_ptau.RDS")

gg_res <- ggplot(res, aes(x=tau, y=cost, fill=`strain 1`, col=`strain 1`))+
  geom_tile()+
  facet_grid(comp_type~p_tau, labeller=label_bquote(cols=p[tau]: .(p_tau), rows=.(comp_type)))+
  scale_fill_viridis(guide= guide_colourbar(barwidth=8,ticks=FALSE))+theme(legend.title=element_text(size=14))+
 scale_colour_viridis(guide= "none")+
  labs(x=expression(Treatment~rate~(tau)), y="Cost Resistance (c)", fill="Relative Freq Resistance")+
  scale_y_continuous(expand=c(0,0))+
  scale_x_continuous(expand=c(0,0))



##############################
# AS FUNCTION OF ALPHA
##############################
alpha_vec <- c(0.02, 0.025, 0.03)
tau_vec <- seq(0.1,1.0, by=0.2)
cost_vec <- seq(0.1,0.5,by=0.1)


sero_file <- "sero_pars.csv"

pars <- expand.grid(alpha1=alpha_vec, 
  alpha2=0, #not needed
  comp_type=c("Antigenic", "Both"),
  rt1=1, 
  rt2=0, 
  p_tau=0.4, 
  cost=cost_vec,
  tau=tau_vec,
  t_max=T_MAX, 
  t_vax=T_MAX)


#res_df <- pars %>% pmap(get_rel_freq_2strain)


#res_df <- do.call("rbind", res_df)
#res <- cbind(pars, res_df)

res <- readRDS("../generated/res_alpha.RDS")
gg_res_alpha <- ggplot(res, aes(x=tau, y=cost, fill=`strain 1`, col=`strain 1`))+
  geom_tile()+
  facet_grid(comp_type~alpha1, labeller=label_bquote(cols=alpha: .(alpha1), rows=.(comp_type)))+
  scale_fill_viridis(guide= guide_colourbar(barwidth=8,ticks=FALSE))+theme(legend.title=element_text(size=14))+
 scale_colour_viridis(guide= "none")+
  labs(x=expression(Treatment~rate~(tau)), y="Cost Resistance (c)", fill="Relative Freq Resistance")+
  scale_y_continuous(expand=c(0,0))+
  scale_x_continuous(expand=c(0,0))



