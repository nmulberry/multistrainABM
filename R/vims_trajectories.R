###############################################
# Investigate Vaccine-Induced Metabolic Shifts
###############################################
# set palette
three_pal <- c("strain 1"="#273046", "strain 2"="#00A08A","strain 3"="#FAD510")
two_pal <- c("strain 1"="#273046", "strain 2"="#00A08A")
t_max <- 15*365
t_vax <- 10*365
#############################
# ALL SENS
#############################

res_file <- "res_pars.csv"
sero_file <- "sero_pars.csv"

setup_sim_dir(
    nstrain=3, 
    nhost=10000, 
    strain_id=c(1,2,3), 
    st=c(1,2,2), 
    rt=c(0,0,0), 
    mt=c(1,1,2),
    vt=c(1,0,0), 
    t_max=t_max,
    t_vax=t_vax,
    beta=BETA)

write_sero_pars(alpha=c(0.025,0.02),fname_sero=sero_file) 
write_res_pars(p_tau=P_TAU,cost_res=0.2,tau=0.3,fname_res=res_file)


## Run
v0_res <- run_simulations(res_file, sero_file, niter=NITER, label="Metabolic")
# first add strain_pars
strain_pars <- read.csv("strain_pars.csv")
strain_pars$id <- paste0("strain ", c(1,2,3))

v0_res_summ <- v0_res %>%
  reshape2::melt(id=c("time", "iter","label")) %>%
  rename(id=variable) %>%
  inner_join(strain_pars,by="id")%>%
  mutate(era = case_when(time >= t_vax-10  & time < t_vax ~ "pre-vax", 
                         time >= t_max-10 ~ "post-vax",
                          TRUE ~ "NA")) %>%
  filter(era != "NA") %>%
  group_by(era, id, id_m, id_s, id_r, id_v, iter,label) %>%
  summarize(freq=mean(value)) %>%
  group_by(era, iter,label) %>%
  nest() %>%
  summarize(map_dfr(data, function(df){
    df_new <- data.frame(
      "NVT" = sum(filter(df, id_v==0)$freq),
      "MT1" = sum(filter(df, id_m==1)$freq),
      "Resistance" = sum(filter(df, id_r==1)$freq))
    return(df_new)
  })) %>%
  reshape2::melt(id=c("era","iter","label")) %>%
  group_by(variable,era,label) %>%
  summarize(p50 = quantile(value, prob=0.5),
          p05 = quantile(value, prob=0.05),
          p95 = quantile(value, prob=0.95))



## Compare to no comp
setup_sim_dir(
  nstrain=3, 
    nhost=10000, 
    strain_id=c(1,2,3), 
    st=c(1,2,2), 
    rt=c(0,0,0), 
    mt=c(1,3,2),
    vt=c(1,0,0), 
    t_max=t_max,
    t_vax=t_vax,
    beta=BETA)

## Run
v0_res2 <- run_simulations(res_file, sero_file, niter=NITER, label="No Competition")

## Plot change in frequency pre-post vax
# first add strain_pars
strain_pars <- read.csv("strain_pars.csv")
strain_pars$id <- paste0("strain ", c(1,2,3))


## PUT TOGETHER FOR PLOTTING
v0_res <- rbind(v0_res, v0_res2)
v0_res$label <- factor(v0_res$label, levels=c("No Competition", "Metabolic"))

## PUT TOGETHER FOR PLOTTING
v0_res <- rbind(v0_res, v0_res2)
v0_res$label <- factor(v0_res$label, levels=c("No Competition", "Metabolic"))
v0_res <- v0_res %>% filter(time > 5*365)
#v0_res <- filter(v0_res, time > 350)
## Plot relative strain frequencies
gg_traj1 <- plot_trajectories(v0_res,plot_rel_freq=FALSE)+
  labs(y="Prevalence", x="Time (years)", col="")+
   scale_color_manual(values=c("#B6854D", "#02401B", "#A2A475"))+
   scale_fill_manual(values=c("#B6854D", "#02401B", "#A2A475"))+
  facet_wrap(~label, ncol=1,strip.position="right")+
  labs(y="Prevalence")+
  ggtitle("")



gg_alpha <- gg_traj1+theme(strip.text.y = element_blank())  


######################
#  ADD RESISTANCE
#####################


setup_sim_dir(
    nstrain=3, 
    nhost=10000, 
    strain_id=c(1,2,3), 
    st=c(1,2,2), 
    rt=c(0,0,1), 
    mt=c(1,1,2),
    vt=c(1,0,0), 
    t_max=t_max,
    t_vax=t_vax,
    beta=BETA)

write_sero_pars(alpha=c(0.025,0.024),fname_sero=sero_file)
write_res_pars(p_tau=P_TAU,cost_res=0.17,tau=0.3,fname_res=res_file)


## Run
v0_res <- run_simulations(res_file, sero_file, niter=NITER, label="Metabolic")
# first add strain_pars
strain_pars <- read.csv("strain_pars.csv")
strain_pars$id <- paste0("strain ", c(1,2,3))

v0_res_summ <- v0_res %>%
  reshape2::melt(id=c("time", "iter","label")) %>%
  rename(id=variable) %>%
  inner_join(strain_pars,by="id")%>%
  mutate(era = case_when(time >= t_vax-10  & time < t_vax ~ "pre-vax", 
                         time >= t_max-10 ~ "post-vax",
                          TRUE ~ "NA")) %>%
  filter(era != "NA") %>%
  group_by(era, id, id_m, id_s, id_r, id_v, iter,label) %>%
  summarize(freq=mean(value)) %>%
  group_by(era, iter,label) %>%
  nest() %>%
  summarize(map_dfr(data, function(df){
    df_new <- data.frame(
      "NVT" = sum(filter(df, id_v==0)$freq),
      "MT1" = sum(filter(df, id_m==1)$freq),
      "Resistance" = sum(filter(df, id_r==1)$freq))
    return(df_new)
  })) %>%
  reshape2::melt(id=c("era","iter","label")) %>%
  group_by(variable,era,label) %>%
  summarize(p50 = quantile(value, prob=0.5),
          p05 = quantile(value, prob=0.05),
          p95 = quantile(value, prob=0.95))



## Compare to no comp
setup_sim_dir(
  nstrain=3, 
    nhost=10000, 
    strain_id=c(1,2,3), 
    st=c(1,2,2), 
    rt=c(0,0,1), 
    mt=c(1,3,2),
    vt=c(1,0,0), 
    t_max=t_max,
    t_vax=t_vax,
    beta=BETA)

## Run
v0_res2 <- run_simulations(res_file, sero_file, niter=NITER, label="No Competition")

## Plot change in frequency pre-post vax
# first add strain_pars
strain_pars <- read.csv("strain_pars.csv")
strain_pars$id <- paste0("strain ", c(1,2,3))


## PUT TOGETHER FOR PLOTTING
v0_res <- rbind(v0_res, v0_res2)
v0_res$label <- factor(v0_res$label, levels=c("No Competition", "Metabolic"))

## PUT TOGETHER FOR PLOTTING
v0_res <- rbind(v0_res, v0_res2)
v0_res$label <- factor(v0_res$label,  levels=c("No Competition", "Metabolic"))

v0_res <- v0_res %>% filter(time > 5*365)
## Plot relative strain frequencies
gg_traj1 <- plot_trajectories(v0_res,plot_rel_freq=FALSE)+
  labs(y="Prevalence", x="Time (years)", col="")+
   scale_color_manual(values=c("#B6854D", "#02401B", "#A2A475"))+
   scale_fill_manual(values=c("#B6854D", "#02401B", "#A2A475"))+
  facet_wrap(~label, ncol=1,strip.position="right")+
  labs(y="Prevalence")+
  ggtitle("")

gg_resistance <- gg_traj1 +  theme(axis.title.y = element_blank())
### put together

gg_traj <- ggarrange(plotlist=list(
  ggarrange(plotlist=list(ggplot(data.frame(type=c("i.","ii.")))+facet_wrap(~type),
    ggplot()+theme_minimal()), widths=c(1,0.4)),
ggarrange(plotlist=list(gg_alpha+ggtitle("i.")+theme(plot.title = element_text(hjust = 0.5)), gg_resistance+ggtitle("ii.")+theme(plot.title = element_text(hjust = 0.5))), nrow=1, common.legend=TRUE, legend="bottom")),  nrow=2, labels=c("A","B"), font.label=c(size=16), heights=c(1,2.5))

