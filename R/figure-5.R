#############################
# CODE FOR GENERATING FIGURE 5
# Run 3 ATs, 3 MTs
# AT1 is VT
#############################
sc_file <- "sc_pars.csv"
sero_file <- "sero_pars.csv"
res_file <- "res_pars.csv" 

t_vax <- 15*365
t_max <- 30*365
# init = c(0,1,0,0,0,0,0,1)
setup_sim_dir(
  nstrain=12, 
  nhost=NHOST, 
  strain_id=c(1,2,3,4,5,6,7,8,9,10,11,12), 
  A=c(1,1,1,1,2,2,2,2,3,3,3,3), 
  G=c(1,2,3,4,1,2,3,4,1,2,3,4), 
  R=c(0,0,0,0,0,0,0,0,0,0,0,0),
  V=c(1,1,1,1,0,0,0,0,0,0,0,0), 
  init = c(1,0,0,1,0,1,0,0,0,0,1,0),
  t_max=t_max,
  t_vax=t_vax,
  beta=0.09,
  t_g=0.00005,
  t_l=0.00005)

write_sero_pars(alpha=c(0.02,0.025,0.03),fname=sero_file)
write_sc_pars(kappa=c(1.1,1.1,1.1,1.1),fname=sc_file)
write_res_pars(p_tau=0.2,cost_res=0.12,tau=0.8,fname_res=res_file)
if (load_data) {
  if (file.exists(paste0(data_dir, "/fig-5-data.RDS"))) {
    message("loading data for fig 5")
    model_res <- readRDS(paste0(data_dir, "/fig-5-data.RDS"))
    } else {
    message("could not load data for fig 5")
    model_res <- run_simulations(sc_file, sero_file, res_file, niter=NITER, label="sim")
    }
  } else {
  model_res <- run_simulations(sc_file, sero_file, res_file, niter=NITER, label="sim")
}

if (save_data){
  saveRDS(model_res, paste0(data_dir, "/fig-5-data.RDS"))
}

#####################
## view trajectories 
p1 <- plot_trajectories(model_res,plot_rel_freq=FALSE)

## summarize
strain_pars <- read.csv("strain_pars.csv")
strain_pars$id <- paste0("strain ", seq(1,nrow(strain_pars)))

dat <- model_res %>%
      reshape2::melt(id=c("time", "iter","label"), 
        variable.name="id", 
        value.name="freq") %>%
      inner_join(strain_pars,by="id") %>%
      mutate(era = case_when(time < t_vax ~ "PREPCV", TRUE ~ "POSTPCV")) %>%
      group_by(era) %>%
      filter(time == max(time))# %>%
      #group_by(label, id_m, id_s,id_r,era) %>%
      #summarize(freq=median(freq))

dat$id_r <- as.factor(dat$id_r)
dat$id_m <- as.factor(dat$id_m)
dat$id_s <- as.factor(dat$id_s)
dat$era <- factor(dat$era, levels=c("PREPCV", "POSTPCV"))

dat <- dat %>%
  mutate(sero = case_when(id_s == 1 ~ "VT", id_s==2 ~ "NVT1", id_s==3 ~ "NVT2"),
  id_m = paste0("MT", id_m))

g1 <- ggplot(dat, aes(x=id_m, y=freq, fill=era))+
  geom_boxplot(outlier.shape=NA)+
  facet_wrap(~sero)+
  labs(x="", y="Prevalence", fill="")+
  scale_fill_manual(values=freq_pal)


dat2 <- dat %>%
  filter(id_v==0) %>%
  group_by(time, iter, id_s,sero,era) %>%
  summarize(freq=sum(freq)) %>%
  group_by(time, iter) %>%
  mutate(freq=freq/sum(freq))
  

g2 <- ggplot(dat2, aes(y=freq, x= sero, fill=era))+
  geom_boxplot(outlier.shape=NA)+
  labs(x="", y="Relative freq among NVTs", fill="")+
  scale_fill_manual(values=freq_pal)+
  theme(legend.position="none")


dat3 <- dat %>%
  group_by(time, iter, id_s,sero,era) %>%
  summarize(freq=sum(freq)) %>%
  group_by(time, iter) %>%
  mutate(freq=freq/sum(freq)) %>%
  filter(sero != "VT")

g3 <- ggplot(dat3, aes(y=freq, x= sero, fill=era))+
  geom_boxplot(outlier.shape=NA)+
  labs(x="", y="Relative freq", fill="")+
  scale_fill_manual(values=freq_pal)+
  theme(legend.position="none")


ggarrange(plotlist=list(g1,ggarrange(plotlist=list(g3,g2), nrow=1)), nrow=2, labels=c("A","B"), font.label=c(size=20), common.legend=TRUE, widths=c(2,1), legend="bottom")

ggsave(paste0(out_dir, "/3sero-boxplots.pdf"), height=8, width=8)
