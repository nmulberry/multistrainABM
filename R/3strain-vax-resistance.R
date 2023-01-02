
res_file <- "res_pars.csv"
sero_file <- "sero_pars.csv"
t_vax <- 10*365 #don't run as long here
t_max <- 20*365

three_pal <- c("strain 1"="#273046", "strain 2"="#00A08A","strain 3"="#FAD510")

#### system inspired by 11A and 9V

setup_sim_dir(
    nstrain=3, 
    nhost=10000, 
    strain_id=c(1,2,3), 
    st=c(1,2,2), 
    rt=c(1,1,0), 
    mt=c(1,1,2),
    vt=c(1,0,0), 
    t_max=t_max,
    t_vax=t_vax,
    beta=BETA)

write_sero_pars(alpha=c(0.035,0.027),fname_sero=sero_file)
write_res_pars(p_tau=0.4,cost_res=0.16,tau=0.35,fname_res=res_file)


## Run
v0_res <- run_simulations(res_file, sero_file, niter=NITER, label="sim")

p1 <- plot_trajectories(filter(v0_res, time > 2*365), plot_rel_freq=FALSE)+
  scale_color_manual(values=three_pal)+
  scale_fill_manual(values=three_pal)+
  labs(x="Time (yrs)", y="Prevalence", col="")+
  theme(strip.text.x=element_blank())
## rel freqs
strain_pars <- read.csv("strain_pars.csv")
strain_pars$id <- paste0("strain ", c(1,2,3))

strain_pars <- strain_pars %>%
  mutate(ST=case_when(id_s==1 ~ "9V", TRUE ~ "11A"), RT=id_r)

sim_data <- v0_res %>%
  reshape2::melt(id=c("time", "iter","label")) %>%
  rename(id=variable) %>%
  inner_join(strain_pars,by="id")%>%
  mutate(era = case_when(time < t_vax ~ "pre-vax", TRUE ~ "post-vax")) %>%
  group_by(era) %>%
  filter(time==max(time)) %>%
  group_by(iter, era) %>%
  mutate(freq = value/sum(value)) %>%
  group_by(ST, RT, era) %>%
  summarize(Rel_Freq = quantile(freq, prob=0.5), 
    Rel_Freq_95 = quantile(freq, prob=0.95),
    Rel_Freq_5 = quantile(freq, prob=0.05)) %>%
  mutate(label="Simulated") %>%
  ungroup() 



## DATA
freqs_raw <- readRDS("../data/strain_freqs_simplified.RDS") %>% 
  filter(ST %in% c("9V", "11A"))
# renormalize

freqs <- freqs_raw %>%
  group_by(Time) %>%
  mutate(Rel_Freq=Freq/sum(Freq)) %>%
  group_by(Time, ST, RT) %>%
  summarize(Rel_Freq=sum(Rel_Freq)) %>%
  ungroup() %>%
  filter(Time != "Y3") %>%
  mutate(label="Data") %>%
  mutate(era = case_when(Time=="Y0" ~ "pre-vax", TRUE ~ "post-vax")) %>%
  dplyr::select(-c(Time)) %>%
  mutate(Rel_Freq_95=NA, Rel_Freq_5=NA)


fake_dat <- data.frame(ST=c("9V","9V"), RT=c(0,0), era=c("post-vax", "pre-vax"),
  Rel_Freq=c(0,0), Rel_Freq_95=c(0,0), Rel_Freq_5=c(0,0), label="Simulated")
all_freqs <- rbind(sim_data, freqs, fake_dat)
all_freqs <- all_freqs %>% mutate(RT=case_when(RT==0~ "Sensitive", TRUE ~ "Resistant"))
all_freqs$era <- factor(all_freqs$era, levels=c('pre-vax', 'post-vax'))
# plot
gg <- ggplot(all_freqs, aes(x=era, y=Rel_Freq, fill=label))+
  geom_col(col="black", position="dodge")+
  facet_grid(RT~ST, scales="free")+
  geom_errorbar(aes(x=era, ymin=Rel_Freq_5, ymax=Rel_Freq_95, group=label), size=.5,
                      width=.2,
                      position=position_dodge(0.9))+
  labs(x="", y="Relative Frequency", fill="")+
  scale_fill_brewer(palette="Blues")


ggarrange(plotlist=list(
  ggarrange(plotlist=list(ggplot()+theme_minimal(), gg),labels=c("A","B", "C"), font.label=list(size=16)), 
  p1), nrow=2, labels=c("", "C"), font.label=list(size=16))
