#############################
# SIMPLE EXAMPLE ON USING MODEL
#############################
source("setup.R")

sc_file <- "sc_pars.csv"
sero_file <- "sero_pars.csv"
res_file <- "res_pars.csv" 

NITER <- 20


run_over_pars <- function(a, k){
  write_sero_pars(alpha=c(0.02, a),fname=sero_file)
  write_sc_pars(kappa=c(k, 1.1),fname=sc_file)
  write_res_pars(p_tau=0.22,cost_res=0.11,tau=0.8,fname_res=res_file)
  model_res1 <- run_simulations(sc_file, sero_file, res_file, niter=NITER, label=paste0(a,",",k))
  return(model_res1)
}

t_max <- 30*365
t_vax <- 15*365
# init = c(0,1,0,0,0,0,0,1)
setup_sim_dir(
  nstrain=4, 
  nhost=5000, 
  strain_id=c(1,2,3,4), 
  A=c(1,1,2,2), 
  G=c(1,2,1,2), 
  R=c(0,0,0,0),
  V=c(1,1,0,0), 
  init = c(1,0,0,1),
  t_max=t_max,
  t_vax=t_vax,
  beta=0.09,
  t_g=0.00005,
  t_l=0.00005)

# param sweep
alpha2 <- c(0.022, 0.025, 0.03)
kappa1 <- c(1.1, 1.15, 1.2)

params <- expand.grid(a=alpha2, k=kappa1)

model_res <- params %>%
  pmap_dfr(run_over_pars)

p1 <- plot_trajectories(model_res,plot_rel_freq=FALSE)+
  labs(x="Time (yrs)", y="Prevalence", col="")+
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank(),
 # panel.border = element_rect(colour = "black", fill=NA, size=3),
  plot.margin = margin(2,1,1,2, "cm")
  )+
  scale_color_manual(values=wes_palette("Moonrise2"), labels=c('AT1MT1', 'AT1MT2', 'AT2MT1', 'AT2MT2'))+
  scale_fill_manual(values=wes_palette("Moonrise2"))


##summarize
strain_pars <- read.csv("strain_pars.csv")
strain_pars$id <- paste0("strain ", seq(1,nrow(strain_pars)))

dat <- model_res %>%
      reshape2::melt(id=c("time", "iter","label"), 
        variable.name="id", 
        value.name="freq") %>%
      inner_join(strain_pars,by="id") %>%
      filter(time < t_vax) %>%
      filter(time == max(time)) %>%
      group_by(label, id_m, id_s,id_r) %>%
      summarize(freq=median(freq)) %>%
      mutate(id_r= case_when(id_r==0 ~ "S", TRUE ~"R"),
              id_s = paste0("AT", id_s),
              id_m = paste0("MT", id_m)) 

dat$id_r <- as.factor(dat$id_r)
dat$id_m <- as.factor(dat$id_m)
dat$id_s <- as.factor(dat$id_s)

## pre vax structure
g1 <- ggplot(dat, aes(x=id_m, y=id_s, fill=freq, col=freq))+
  geom_tile()+
  facet_wrap(~label)+  
  scale_y_discrete(expand=c(0,0))+
  scale_x_discrete(expand=c(0,0))+
  guides(fill=guide_colourbar(barwidth=12, ticks=FALSE,title.position = "top"), color="none")+
  labs(x="", y="", fill="Pre-vax median freq ")+
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )+scale_color_viridis()+scale_fill_viridis()


## post-vax emergence
# frac sims with AT2MT1
dat2 <- model_res %>%
      reshape2::melt(id=c("time", "iter","label"), 
        variable.name="id", 
        value.name="freq") %>%
      inner_join(strain_pars,by="id") %>%
      filter(time == max(time)) %>%
      filter(id=='strain 3' & freq > 0) %>%
      group_by(label) %>%
      summarize(p_emerg = n()/NITER) %>%
      ungroup() %>%
      tidyr::separate(col=label, into=c("a2", "kappa1"), sep=",")

dat2$kappa1 <- as.factor(dat2$kappa1)
dat2$a2 <- factor(dat2$a2, levels=c(0.03, 0.025, 0.022))
g2 <- ggplot(dat2, aes(x=kappa1, y=a2, fill=p_emerg, col=p_emerg))+
  geom_tile()+ 
  scale_y_discrete(expand=c(0,0))+
  scale_x_discrete(expand=c(0,0))+
  guides(fill=guide_colourbar(barwidth=12, ticks=FALSE,title.position = "top"), color="none")+
  labs(y=expression(alpha[2]), x=expression(kappa[1]), fill="Post-vax emergence of AT2MT2 ")+
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )+scale_color_viridis()+scale_fill_viridis()

gg <- ggarrange(plotlist=list(g1,g2),labels=c("B","C"),  font.label=c(size=20))
ggarrange(plotlist=list(p1,gg), nrow=2, labels=c("A",""),  font.label=c(size=20))
#ggsave("run_2sero_2meta_both_init1001.pdf", height=10, width=8)






