source("setup.R")
t_max <- 50*365
t_vax <- 25*365


sc_file <- "sc_pars.csv"
sero_file <- "sero_pars.csv"
res_file <- "res_pars.csv" 

run_over_pars <- function(a, p_tau){
  write_sero_pars(alpha=c(0.02, a),fname=sero_file)
  write_sc_pars(kappa=c(1.1, 1.1),fname=sc_file)
  write_res_pars(p_tau=p_tau,cost_res=0.2,tau=0.9,fname_res=res_file)
  model_res1 <- run_simulations(sc_file, sero_file, res_file, niter=no_cores*3, label=paste(p_tau))
  return(model_res1)
}

# init = c(0,1,0,0,0,0,0,1)
setup_sim_dir(
  nstrain=8, 
  nhost=5000, 
  strain_id=c(1,2,3,4,5,6,7,8), 
  A=c(1,1,1,1,2,2,2,2), 
  G=c(1,1,2,2,1,1,2,2), 
  R=c(0,1,0,1,0,1,0,1),
  V=c(1,1,1,1,0,0,0,0), 
  init = c(0,1,0,0,0,0,1,0),
  t_max=t_max,
  t_vax=t_vax,
  beta=0.09,
  t_g=0.00005,
  t_l=0.00005)
 

## Run
params <- expand.grid(a=c(0.023), p_tau = c(0.33, 0.34, 0.35, 0.36))

model_res <- params %>%
  pmap_dfr(run_over_pars)


###################################
## Plot relative strain frequencies
p1 <- plot_trajectories(model_res,plot_rel_freq=FALSE)+
  scale_color_manual(
    labels=c('AT1MT1S', 'AT1MT1R', 'AT1MT2S', 'AT1MT2R', 'AT2MT1S', 'AT2MT1R', 'AT2MT2S', 'AT2MT2R'),  
    values=rainbow(8))+
    facet_wrap(~label, nrow=4)

strain_pars <- read.csv("strain_pars.csv")
strain_pars$id <- paste0("strain ", seq(1,8))

## pre & post vax relative freq of resistance types and meta types
# look at cases where replacement did happen

res <- model_res %>%
  mutate(era=case_when(time <= t_vax ~ "pre-vax", TRUE ~ "post-vax")) %>%
  group_by(era) %>%
  filter(time==max(time)) %>%
  group_by(era, iter, label) %>%
  mutate(tot_freq = sum(across(starts_with("strain")), na.rm = T)) %>%
  summarize(MT1 = (`strain 1` + `strain 2` + `strain 5` + `strain 6`)/tot_freq,
      R = (`strain 2` + `strain 4` + `strain 6` + `strain 8`)/tot_freq) %>%
  group_by(iter, label) %>% nest() %>%
  mutate(did_switch = map_dbl(data, function(df) {
      x <- filter(df, era=="post-vax")$MT1
      if (x > 0) {return(1)}
      else {return(0)}
    })) %>%
  filter(did_switch > 0) %>%
  unnest(cols=c(data)) %>%
  reshape2::melt(id=c("era", "iter", "label", "did_switch")) 

res$era <- factor(res$era, levels=c("pre-vax", "post-vax"))

ggplot(res, aes(x=label, y=value, fill=era))+
  geom_boxplot(outlier.shape=NA)+
  facet_wrap(~variable)

