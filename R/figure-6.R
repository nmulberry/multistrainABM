#############################
# CODE FOR GENERATING FIGURE 6
# Run 2 ATs, 2 MTs, with resistance
#############################

t_max <- 20*365
t_vax <- 10*365


sc_file <- "sc_pars.csv"
sero_file <- "sero_pars.csv"
res_file <- "res_pars.csv" 
#0.2, 0.9, 0.34
run_over_pars <- function(a, p_tau){
  write_sero_pars(alpha=c(0.025, a),fname=sero_file)
  write_sc_pars(kappa=c(1.1, 1.1),fname=sc_file)
  write_res_pars(p_tau=p_tau,cost_res=0.2,tau=0.9,fname_res=res_file)
  model_res1 <- run_simulations(sc_file, sero_file, res_file, niter=NITER, label=paste(a))
  return(model_res1)
}

# init = c(0,1,0,0,0,0,0,1)
setup_sim_dir(
  nstrain=8, 
  nhost=NHOST, 
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
params <- expand.grid(a=c(0.026, 0.027, 0.028), p_tau = c(0.34))


if (load_data) {
  if (file.exists(paste0(data_dir, "/fig-6-data.RDS"))) {
    message("loading data for fig 6")
    model_res <- readRDS(paste0(data_dir, "/fig-6-data.RDS"))
    } else {
    message("could not load data for fig 6")
    model_res <- params %>%
      pmap_dfr(run_over_pars)
    }
  } else {
  model_res <- params %>%
    pmap_dfr(run_over_pars)
}

if (save_data){
  saveRDS(model_res, paste0(data_dir, "/fig-6-data.RDS"))
}


## get sims where there is escape
model_res2 <- model_res %>%
  group_by(iter,label) %>%
  nest() %>%
  mutate(keep=map_chr(data, function(df) {
    last_time <- filter(df, time==max(time))
    if (last_time$`strain 5` +last_time$`strain 6`  > 0) {
      return("Y") } else {return("N")}})) %>%
  filter(keep=="Y") %>%
  unnest(cols=c(data))
###################################
## Plot relative strain frequencies
p1 <- plot_trajectories(model_res2,plot_rel_freq=FALSE)+
  scale_color_manual(
    labels=c('AT1-MT1-S', 'AT1-MT1-R', 'AT1-MT2-S', 'AT1-MT2-R', 'AT2-MT1-S', 'AT2-MT1-R', 'AT2-MT2-S', 'AT2-MT2-R'),  
    values=rainbow(8))+
   # facet_grid(label~., labeller=label_bquote(alpha[2]:~.(label)))+
    facet_wrap(~label,nrow=1)+
theme(
  strip.background = element_blank(),
  strip.text.x = element_blank())+
    labs(x="Time (yrs)", y="Prevalence", col="Strain")
    
ggsave(paste0(out_dir, "/resistance-figure-main.pdf"), height=5, width=8)
