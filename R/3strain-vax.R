
res_file <- "res_pars.csv"
sero_file <- "sero_pars.csv"
t_vax <- 5*365 #don't run as long here
t_max <- 10*365



rel_ch_nvt <- function(dat, t_vax){
 dat <- dat %>%
    mutate(era=case_when(time < t_vax ~ "pre-vax", TRUE ~ "post-vax")) %>%
    group_by(era) %>%
    filter(time == max(time)) %>%
    mutate(value=`strain 2`+`strain 3`) %>%
    dplyr::select(value,iter, label, era) %>%
    group_by(era, iter, label, value) %>%
    summarize(value=mean(value)) %>% #??? doubled values why
    tidyr::pivot_wider(names_from="era", values_from="value") %>%
    mutate(delta = (`post-vax`-`pre-vax`)/`pre-vax`)     
  return(dat)
}


rel_ch_mt <- function(dat, t_vax){
 dat <- dat %>%
    mutate(era=case_when(time < t_vax ~ "pre-vax", TRUE ~ "post-vax")) %>%
    group_by(era) %>%
    filter(time == max(time)) %>%
    mutate(value=`strain 2`) %>%
    dplyr::select(value,iter, label, era) %>%
    group_by(era, iter, label, value) %>%
    summarize(value=mean(value)) %>% #??? doubled values why
    tidyr::pivot_wider(names_from="era", values_from="value") %>%
    mutate(delta = (`post-vax`-`pre-vax`)/`pre-vax`)   
  return(dat)
}

#######################
# system A
#######################

comp_type <- c("Metabolic", "None")

ratio <- seq(0.8, 1.2, by=0.1)
pars <- tidyr::crossing(
    alpha1=0.025, 
    alpha2=ratio*0.025,
    comp_type=comp_type, 
    rt1=0, 
    rt2=0,
    rt3=1, 
    p_tau=P_TAU, 
    cost=COST,
    tau=TAU,
    t_max=t_max, 
    t_vax=t_vax,
    vt1=1,
    vt2=0,
    vt3=0)

res <- pars %>%
    pmap(run_3strain)



deltas <- lapply(res, rel_ch_nvt, t_vax=t_vax)

res_df <- rbindlist(deltas, idcol = "index") %>%
  mutate(alpha2 = pars$alpha2[index]) %>%
  mutate(ratio = alpha2/pars$alpha1)

res_df$ratio <- factor(res_df$ratio)
res_df_st1 <- res_df
res_df_st1$system <- "A"
res_df_st1$type <- "Serotype"

deltas <- lapply(res, rel_ch_mt, t_vax=t_vax)

res_df <- rbindlist(deltas, idcol = "index") %>%
  mutate(alpha2 = pars$alpha2[index]) %>%
  mutate(ratio = alpha2/pars$alpha1)

res_df$ratio <- factor(res_df$ratio)
res_df_mt1 <- res_df
res_df_mt1$system <- "A"
res_df_mt1$type <- "Metabolic Type"



#######################
# system B
#######################

pars <- tidyr::crossing(
    alpha1=0.025, 
    alpha2=ratio*0.025,
    comp_type=comp_type, 
    rt1=0, 
    rt2=0,
    rt3=0, 
    p_tau=P_TAU, 
    cost=COST,
    tau=TAU,
    t_max=t_max, 
    t_vax=t_vax,
    vt1=1,
    vt2=0,
    vt3=0)

res <- pars %>%
    pmap(run_3strain)


deltas <- lapply(res, rel_ch_nvt, t_vax=t_vax)

res_df <- rbindlist(deltas, idcol = "index") %>%
  mutate(alpha2 = pars$alpha2[index]) %>%
  mutate(ratio = alpha2/pars$alpha1)

res_df$ratio <- factor(res_df$ratio)
res_df_st3 <- res_df
res_df_st3$system <- "C"
res_df_st3$type <- "Serotype"



deltas <- lapply(res, rel_ch_mt, t_vax=t_vax)

res_df <- rbindlist(deltas, idcol = "index") %>%
  mutate(alpha2 = pars$alpha2[index]) %>%
  mutate(ratio = alpha2/pars$alpha1)

res_df$ratio <- factor(res_df$ratio)
res_df_mt3 <- res_df
res_df_mt3$system <- "C"
res_df_mt3$type <- "Metabolic Type"

#######################################
# put together
res_df <- rbind(res_df_st1,res_df_mt1,
    res_df_mt3,res_df_st3)

res_df$ratio <- factor(res_df$ratio)
g1 <- ggplot(filter(res_df, !(is.na(delta))), aes(x=ratio, y=`pre-vax`, fill=label))+
  geom_boxplot(notch=FALSE, outlier.shape=NA)+
  facet_grid(type~system, scales="free_y")+
  labs(x=expression(alpha[NVT]/alpha[VT]), y="Pre-vaccination prevalence", fill="Competition")+
  scale_fill_manual(values=c("#C93312", "#899DA4"))

g2 <- ggplot(filter(res_df, !(is.na(delta))), aes(x=ratio, y=`post-vax`, fill=label))+
  geom_boxplot(notch=FALSE, outlier.shape=NA)+
  facet_grid(type~system, scales="free_y")+
  labs(x=expression(alpha[NVT]/alpha[VT]), y="Post-vaccination prevalence", fill="Competition")+
  scale_fill_manual(values=c("#C93312", "#899DA4"))


gg_res <- ggplot(filter(res_df, !(is.na(delta)), `pre-vax` > 0.01), aes(x=ratio, y=delta, fill=label))+
  geom_boxplot(notch=FALSE, outlier.shape=NA)+
  facet_grid(type~system, scales="free_y")+
  labs(x=expression(alpha[NVT]/alpha[VT]), y="Relative Change", fill="Competition")+
  scale_fill_manual(values=c("#C93312", "#899DA4"))



gg_3strainvax <- ggarrange(plotlist=list(
  ggarrange(plotlist=list(
   ggplot()+theme_minimal(),
    ggplot()+theme_minimal(), 
    ggplot()+theme_minimal()),nrow=1), 
  ggarrange(plotlist=list(
    gg_res+theme(strip.text.x = element_blank())), 
    nrow=1)),nrow=2,
    font.label=list(size=16), heights=c(1,1.75,2))


gg_freqs <- ggarrange(plotlist=list(
  ggarrange(plotlist=list(
   ggplot()+theme_minimal(),
    ggplot()+theme_minimal(), 
    ggplot()+theme_minimal()),nrow=1), 
  ggarrange(plotlist=list(
    g1+theme(strip.text.x = element_blank()), 
    g2+theme(strip.text.x = element_blank())), 
    nrow=2, common.legend=TRUE)),nrow=2,
    font.label=list(size=16), heights=c(1,3,3))


