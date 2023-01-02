
res_file <- "res_pars.csv"
sero_file <- "sero_pars.csv"
t_vax <- 5*365 #don't run as long here
t_max <- 10*365

#######################
# Both sensitive
#######################

## plots
make_plots <- function(dat){
    dat <- dat %>% dplyr::select(-c(`strain 1`)) %>%
      rename(y=`strain 2`) %>%
      filter(time > 2*365)
  gg <- dat %>%
      group_by(time,label) %>%
      summarize(p50 = quantile(y, prob=0.5),
          p05 = quantile(y, prob=0.05),
          p95 = quantile(y, prob=0.95)) %>%
      ggplot(aes(x=time/365))+
        geom_line(aes(y=p50, col=label), size=1.5)+
        geom_ribbon(aes(ymin=p05, ymax=p95, fill=label), alpha=0.3)+
        guides(fill="none")+
        scale_x_continuous(breaks = c(t_vax/365), labels=expression(t[vax]))+
        labs(y="NVT Prevalence", x="", col="")+
        scale_color_manual(values=c("#C93312", "#899DA4"))+
        scale_fill_manual(values=c("#C93312", "#899DA4")) 
      return(gg)
}



comp_type <- c("Metabolic", "None")
alpha_ratio <- seq(1.2, 2.4, by=0.2)

pars <- tidyr::crossing(alpha1=0.02, 
    alpha2=0.02*alpha_ratio,comp_type=comp_type,
    rt1=0, rt2=1,
    p_tau=P_TAU, cost=COST, tau=TAU,
    t_max=t_max, t_vax=t_vax,
    vt1=1, vt2=0)

res <- pars %>%
    pmap(run_2strain)

rel_ch <- function(dat, t_vax){
  dat <- dat %>%
    mutate(era=case_when(time < t_vax ~ "pre-vax", TRUE ~ "post-vax")) %>%
    group_by(era) %>%
    filter(time == max(time)) %>%
    ungroup() %>%
    dplyr::select(value=`strain 2`, iter, label, era) %>%
    tidyr::pivot_wider(names_from="era", values_from="value") %>%
    mutate(delta = `post-vax`/`pre-vax`,
        delta1 = `post-vax`-`pre-vax`,
        delta2 = (`post-vax`-`pre-vax`)/`pre-vax`)
  return(dat)
}

deltas <- lapply(res, rel_ch, t_vax=t_vax)

res_df1 <- rbindlist(deltas, idcol = "index") %>%
  mutate(alpha2 = pars$alpha2[index]) %>%
  mutate(ratio = alpha2/pars$alpha1)

res_df1$ratio <- factor(res_df1$ratio)
res_df1$system <- "A"
#####################


pars <- tidyr::crossing(alpha1=0.02, 
    alpha2=0.02*alpha_ratio,comp_type=comp_type,
    rt1=0, rt2=0,
    p_tau=P_TAU, cost=COST, tau=TAU,
    t_max=t_max, t_vax=t_vax,
    vt1=1, vt2=0)

res <- pars %>%
    pmap(run_2strain)

deltas <- lapply(res, rel_ch, t_vax=t_vax)

res_df2 <- rbindlist(deltas, idcol = "index") %>%
  mutate(alpha2 = pars$alpha2[index]) %>%
  mutate(ratio = alpha2/pars$alpha1)

res_df2$ratio <- factor(res_df2$ratio)

res_df2$system <- "B"



res_df <- rbind(res_df1, res_df2)
res_df <- filter(res_df, `pre-vax` > 0.0001)
 g1 <- ggplot(res_df, aes(x=ratio, y=`pre-vax`, fill=label))+
  geom_boxplot(notch=FALSE, outlier.shape=NA)+
  labs(x=expression(alpha[NVT]/alpha[VT]), y="Pre-vaccination prevalence", fill="Competition")+
  scale_fill_manual(values=c("#C93312", "#899DA4"))+
  facet_wrap(~system)+
  facet_wrap(~system)+theme( strip.text.x = element_blank())

 g2 <- ggplot(res_df, aes(x=ratio, y=`post-vax`, fill=label))+
  geom_boxplot(notch=FALSE, outlier.shape=NA)+
  labs(x=expression(alpha[NVT]/alpha[VT]), y="Post-vaccination prevalence", fill="Competition")+
  scale_fill_manual(values=c("#C93312", "#899DA4"))+
  facet_wrap(~system)+
  facet_wrap(~system)+theme( strip.text.x = element_blank())

g3 <- ggplot(filter(res_df, !(is.na(delta))), aes(x=ratio, y=delta, fill=label))+
  geom_boxplot(notch=FALSE,outlier.shape=NA)+
  ylim(c(0,3))+
  labs(x=expression(alpha[NVT]/alpha[VT]), y="Ratio", fill="Competition")+
  scale_fill_manual(values=c("#C93312", "#899DA4"))+
  facet_wrap(~system)+theme( strip.text.x = element_blank())



#################################
# ALSO RUN EXAMPLE TRAJECTORIES
##################################

alpha_ratio <- 1.2
alpha1 <- 0.02
pars <- tidyr::crossing(comp_type=comp_type,
    alpha1=alpha1, 
    alpha2=alpha1*alpha_ratio,rt1=0, rt2=0, p_tau=P_TAU, cost=COST, tau=TAU,
    t_max=t_max, t_vax=t_vax,
    vt1=1, vt2=0)

res0 <- pars %>%
    pmap(run_2strain)
gg_ex0 <- make_plots(rbindlist(res0))+
  ggtitle(as.expression(bquote(alpha[NVT]~"="~.(alpha_ratio)~alpha[VT])))



alpha_ratio <- 1.8
alpha1 <- 0.02
pars <- tidyr::crossing(comp_type=comp_type,
    alpha1=alpha1, 
    alpha2=alpha1*alpha_ratio,rt1=0, rt2=0, p_tau=P_TAU, cost=COST, tau=TAU,
    t_max=t_max, t_vax=t_vax,
    vt1=1, vt2=0)

res1 <- pars %>%
    pmap(run_2strain)
gg_ex1 <- make_plots(rbindlist(res1))+
  ggtitle(as.expression(bquote(alpha[NVT]~"="~.(alpha_ratio)~alpha[VT])))

#######
alpha_ratio <- 2.0
alpha1 <- 0.02
pars <- tidyr::crossing(comp_type=comp_type,
    alpha1=alpha1, 
    alpha2=alpha1*alpha_ratio)
res2 <- pars %>%
    pmap(run_2strain,rt1=0, rt2=0, p_tau=P_TAU, cost=COST, tau=TAU,
    t_max=t_max, t_vax=t_vax,
    vt1=1, vt2=0)
gg_ex2 <- make_plots(rbindlist(res2))+
  ggtitle(as.expression(bquote(alpha[NVT]~"="~.(alpha_ratio)~alpha[VT])))



gg_2strainvax <- ggarrange(plotlist=list(
  ggarrange(plotlist=list(
   ggplot()+theme_minimal(),
    ggplot()+theme_minimal()), nrow=1), 
  g3),nrow=2,  
    font.label=list(size=16), heights=c(1,3))


gg_freqs <- ggarrange(plotlist=list(
  ggarrange(plotlist=list(
   ggplot()+theme_minimal(),
    ggplot()+theme_minimal()), nrow=1), 
  ggarrange(plotlist=list(g1,g2), nrow=2, common.legend=TRUE, legend="bottom")),nrow=2,  
    font.label=list(size=16), heights=c(1,3))

gg_traj <- ggarrange(plotlist=list(
    gg_ex0,
    gg_ex1,
    gg_ex2+labs(x="Time (yrs)")), 
    nrow=3, common.legend=TRUE, legend="right")


