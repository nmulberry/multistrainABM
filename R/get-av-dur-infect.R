

run_over_alpha <- function(a, p=0, r=0){
  dur_infect <- run_one_host(
    nstrain=1,
    strain_id=c(1),
    st=c(1), 
    rt=r, 
    mt=c(1),
    vt=c(0), 
    alpha=c(a), 
    p_tau=p, 
    cost_res=COST, 
    tau=TAU, 
    t_max=3000,
    t_vax=3000,
    t_inv=0,
    init_strain=c(1)
  ) %>%
    filter(time > 5, X < 0.01) %>%
    summarize(t = min(time))

  dur_immune <- run_one_host(
    nstrain=1,
    strain_id=c(1),
    st=c(1), 
    rt=r, 
    mt=c(1),
    vt=c(0), 
    alpha=c(a), 
    p_tau=p, 
    cost_res=COST, 
    tau=TAU, 
    t_max=3000,
    t_vax=3000,
    t_inv=0,
    init_strain=c(1)
  ) %>%
    filter(time > 5, I < 1e-3) %>%
    summarize(t = min(time))

  return(data.frame(dur_infect=dur_infect$t, dur_immune=dur_immune$t))
}

alpha_vec <- seq(0.02, 0.08, by=0.001)


res1 <- alpha_vec %>% map_dfr(run_over_alpha, p=0, r=0)
res1$alpha <- alpha_vec
res1$type <- "Sensitive"
res1$treat <- "Untreated"

res2 <- alpha_vec %>% map_dfr(run_over_alpha, p=1, r=0)
res2$alpha <- alpha_vec
res2$type <- "Sensitive"
res2$treat <- "Treated"

res3 <- alpha_vec %>% map_dfr(run_over_alpha, p=0, r=1)
res3$alpha <- alpha_vec
res3$type <- "Resistant"
res3$treat <- "Untreated"

res <- rbind(res1,res2,res3)

g1 <- ggplot(res, aes(x=alpha, y=dur_infect/30, col=type, linetype=treat))+
  geom_line(size=1.5)+
  labs(y="Duration of Infection (months)", x=expression(alpha), col="", linetype="")+
  scale_color_manual(values=c("Resistant"="#E69F00","Sensitive"="#56B4E9"))


g2 <- ggplot(res, aes(x=alpha, y=dur_immune/365, col=type, linetype=treat))+
  geom_line(size=1.5)+
  labs(y="Duration of Immunity (yrs)", x=expression(alpha), col="", linetype="")+
  scale_color_manual(values=c("Resistant"="#E69F00","Sensitive"="#56B4E9"))

ggarrange(plotlist=list(g1,g2), labels=c("A","B"), font.label=list(size=16),
  common.legend=TRUE, legend="bottom")



