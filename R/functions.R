#=====================================#
# SIMULATION FUNCTIONS                #
#=====================================#
# set up const par files
setup_sim_dir <- function(nstrain, # number of strains
                          nhost, # number of hosts
                          strain_id, # nstrain x1 array of strain ids 
                          st, # nstrain x 1 array of serotypes (integer-values starting at 1)
                          rt, # nstrain x 1 array of resistant types (binary-values)
                          mt,  # nstrain x 1 array of metabolic types (integer-values)
                          vt, # nstrain x 1 array of vaccine types (*BY SERO*)
                          gt, #nstrain x1 array of sero groups (*BY SERO*)
                          t_max, # max simulation time
                          t_vax, # vaccination time
                          beta=0.5, # bwh transmission rate
                          kappa=1.1, # growth rate of sensitive strain
                          eps_x=0.03, # wh death rate
                          eps_a=0.02, # immune decay rate
                          rho=0.01, # minimum wh pathogen density
                          mu=0.00, # host migration rate
                          k_g=0.1, # within sero-group cross-immunity
                          theta=70
  ){
  #----- check input-----#
  nsero <- length(unique(st))

  if (!(length(st)==nstrain)){
    stop("Check dimensions of st")
  }
  if (!(length(mt)==nstrain)){
    stop("Check dimensions of mt")
  }
  if (!(length(rt)==nstrain)){
    stop("Check dimensions of rt")
  } 
  if (!(length(vt)==nstrain)){
    stop("Check dimensions of vt")
  }
  # ------SETUP-------#
  write_constant_pars(nstrain=nstrain,nhost=nhost,t_max=t_max,t_vax=t_vax,
    kappa=kappa,eps_x=eps_x, eps_a=eps_a, rho=rho,theta=theta,mu=mu,k_g=k_g,beta=beta)
  write_strain_pars(nstrain=nstrain,strain_id=strain_id,st=st,mt=mt,rt=rt,vt=vt, gt=gt)
}

## Run for a single host
# to illustrate wh dynamics
run_one_host <- function(nstrain, 
                          strain_id, # nstrain x1 array of strain ids 
                          st, # nstrain x 1 array of serotypes (integer-values)
                          rt, # nstrain x 1 array of resistant types (binary-values)
                          mt,  # nstrain x 1 array of metabolic types (integer-values)
                          vt, # nstrain x 1 array of vaccine types (*BY SERO*)
                          alpha, # nsero x 1 array of alpha values (doubles)
                          p_tau, # f64 in (0,1), treatment coverage
                          cost_res, # f64 in (0,1), wh growth cost of resistance
                          tau, # f64 >=0
                          t_max, # max simulation time
                          t_vax, # vaccination time
                          t_inv, # invasion time of second strains (if applicable)
                          init_strain, # initial strain(s) by index
                          kappa=1.1, # growth rate of sensitive strain
                          eps_x=0.03, # wh death rate
                          eps_a=0.02, # immune decay rate
                          rho=0.01, # minimum wh pathogen density
                          theta=70
  ){ 

  if (t_inv > 0) {
      X0 <- rep(0, nstrain) 
      I0 <- rep(0, nstrain)

      X0[init_strain] <- rho
      y0 <- c(X=X0, I=I0)

      parameters <- c(nstrain=nstrain,eps_a=eps_a,
        eps_x=eps_x,rho=rho, theta=theta,
        kappa=kappa, alpha=alpha, 
        tau=tau, st=st, mt=mt,rt=rt,
        p_tau=p_tau,cost_res=cost_res)

      times <- seq(0, t_inv, by = 1)
      out1 <- ode(y=y0, times=times, func=withinhost_model, 
        parms=parameters)  

      times <- seq(t_inv, t_max, by=1)
      y0 <- tail(out1, n=1)[2:ncol(out1)]

      remaining_strain <- setdiff(1:nstrain, init_strain)
      y0[remaining_strain] <- rho

      out2 <- ode(y=y0, times=times, func=withinhost_model, 
        parms=parameters)

      out2 <- out2[2:(t_max-t_inv),]
      out <- rbind(out1, out2)
    } else {
      X0 <- rep(0, nstrain) 
      I0 <- rep(0, nstrain)

      X0[init_strain] <- rho
      y0 <- c(X=X0, I=I0)

      parameters <- c(nstrain=nstrain,eps_a=eps_a,
        eps_x=eps_x,rho=rho, theta=theta,
        kappa=kappa, alpha=alpha, 
        tau=tau, st=st, mt=mt,rt=rt,
        p_tau=p_tau,cost_res=cost_res)

      times <- seq(0, t_max, by = 1)
      out <- ode(y=y0, times=times, func=withinhost_model, 
        parms=parameters)
    }
   return(as.data.frame(out)) 
}

## Run over multiple hosts
# main model simulation function
# ASSUME: being called from working directoy
# with access to const_pars and strain_pars
run_simulations <- function(res_file,sero_file,niter,label="None"){
  if (!(file.exists(res_file))){
    "resistance parameter file not found"
    return(NULL)
  } else if (!(file.exists(sero_file))){
    "sero parameter file not found"
    return(NULL)
  } else {
    sim <- function(iter) {
             tmp <- paste0("tmp_out_", iter, ".csv")
             system(paste(multiabm, res_file, sero_file, ">", tmp), wait=TRUE)
             res <- read.csv(tmp, header=FALSE)
             names(res) <- c("time", paste("strain", seq(1:(length(names(res))-1))))
             res$iter <- iter
             file.remove(tmp)
             return(res) 
      }
    res <- future_map_dfr(
       1:niter, 
       ~ sim(.x),   
       .progress=TRUE) 
    res$label <- label
    return(res)
  }
}


run_trajectories <- function(res_file,sero_file){
  if (!(file.exists(res_file))){
    "resistance parameter file not found"
    return(NULL)
  } else if (!(file.exists(sero_file))){
    "sero parameter file not found"
    return(NULL)
  } else {

   tmp <- paste0("tmp_out.csv")
   system(paste(multiabm_trajectories, res_file, sero_file, ">", tmp), wait=TRUE)
   res <- read.csv(tmp, header=FALSE)
   names(res) <- c("time", paste("strain", seq(1:(length(names(res))-1))))
   file.remove(tmp)

    return(res)
  }
}
###########################
# 2-strain parameter sweep
###########################
run_2strain <- function(
  alpha1, 
  alpha2, 
  comp_type, 
  rt1, 
  rt2, 
  p_tau, 
  cost,
  tau,
  t_max, 
  t_vax,
  vt1=0,
  vt2=0,
  gt1=1,
  gt2=1
){

  if (comp_type == "Antigenic" | comp_type == "Both"){
    st <- c(1,1)
    alpha <- c(alpha1)
  } else {
    alpha <- c(alpha1, alpha2)
    st <- c(1,2)
  }

  if (comp_type == "Metabolic" | comp_type == "Both"){
    mt <- c(1,1)
  } else {
    mt <- c(1,2)
  }

  setup_sim_dir(
      nstrain=2, 
      nhost=NHOST, 
      strain_id=seq(1,2), 
      st=st, 
      rt=c(rt1,rt2), 
      mt=mt,
      vt=c(vt1,vt2), 
      gt=c(gt1,gt2),
      t_max=t_max,
      t_vax=t_vax,
      beta=BETA,
      mu=MU,
      k_g=CROSS)
  
  write_res_pars(p_tau=p_tau,cost_res=cost,tau=tau,fname_res=res_file)
  write_sero_pars(alpha=alpha,fname_sero=sero_file)
  ## Run
  v0_res <- run_simulations(res_file, sero_file, niter=NITER, label=comp_type)
  return(v0_res)
}

get_rel_freq_2strain <- function(alpha1, 
  alpha2, 
  comp_type, 
  rt1, 
  rt2, 
  p_tau, 
  cost,
  tau,
  t_max, 
  t_vax){  

  v0_res <- run_2strain(alpha1, alpha2, comp_type,
    rt1,rt2,p_tau, cost, tau, t_max, t_vax)

  v0_res <- filter(v0_res, time == max(time)) %>%
    rowwise() %>%
    mutate(tot=sum(c_across(starts_with("strain")))) %>%
    mutate(across(starts_with("strain"), ~ .x/tot)) %>%
    ungroup() %>%
    summarize(across(starts_with("strain"), ~ mean(.x)))    

  return(v0_res)    
}

##################################
# 3-strain, 2 sero parameter sweep
##################################
run_3strain <- function(
  alpha1, 
  alpha2, 
  comp_type, 
  rt1, 
  rt2,
  rt3, 
  p_tau, 
  cost,
  tau,
  t_max, 
  t_vax,
  vt1=0,
  vt2=0,
  vt3=0,
  gt1=1,
  gt2=1,
  gt3=1
){

  if (comp_type == "Metabolic"){
    mt <- c(1,1,2)
  } else {
    mt <- c(1,2,3)
  }

  setup_sim_dir(
      nstrain=3, 
      nhost=NHOST, 
      strain_id=seq(1,2,3), 
      st=c(1,2,2), 
      rt=c(rt1,rt2,rt3), 
      mt=mt,
      vt=c(vt1,vt2,vt3), 
      gt=c(gt1,gt2,gt3),
      t_max=t_max,
      t_vax=t_vax,
      beta=BETA,
      mu=MU,
      k_g=CROSS)
  
  write_res_pars(p_tau=p_tau,cost_res=cost,tau=tau,fname_res=res_file)
  write_sero_pars(alpha=c(alpha1,alpha2),fname_sero=sero_file)
  ## Run
  v0_res <- run_simulations(res_file, sero_file, niter=NITER, label=comp_type)
  return(v0_res)
}


#=====================================#
# PLOTTING FUNCTIONS                  #
#=====================================#
## Plot output of run_simulations
plot_trajectories <- function(dat, plot_by="strain", with_total=FALSE, plot_rel_freq=TRUE){

  dat <- dat %>%
      dplyr::select(time, iter, starts_with("strain"), label) %>%
      reshape2::melt(id=c("time", "iter","label"), 
      variable.name="id", 
      value.name="freq") %>%
      group_by(time, iter,label) %>%
      mutate(tot_freq=sum(freq)) %>%
      mutate(rel_freq=freq/tot_freq) %>%
      ungroup() 

  if (plot_rel_freq) { 
    dat$y <- dat$rel_freq
  } else {
    dat$y <- dat$freq
  }  
  strain_pars <- read.csv("strain_pars.csv")
  strain_pars$id <- paste0("strain ", seq(1,nrow(strain_pars)))
  if (plot_by == "serotype") {
    # Group strains by ST
  dat <- dat %>%
        inner_join(strain_pars,by="id") %>%
        group_by(id_s, time, label, iter) %>%
        summarize(y=sum(y)) %>%
        group_by(id_s,time,label) %>%
        summarize(p50 = quantile(y, prob=0.5),
            p05 = quantile(y, prob=0.05),
            p95 = quantile(y, prob=0.95))
   dat$id_s <- factor(dat$id_s)
    gg <- ggplot(dat, aes(x=time/365))+
          geom_line(aes(y=p50, col=id_s, group=id_s), linewidth=1.5)+
          geom_ribbon(aes(ymin=p05, ymax=p95, fill=id_s, group=id_s), alpha=0.3)+
          facet_wrap(~label)+
        guides(fill="none")
  }  else if (plot_by == "resistance"){
     dat <- dat %>%
        inner_join(strain_pars,by="id") %>%
        group_by(id_r, time, label, iter) %>%
        summarize(y=sum(y)) %>%
        group_by(id_r,time,label) %>%
        summarize(p50 = quantile(y, prob=0.5),
            p05 = quantile(y, prob=0.05),
            p95 = quantile(y, prob=0.95))
   dat$id_r <- factor(dat$id_r)
    gg <- ggplot(dat, aes(x=time/365))+
          geom_line(aes(y=p50, col=id_r, group=id_r), linewidth=1.5)+
          geom_ribbon(aes(ymin=p05, ymax=p95, fill=id_r, group=id_r), alpha=0.3)+
          facet_wrap(~label)+
        guides(fill="none")
  } else {
    gg <- dat %>%
      group_by(id, time,label) %>%
      summarize(p50 = quantile(y, prob=0.5),
          p05 = quantile(y, prob=0.05),
          p95 = quantile(y, prob=0.95)) %>%
      ggplot(aes(x=time/365))+
        geom_line(aes(y=p50, col=id), linewidth=1.5)+
        geom_ribbon(aes(ymin=p05, ymax=p95, fill=id), alpha=0.3)+
        facet_wrap(~label)+
      guides(fill="none")

    if (with_total) {
      dat <- dat %>%
        group_by(time,iter,label) %>%
        summarize(tot_freq=sum(freq)) %>%
        group_by(time,label) %>%
        summarize(p50 = quantile(tot_freq, prob=0.5),
            p05 = quantile(tot_freq, prob=0.05),
            p95 = quantile(tot_freq, prob=0.95))
      gg <- gg + 
        geom_line(data=dat, aes(x=time/365, y=p50), linewidth=1.5, linetype="dashed")
    }
  }
  
  return(gg)
  
}

## Plot strain competition network
plot_network <- function(strain_pars,fname){
  pal <- categorical_pal(2)
  sero_pal <- diverging_pal(8)
  nodes <- strain_pars %>%
    mutate(id=id, label=parse_number(id), 
      color = case_when(
        id_r==1~pal[1],
        id_r==0~pal[2],
        TRUE ~ pal[3]),
      label.font=2) %>%
    dplyr::select(id, id_s, id_v, label, color, label.font) %>%
    mutate(r=row_number())


  st_list <- lapply(unique(nodes$id_s), function(s) 
    {filter(nodes, id_s == s)$r} )

  st_col_list <- sapply(unique(nodes$id_s), 
    function(s){
      if (all(filter(nodes, id_s == s)$id_v == 0)){
        return(sero_pal[4])}
      else { return(sero_pal[7])}
    })


  edges <- make_edges_by_strain(strain_pars) %>%
    mutate(color=case_when(type=="meta" ~ "black"))
  graph <- graph_from_data_frame(edges, directed=FALSE, vertices=nodes)
  pdf(fname)
  plot(graph,layout=layout_with_fr(graph), 
    mark.groups=st_list,
    edge.width=2,
    mark.border=NA,mark.col=st_col_list, mark.shape=1/2, mark.expand=50,
    margin=1)

  legend("topleft", legend= c("Resistant", "Sensitive"), 
    col=pal, bty="n", pch=c(19,19), horiz=F)
  legend("topright", legend= c("VT", "NVT"), 
   col=c(sero_pal[7], sero_pal[4]), bty="n", pch=c(15,15), horiz=F)
  dev.off()
  return()
}

make_edges_by_strain <- function(df){
  edges <- data.frame(from=c(), to=c(), type=c())  
  for (s in df$id) {
    # get all edges for that sero based on meta type
    mts <- unique(filter(df, id==s)$id_m)
    df2 <- filter(df, id_m %in% mts)
    from_nodes <- unique(filter(df2, id != s)$id)
    # remove nodes that have already been checked
    from_nodes <- from_nodes[!(from_nodes %in% edges$to)]
    if (length(from_nodes) > 0){
      edges <- rbind(edges, data.frame(to=s, from=from_nodes, type="meta"))
    }

   # also ST
    st <- unique(filter(df, id==s)$id_s)
    df2 <- filter(df, id_s == st)
    from_nodes <- unique(filter(df2, id_s != s)$id)
    # remove nodes that have already been checked
    from_nodes <- from_nodes[!(from_nodes %in% edges$to)]
    if (length(from_nodes) > 0){
      edges <- rbind(edges, data.frame(to=s, from=from_nodes, type="sero"))
    }
  }
  return(edges)
}
#=====================================#
# Helpers for writing parameter files #
#=====================================#
write_constant_pars <- function(nstrain, # number of strains
                                nhost, # number of hosts
                                beta, # bwh transmission rate
                                t_max, # max simulation time
                                t_vax, # vaccination time
                                kappa, # growth rate of sensitive strain
                                eps_x, # wh death rate
                                eps_a, # immune decay rate
                                rho, # minimum wh pathogen density
                                theta,
                                mu,#migration rate
                                k_g){ #cross-immunity
    const_pars <- data.frame(t_max=t_max, t_vax=t_vax,
      nhost=nhost, nstrain=nstrain, kappa=kappa, beta=beta, theta=theta, eps_x=eps_x,
      eps_a=eps_a, rho=rho, mu=mu, k_g=k_g)
    # Note: multiabm_samplestrains will look for a file called "const_pars.csv" in the working directory
    write.table(const_pars, "const_pars.csv", row.names=FALSE, sep=",")
}

write_res_pars <- function(p_tau, # treatment coverage (pop level)
                           cost_res, # wh growth cost of resistance
                           tau, # wh treatment rate
                           fname_res="res_pars.csv"){#default name for res_pars
  # Note: multiabm_samplestrains takes fname_res as its 2nd argument
  write.table(t(c(p_tau, cost_res,tau)), fname_res, row.names=FALSE, col.names=FALSE, sep=",") 
}

write_sero_pars <- function(alpha, #nsero x 1 array 
                            fname_sero="sero_pars.csv"){
  # Note: multiabm_samplestrains takes fname_sero as its 3rd argument
  write.table(alpha, fname_sero, row.names=FALSE, col.names=FALSE, sep=",") 
}

write_strain_pars <- function(nstrain,
                              strain_id,
                              st, #serotype
                              rt, #resistance type
                              mt, #metabolic type/sequence type
                              vt, #vaccine type
                              gt){ #sero group
  # Check that 
  strain_pars <- data.frame(id_m=mt, id_s=st, id_r=rt, id_v=vt, id_g=gt)
  # Note: multiabm_samplestrains will look for a file called "strain_pars.csv" in working directory
  write.table(strain_pars, "strain_pars.csv", row.names=FALSE, sep=",")
}


