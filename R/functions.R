#=====================================#
# SIMULATION FUNCTIONS                #
#=====================================#
# set up const par files
setup_sim_dir <- function(nstrain, # number of strains
                          nhost, # number of hosts
                          strain_id, # nstrain x1 array of strain ids 
                          A, # nstrain x 1 array of antigenic types (integer-values starting at 1)
                          G, # nstrain x 1 array of gpscs (integer-values starting at 1)
                          O,  # nstrain x 1 array of accessory types (integer-values)
                          V, # nstrain x 1 array of vaccine types (*BY SERO*)
                          init, # initial strains (binary values init or not)
                          t_max, # max simulation time
                          t_vax, # vaccination time
                          t_g, # transformation gain rate
                          t_l, #transformation loss rate
                          beta=0.5, # bwh transmission rate
                          eps_x=0.03, # wh death rate
                          eps_a=0.02, # immune decay rate
                          rho=0.01, # minimum wh pathogen density
                          theta=70
  ){
  #----- check input-----#
  nsero <- length(unique(A))

  if (!(length(A)==nstrain)){
    stop("Check dimensions of A")
  }
  if (!(length(G)==nstrain)){
    stop("Check dimensions of G")
  }
  if (!(length(O)==nstrain)){
    stop("Check dimensions of O")
  } 
  if (!(length(V)==nstrain)){
    stop("Check dimensions of V")
  }
  # ------SETUP-------#
  write_constant_pars(nstrain=nstrain,nhost=nhost,t_max=t_max,t_vax=t_vax,
    eps_x=eps_x, eps_a=eps_a, rho=rho,theta=theta,beta=beta, t_g=t_g, t_l=t_l)
  write_strain_pars(nstrain=nstrain,strain_id=strain_id,A=A, G=G, O=O, V=V, init=init)
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
run_simulations <- function(sc_file, sero_file, ag_file, niter,label="None"){
  if (!(file.exists(sc_file))){
    "sc parameter file not found"
    return(NULL)
  } else if (!(file.exists(sero_file))){
    "sero parameter file not found"
    return(NULL)
  } else if (!(file.exists(ag_file))){
    "ag parameter file not found"
    return(NULL)
  } else {
    sim <- function(iter) {
         tmp <- paste0("tmp_out_", iter, ".csv")
         system(paste(multiabm, sc_file, sero_file, ag_file, ">", tmp), wait=TRUE)
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
            p05 = quantile(y, prob=0.25),
            p95 = quantile(y, prob=0.75))
   dat$id_s <- factor(dat$id_s)
    gg <- ggplot(dat, aes(x=time/365))+
          geom_line(aes(y=p50, col=id_s, group=id_s), size=1.5)+
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
            p05 = quantile(y, prob=0.25),
            p95 = quantile(y, prob=0.75))
   dat$id_r <- factor(dat$id_r)
    gg <- ggplot(dat, aes(x=time/365))+
          geom_line(aes(y=p50, col=id_r, group=id_r), size=1.5)+
          geom_ribbon(aes(ymin=p05, ymax=p95, fill=id_r, group=id_r), alpha=0.3)+
          facet_wrap(~label)+
        guides(fill="none")
  } else {
    gg <- dat %>%
      group_by(id, time,label) %>%
      summarize(p50 = quantile(y, prob=0.5),
          p05 = quantile(y, prob=0.25),
          p95 = quantile(y, prob=0.75)) %>%
      ggplot(aes(x=time/365))+
        geom_line(aes(y=p50, col=id), size=1.5)+
        geom_ribbon(aes(ymin=p05, ymax=p95, fill=id), alpha=0.3)+
        facet_wrap(~label)+
      guides(fill="none")

    if (with_total) {
      dat <- dat %>%
        group_by(time,iter,label) %>%
        summarize(tot_freq=sum(freq)) %>%
        group_by(time,label) %>%
        summarize(p50 = quantile(tot_freq, prob=0.5),
            p05 = quantile(tot_freq, prob=0.25),
            p95 = quantile(tot_freq, prob=0.75))
      gg <- gg + 
        geom_line(data=dat, aes(x=time/365, y=p50), size=1.5, linetype="dashed")
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
                                t_g,
                                t_l,
                                eps_x, # wh death rate
                                eps_a, # immune decay rate
                                rho, # minimum wh pathogen density
                                theta){
    const_pars <- data.frame(t_max=t_max, t_vax=t_vax,
      nhost=nhost, nstrain=nstrain, beta=beta, theta=theta, eps_x=eps_x,
      eps_a=eps_a, rho=rho, t_g=t_g, t_l=t_l)
    # Note: multiabm_samplestrains will look for a file called "const_pars.csv" in the working directory
    write.table(const_pars, "const_pars.csv", row.names=FALSE, sep=",")
}

write_ag_pars <- function(kappa, #nsero x 1 array 
                            fname){
  # Note: multiabm_samplestrains takes fname_sero as its 3rd argument
  write.table(kappa, fname, row.names=FALSE, col.names=FALSE, sep=",") 
}

write_sc_pars <- function(kappa, #nsero x 1 array 
                            fname){
  # Note: multiabm_samplestrains takes fname_sero as its 3rd argument
  write.table(kappa, fname, row.names=FALSE, col.names=FALSE, sep=",") 
}

write_sero_pars <- function(alpha, #nsero x 1 array 
                            fname){
  # Note: multiabm_samplestrains takes fname_sero as its 3rd argument
  write.table(alpha, fname, row.names=FALSE, col.names=FALSE, sep=",") 
}

write_strain_pars <- function(nstrain,
                              strain_id,
                              A,
                              G,
                              O,
                              V,
                              init){
  # Check that 
  strain_pars <- data.frame(id_m=G, id_s=A, id_o=O, id_v=V, id_init=init)
  # Note: multiabm_samplestrains will look for a file called "strain_pars.csv" in working directory
  write.table(strain_pars, "strain_pars.csv", row.names=FALSE, sep=",")
}


