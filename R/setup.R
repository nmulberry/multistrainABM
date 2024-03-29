packages <- c("tidyverse", "deSolve", "reshape2", "ggpubr", "RColorBrewer", "viridis",
  "ggridges","igraph", "furrr", "wesanderson") 

for (p in packages){
  if (!(p %in% installed.packages())){
    install.packages(p)
  }
}
library(furrr)
library(tidyverse)
library(deSolve)
library(reshape2)
library(ggpubr)
library(RColorBrewer)
library(viridis)
library(ggridges)
library(igraph)
library(data.table) 
library(wesanderson)
source('functions.R')
# furrr options
no_cores <- availableCores()-1
plan(multicore, workers = no_cores)

# ggplot style
tsize <- 16
theme_set(theme_light(base_size = tsize)+
  theme(strip.text.x = element_text(size=tsize, color="black"),
    strip.text.y = element_text(size=tsize, color="black"),
    legend.position="bottom",
  strip.background = element_rect(
     color="white", fill="white"
     ),
    ))
lsize <- tsize

# some colour palettes
sens_res_pal <- c("strain 2"="#E69F00","strain 1"="#56B4E9")
vax_pal <- c("#FEE0B6", "#8073AC")
freq_pal <- c("PREPCV"="#e5b8ab","POSTPCV"="#972D15")

## R model for illustrating within host dynamics only
withinhost_model <- function(time,state,parameters){
   with(as.list(c(state, parameters)), {      
    #-------DEs----------#
    x <- state[1:nstrain]
    a <- state[(nstrain+1):(2*nstrain)]
    st <-  parameters[grepl("^st", names(parameters))]
    mt <- parameters[grepl("^mt", names(parameters))]
    rt <-  parameters[grepl("^rt", names(parameters))]
    # Compute total carriage by sero and meta for each strain
    s <- rep(0, nstrain)
    m <- rep(0, nstrain)
    for (i in 1:nstrain) {
      for (j in 1:nstrain) {
        if (st[i] == st[j]) {
          # include i=j
          s[i] = s[i] + x[j]
        }
        if (mt[i] == mt[j]) {
          # include i=j
          m[i] = m[i] + x[j]
        }
      }
    }
    dx <- (kappa*(1-cost_res)*(rt)+kappa*(1-rt))*x*(1-m)-alpha*x*a-tau*x*(1-rt)*p_tau-eps_x*x
    da <- s^2/(rho^2 + s^2)-eps_a*a + a^2/(theta^2+a^2)
    # correct carrying capcity
    dx[x<rho] <- -x[x<rho]
return(list(c(dx, da)))
  })
}


