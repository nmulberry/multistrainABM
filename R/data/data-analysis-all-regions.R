setwd("../")
source("setup.R")


`%notin%` <- Negate(`%in%`)


make_edges_by_strain <- function(df){
  # SAME ERA EDGES
  edges <- data.frame(from=c(), to=c(), type=c())  
  for (s in df$profile) {
    # get era
    era <-  unique(filter(df, profile==s)$Era)
    df1 <-  filter(df, Era==era)
    # all nodes in same era get a (low) weighted edge (for plotting purposes)
    from_nodes <- unique(filter(df1, profile != s)$profile)
    # remove nodes that have already been checked
    from_nodes <- from_nodes[!(from_nodes %in% edges$to)]
    if (length(from_nodes) > 0){
      edges <- rbind(edges, data.frame(to=s, from=from_nodes, type="era", weight=.3, time=era))
    }

    # get all edges for that strain based on meta type
    mts <- unique(filter(df1, profile==s)$MT)
    df2 <- filter(df1, MT %in% mts)
    from_nodes <- unique(filter(df2, profile != s)$profile)
    # remove nodes that have already been checked
    from_nodes <- from_nodes[!(from_nodes %in% edges$to)]
    if (length(from_nodes) > 0){
      edges <- rbind(edges, data.frame(to=s, from=from_nodes, type="meta", weight=5, time=era))
    }

   # by serotype
    st <- unique(filter(df1, profile==s)$ST)
    df2 <- filter(df1, ST == st)
    from_nodes <- unique(filter(df2, profile != s)$profile)
    # remove nodes that have already been checked
    from_nodes <- from_nodes[!(from_nodes %in% edges$to)]
    if (length(from_nodes) > 0){
      edges <- rbind(edges, data.frame(to=s, from=from_nodes, type="sero", weight=5.5, time=era))
    }


  # NOW LOOK AT PUTATIVE "SWITCHES"

    if (era == "POSTPCV") {
        df2 <- filter(df, GPSC %in% mts & Era=="PREPCV")
        # make sure original sero not one of these
        if (!(st %in% df2$ST)) {
          from_nodes <- unique(filter(df2, profile != s)$profile)
          if (length(from_nodes) > 0){
            edges <- rbind(edges, data.frame(to=s, from=from_nodes, type="switch", weight=5, time="both"))
          }
        } else { # add original sero as an edge
          from_nodes <- unique(filter(df2, ST==st, profile != s)$profile)
          if (length(from_nodes) > 0){
            edges <- rbind(edges, data.frame(to=s, from=from_nodes, type="expand", weight=5, time="both"))
          }
        }
      }
    }  
  
  return(edges)
}

pal <- categorical_pal(8)[7:8]
sero_pal <- diverging_pal(8)

freq_pal <- c("PREPCV"="#e5b8ab","POSTPCV"="#972D15")
#######################################################
####################################################### 
raw_data <- read.csv("data/monocle-metadata.csv", sep=",")
results <- data.frame()
####---------- UNITED STATES-------------------####
regions <- c("CONNECTICUT", "MARYLAND", "MINNESOTA", "GEORGIA", "TENNESSEE") #at least 200 samples

## Some plots
data <- filter(raw_data, Country=="UNITED STATES" & Region %in% regions) %>%
 filter(Vaccine_period == "PREPCV" | Vaccine_period == "POSTPCV7-9YR")

data <- data %>%
  mutate(Vaccinetype = case_when(
    In_silico_serotype %in% c("4","6B","9V", "14","18C","19F","23F", "6A") ~ "VT", # pcv7 and 13
    TRUE ~ "NVT")) %>%
  mutate(Vaccine_period = case_when(
    Vaccine_period == "POSTPCV7-9YR" ~ "POSTPCV",
    TRUE ~ "PREPCV"))

data$Vaccine_period <- factor(data$Vaccine_period,
  levels=c("PREPCV", "POSTPCV"))


## plot pre-post serotype freqs by region
stdat <- data %>%
  group_by(Region, Vaccine_period, Serotype=In_silico_serotype, Vaccinetype) %>%
  summarize(count=n()) %>%
  group_by(Region, Vaccine_period) %>%
  mutate(freq=count/sum(count)) %>%
  ungroup() %>%
  group_by(Region, Vaccinetype) %>%
  tidyr::complete(Vaccine_period, Serotype, fill=list(freq=0, count=0))

ggplot(stdat, aes(x=Serotype, y=freq, fill=Vaccine_period))+
  geom_col(position="dodge")+
  facet_grid(Region ~ Vaccinetype, scales = "free_x", space='free')+
  scale_fill_manual(values=freq_pal)+
  labs(x="Serotype", y="Relative Frequency", fill="")+
  theme(legend.position="bottom")+
  theme(axis.text.x = element_text(angle = 66, vjust=0.5))

ggsave(paste0(out_dir, "/SERO_FREQS.pdf"), height=12, width=14)

## Just sc freqs too for comparison
scdat <- data %>%
  group_by(Region, Vaccine_period, GPSC) %>%
  summarize(count=n()) %>%
  group_by(Region, Vaccine_period) %>%
  mutate(freq=count/sum(count)) %>%
  ungroup() %>%
  group_by(Region) %>%
  tidyr::complete(Vaccine_period, GPSC, fill=list(freq=0, count=0))

scdat$GPSC <- as.factor(scdat$GPSC)
ggplot(scdat, aes(x=GPSC, y=freq, fill=Vaccine_period))+
  geom_col(position="dodge")+
  facet_wrap(~ Region, nrow=length(regions))+
  scale_fill_manual(values=freq_pal)+
  labs(x="GPSC", y="Relative Frequency", fill="")+
  theme(legend.position="bottom")+
  theme(axis.text.x = element_text(angle = 66, vjust=0.5))

ggsave(paste0(out_dir, "/GPSC_FREQS.pdf"), height=12, width=14)


### pre-post vax GPSCs across entire US
data2 <- data %>%
  group_by(Vaccine_period, GPSC, Region) %>%
  summarize(freq=n()) %>%
  ungroup() %>%
  group_by(Vaccine_period,Region) %>%
  mutate(rel_freq=freq/sum(freq)) %>%
  ungroup() %>%
  dplyr::select(-c(freq)) %>%
  pivot_wider(names_from=Vaccine_period, values_from=rel_freq) %>%
  mutate(PREPCV=replace_na(PREPCV, 0)) %>%
  mutate(POSTPCV=replace_na(POSTPCV, 0))
#
data2$GPSC <- as.factor(data2$GPSC)
ggplot(data2, aes(x=PREPCV, y=POSTPCV))+
  geom_point(col="black", size=2)+
  geom_point(aes(col=GPSC), size=1.5)+
  theme(legend.position="none")+
  facet_wrap(~Region)+
  theme(legend.position="none")+
  geom_abline(slope=1, intercept=0)+
  xlim(c(0,0.25))+ylim(c(0,0.25))


ggsave(paste0(out_dir, "/GPSC_SCATTER_FREQS.pdf"), height=12, width=14)

#################################################################
#################################################################

for (region in regions) {
  data <- filter(raw_data, Country=="UNITED STATES" & Region==region) %>%
   filter(Vaccine_period == "PREPCV" | Vaccine_period == "POSTPCV7-9YR")

  # add Vaccinetypes
  data <- data %>%
    mutate(Vaccinetype = case_when(
      In_silico_serotype %in% c("4","6B","9V", "14","18C","19F","23F", "6A") ~ "VT", # pcv7 and 13
      TRUE ~ "NVT")) %>%
    mutate(Vaccine_period = case_when(
      Vaccine_period == "POSTPCV7-9YR" ~ "POSTPCV",
      TRUE ~ "PREPCV"))

  data$Vaccine_period <- factor(data$Vaccine_period,
    levels=c("PREPCV", "POSTPCV"))


  compdata <- data %>%
   dplyr::select(Era=Vaccine_period,
    MT=GPSC,
    GPSC=GPSC,
    ST=In_silico_serotype,
    VT=Vaccinetype) %>%
    group_by(ST,VT,MT,GPSC,Era) %>%
    summarize(samples=n()) %>%
    ungroup() %>%
    mutate(profile=row_number())

  nodes <- compdata %>%
    group_by(Era) %>%
    mutate(rel_freq = samples/sum(samples)) %>%
    ungroup() %>%
    dplyr::select(id=profile, ST, VT, GPSC, rel_freq, Era) %>%
    mutate(col = case_when(Era=="PREPCV" ~ freq_pal[["PREPCV"]], TRUE ~ freq_pal[["POSTPCV"]]))

  edges <- make_edges_by_strain(compdata) %>%
      mutate(color=case_when(type=="meta" ~ "black", type == "switch" ~ "black", type== "expand" ~ "black"))%>%
      mutate(lty=case_when(type=="meta" ~ 1, type == "switch" ~ 2, type== "expand" ~ 3, TRUE ~ 0))

  pre_nodes <- filter(nodes, Era=="PREPCV")
  post_nodes <- filter(nodes, Era=="POSTPCV")

  st_list <- c( lapply(unique(pre_nodes$ST), function(s) 
    {filter(pre_nodes, ST == s)$id} ), 
   lapply(unique(post_nodes$ST), function(s) 
    {filter(post_nodes, ST == s)$id} ))

  st_col_list <- c( sapply(unique(pre_nodes$ST), 
    function(s){
      if (all(filter(pre_nodes, ST == s)$VT == "VT")){
        return(sero_pal[7])}
      else { return(sero_pal[4])}
    }),
  sapply(unique(post_nodes$ST), 
    function(s){
      if (all(filter(post_nodes, ST == s)$VT == "VT")){
        return(sero_pal[7])}
      else { return(sero_pal[4])}
    }))
      
     
  graph <- graph_from_data_frame(edges, directed=FALSE, vertices=nodes)
  #v0 = which(degree(graph)==0)
  #graph = delete.vertices(graph, v0)
  set.seed(123)
  pdf(paste0(out_dir, "/network_", region, ".pdf"))
  #op <- par(family = "sans")
  plot(graph, layout = layout_with_fr,
     edge.width=2, 
      vertex.label=NA,  
      vertex.color=nodes$col, 
      vertex.size=log(nodes$rel_freq*1000), 
      mark.groups=st_list,
      mark.border=NA, 
      mark.col=st_col_list,
      mark.shape=1/2, 
      mark.expand=10, 
      edge.weight=edges$weight)
  dev.off()


  #######################
  # SOME NETWORK MEASURES
  pre_vax_samp <- sum(filter(compdata, Era=="PREPCV")$samples)
  post_vax_samp <-sum(filter(compdata, Era=="POSTPCV")$samples)
  pre_vax_meta <-  nrow(filter(edges, type=="meta" & time =="PREPCV"))
  post_vax_meta <-  nrow(filter(edges, type=="meta" & time == "POSTPCV"))
  nswitch <- nrow(filter(edges, type=="switch"))

  sero_freqs <- nodes %>% 
    group_by(ST, Era) %>% 
    summarize(rel_freq=sum(rel_freq)) 

  pre_vax_dom_sero <- filter(sero_freqs, Era=="PREPCV") %>% arrange(desc(rel_freq))
  pre_vax_dom_sero <- pre_vax_dom_sero$ST[1:5]
  post_vax_dom_sero <- filter(sero_freqs, Era=="POSTPCV") %>% arrange(desc(rel_freq))
  post_vax_dom_sero <- post_vax_dom_sero$ST[1:5]


  sc_freqs <- nodes %>% 
    group_by(GPSC, Era) %>% 
    summarize(rel_freq=sum(rel_freq)) 

  pre_vax_dom_sc <- filter(sc_freqs, Era=="PREPCV") %>% arrange(desc(rel_freq))
  pre_vax_dom_sc <- pre_vax_dom_sc$GPSC[1:5]
  post_vax_dom_sc <- filter(sc_freqs, Era=="POSTPCV") %>% arrange(desc(rel_freq))
  post_vax_dom_sc <- post_vax_dom_sc$GPSC[1:5]

  res1 <- data.table(location=region,
    `Pre-PCV Samples`=pre_vax_samp,
    `Post-PCV Samples`=post_vax_samp,
    `Pre-PCV GPSC-pairs`=pre_vax_meta/pre_vax_samp,
    `Post-PCV GPSC-pairs`=post_vax_meta/post_vax_samp,
    `Putative switches`=nswitch/(pre_vax_samp+post_vax_samp),
    `Pre-PCV Dominant Serotypes`= list(list(pre_vax_dom_sero)),
    `Post-PCV Dominant Serotypes`= list(list(post_vax_dom_sero)),
    `Pre-PCV Dominant GPSCs`= list(list(pre_vax_dom_sc)),
    `Post-PCV Dominant GPSCs`= list(list(post_vax_dom_sc)))

  results <- rbind(res1, results)

}


df <- knitr::kable(results, format = 'latex')
writeLines(df, paste0(out_dir, "/results_data_analysis.tex"))





