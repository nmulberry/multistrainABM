setwd("../")
source("setup.R")

make_edges_by_strain <- function(df){
  # SAME ERA EDGES
  edges <- data.frame(from=c(), to=c(), type=c())  
  for (s in df$profile) {
    # get all edges for that strain based on meta type
    mts <- unique(filter(df, profile==s)$MT)
    df2 <- filter(df, MT %in% mts)
    from_nodes <- unique(filter(df2, profile != s)$profile)
    # remove nodes that have already been checked
    from_nodes <- from_nodes[!(from_nodes %in% edges$to)]
    if (length(from_nodes) > 0){
      edges <- rbind(edges, data.frame(to=s, from=from_nodes, type="meta"))
    }

   # by serotype
    at <- unique(filter(df, profile==s)$AT)
    df2 <- filter(df, AT == at)
    from_nodes <- unique(filter(df2, profile != s)$profile)
    # remove nodes that have already been checked
    from_nodes <- from_nodes[!(from_nodes %in% edges$to)]
    if (length(from_nodes) > 0){
      edges <- rbind(edges, data.frame(to=s, from=from_nodes, type="sero"))
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
regions <- c("CONNECTICUT", "MARYLAND", "GEORGIA", "TENNESSEE") #at least 200 samples

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

st_gpsc_comp <- data %>%
  group_by(GPSC, Vaccine_period, ST=In_silico_ST, Region) %>%
  summarize(count=n())  %>%
  group_by(Vaccine_period) %>%
  mutate(freq=count/sum(count)) %>%
  ungroup() %>%
  group_by(Region)# %>%
 # tidyr::complete(Vaccine_period, GPSC, ST, fill=list(freq=0, count=0))

st_gpsc_comp$GPSC <- as.factor(st_gpsc_comp$GPSC)


test <- st_gpsc_comp %>%
  group_by(ST, Vaccine_period) %>%
  summarize(num_gpscs=n())

ggplot(st_gpsc_comp, aes(x=ST, y=freq, fill=GPSC))+
  geom_col()+
  facet_grid(Vaccine_period ~ Region)+
  theme(legend.position="none")


## and seq types freqs too for comparison
scdat <- data %>%
  group_by(Region, Vaccine_period, ST=In_silico_ST) %>%
  summarize(count=n()) %>%
  group_by(Region, Vaccine_period) %>%
  mutate(freq=count/sum(count)) %>%
  ungroup() %>%
  group_by(Region) %>%
  tidyr::complete(Vaccine_period, ST, fill=list(freq=0, count=0))

scdat$GPSC <- as.factor(scdat$ST)
ggplot(scdat, aes(x=ST, y=freq, fill=Vaccine_period))+
  geom_col(position="dodge")+
  facet_wrap(~ Region, nrow=length(regions))+
  scale_fill_manual(values=freq_pal)+
  labs(x="Sequence Type", y="Relative Frequency", fill="")+
  theme(legend.position="bottom")+
  theme(axis.text.x = element_text(angle = 66, vjust=0.5))

ggsave(paste0(out_dir, "/ST_FREQS.pdf"), height=12, width=14)

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



### pre-post vax GPSCs across entire US
data2 <- data %>%
  group_by(Vaccine_period, ST=In_silico_ST, Region) %>%
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
data2$ST <- as.factor(data2$ST)
ggplot(data2, aes(x=PREPCV, y=POSTPCV))+
  geom_point(col="black", size=2)+
  geom_point(aes(col=ST), size=1.5)+
  theme(legend.position="none")+
  facet_wrap(~Region)+
  theme(legend.position="none")+
  geom_abline(slope=1, intercept=0)+
  xlim(c(0,0.25))+ylim(c(0,0.25))


ggsave(paste0(out_dir, "/ST_SCATTER_FREQS.pdf"), height=12, width=14)

## all
data_US <- filter(raw_data, Country=="UNITED STATES") %>%
 filter(Vaccine_period == "PREPCV" | Vaccine_period == "POSTPCV7-9YR")

data_US <- data_US %>%
  mutate(Vaccinetype = case_when(
    In_silico_serotype %in% c("4","6B","9V", "14","18C","19F","23F", "6A") ~ "VT", # pcv7 and 13
    TRUE ~ "NVT")) %>%
  mutate(Vaccine_period = case_when(
    Vaccine_period == "POSTPCV7-9YR" ~ "POSTPCV",
    TRUE ~ "PREPCV"))

data$Vaccine_period <- factor(data$Vaccine_period,
  levels=c("PREPCV", "POSTPCV"))

### pre-post vax GPSCs across entire US
data2 <- data_US %>%
  group_by(Vaccine_period, ST=In_silico_ST) %>%
  summarize(freq=n()) %>%
  ungroup() %>%
  group_by(Vaccine_period) %>%
  mutate(rel_freq=freq/sum(freq)) %>%
  ungroup() %>%
  dplyr::select(-c(freq)) %>%
  pivot_wider(names_from=Vaccine_period, values_from=rel_freq) %>%
  mutate(PREPCV=replace_na(PREPCV, 0)) %>%
  mutate(POSTPCV=replace_na(POSTPCV, 0))
#
data2$ST <- as.factor(data2$ST)
ggplot(data2, aes(x=PREPCV, y=POSTPCV))+
  geom_point(col="black", size=2)+
  geom_point(aes(col=ST), size=1.5)+
  theme(legend.position="none")+
  theme(legend.position="none")+
  geom_abline(slope=1, intercept=0)+
  xlim(c(0,0.18))+ylim(c(0,0.18))

#########################
# PREVAX NETWORKS

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

  data <- filter(data, Vaccine_period=="PREPCV")

  compdata <- data %>%
   dplyr::select(Era=Vaccine_period,
  #  MT=In_silico_ST,
    MT=GPSC,
    AT=In_silico_serotype,
    VT=Vaccinetype) %>%
    group_by(AT,VT,MT) %>%
    summarize(samples=n()) %>%
    ungroup() %>%
    mutate(profile=row_number())

  nodes <- compdata %>%
    mutate(rel_freq = samples/sum(samples)) %>%
    ungroup() %>%
    dplyr::select(id=profile, AT, VT, rel_freq) 

  edges <- make_edges_by_strain(compdata) %>%
      mutate(color=case_when(type=="meta" ~ "black"))

  st_list <- lapply(unique(nodes$AT), function(s) 
    {filter(nodes, AT == s)$id} )

  st_col_list <- sapply(unique(nodes$AT), 
    function(s){
      if (all(filter(nodes, AT == s)$VT == "NVT")){
        return(sero_pal[4])}
      else { return(sero_pal[7])}
    })

  nodes$AT_label <- nodes$AT
  nodes$AT_label[duplicated(nodes$AT_label)] <- ""
  graph <- graph_from_data_frame(edges, directed=FALSE, vertices=nodes)
  #v0 = which(degree(graph)==0)
  #graph = delete.vertices(graph, v0)
  set.seed(123)
  pdf(paste0(out_dir, "/prevax_network_GPSC_", region, ".pdf"))
  op <- par(family = "sans")
  plot(graph, 
     edge.width=2, 
      vertex.label=nodes$AT_label,  
      vertex.label.cex = 1.5,
      vertex.label.color="black",
      vertex.label.dist=rep(2, length(nodes$AT)),
      vertex.size=4, 
      mark.groups=st_list,
      mark.border=NA, 
      mark.col=st_col_list,
      mark.shape=1/2, 
      mark.expand=10)
  title(region)
  dev.off()
}



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
    #MT=In_silico_ST,
   MT=GPSC,
    AT=In_silico_serotype,
    VT=Vaccinetype) %>%
    group_by(AT,VT,MT,Era) %>%
    summarize(samples=n()) %>%
    ungroup() %>%
    group_by(Era) %>%
    mutate(rel_freq = samples/sum(samples)) %>%
    ungroup() %>%
    group_by(AT, VT, MT) %>%
    top_n(1, rel_freq) %>%
    ungroup() %>%
    mutate(profile=row_number())

  nodes <- compdata %>%
    dplyr::select(id=profile, AT, VT, Era) %>%
    mutate(col=case_when(Era=="PREPCV" ~ "#e5b8ab", TRUE ~ "#972D15"))

  edges <- make_edges_by_strain(compdata) %>%
      mutate(color=case_when(type=="meta" ~ "black"))

  st_list <- lapply(unique(nodes$AT), function(s) 
    {filter(nodes, AT == s)$id} )

  st_col_list <- sapply(unique(nodes$AT), 
    function(s){
      if (all(filter(nodes, AT == s)$VT == "NVT")){
        return(sero_pal[4])}
      else { return(sero_pal[7])}
    })

  nodes$AT_label <- nodes$AT
  nodes$AT_label[duplicated(nodes$AT_label)] <- ""
  graph <- graph_from_data_frame(edges, directed=FALSE, vertices=nodes)
  #v0 = which(degree(graph)==0)
  #graph = delete.vertices(graph, v0)
  set.seed(123)
  pdf(paste0(out_dir, "/comp_network_GPSC_", region, ".pdf"))
  op <- par(family = "sans")
  plot(graph, 
     edge.width=2, 
      vertex.color=nodes$col,
      vertex.label.cex = 1.5,
      vertex.label.color="black",
      vertex.label=nodes$AT_label,  
      vertex.label.dist=rep(2, length(nodes$AT)),
      vertex.size=4, 
      mark.groups=st_list,
      mark.border=NA, 
      mark.col=st_col_list,
      mark.shape=1/2, 
      mark.expand=10)
  title(region)
  dev.off()
}




