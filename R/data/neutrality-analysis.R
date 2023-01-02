
library("tidyverse")
library("reshape2")
library("igraph")
library("ggpubr")

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

strain_map <- readRDS("../data/strain_map_simplified.RDS")
freqs_raw <- readRDS("../data/strain_freqs_simplified.RDS")
freqs_raw$profile <- factor(freqs_raw$profile)
par_map <- strain_map %>%
  mutate(id_v = case_when(VT == "NVT" ~ 0, TRUE ~ 1))
par_map$id_m <- as.double(factor(par_map$MT), levels = unique(par_map$MT))
par_map$id_s <- as.double(factor(par_map$ST), levels = unique(par_map$ST))
par_map$id_r <- par_map$RT


dat <- freqs_raw %>%
  filter(Time != "Y3") %>%
  mutate(Time=case_when(Time=="Y0"~"pre-vax", TRUE~"post-vax")) %>%
  group_by(ST,VT, Time) %>%
  summarize(Freq=sum(Freq)) 

dat$Time <- factor(dat$Time, levels=c("pre-vax","post-vax"))

gg_nvt <- ggplot(filter(dat,VT=="NVT"), aes(x=ST, y=Freq, group=Time, fill=Time))+
    geom_col(position="dodge",col="black")+
    scale_fill_manual(values=c("#972D15", "#e5b8ab"))+
    labs(x="Serotype", y="Relative Frequency", fill="")+
    theme(legend.position="none")


gg_vt <- ggplot(filter(dat,VT=="VT"), aes(x=ST, y=Freq, group=Time, fill=Time))+
    geom_col(position="dodge",col="black")+
    scale_fill_manual(values=c("#972D15", "#e5b8ab"))+
    labs(x="Serotype", y="Relative Frequency", fill="")+
    theme(legend.position="right")

gg_all <- ggplot(dat, aes(x=ST, y=Freq, group=Time, fill=Time))+
    geom_col(position="dodge",col="black")+
    scale_fill_manual(values=c("#972D15", "#e5b8ab"))+
    labs(x="Serotype", y="Relative Frequency", fill="")+
    theme(legend.position="bottom")+
    facet_grid(~VT, scales="free", space="free")


dat_delta <- dat %>%
  filter(VT=="NVT") %>%
  pivot_wider(names_from="Time", values_from="Freq") %>%
  filter(`pre-vax` > 0) %>%
  mutate(diff=(`post-vax`-`pre-vax`)/`pre-vax`) %>%
  mutate(mark = case_when(
          ST %in% c("23B","23A","6C","34","7C","15B/C","15F","15A","19A","11A") ~ "Y"))

gg_delta <- ggplot(dat_delta, aes(y=ST, x=diff))+
  geom_col(col="black", fill="#972D15")+
  labs(y="Serotype", x="Relative change")+
 geom_text(aes(label = ifelse(mark=="Y", "*", "")), 
            position = position_dodge(width = .5), size=12) 



ggarrange(plotlist=list(
  ggarrange(plotlist=list(ggplot()+theme_minimal(), gg_delta),
   nrow=1, labels=c("A","B"), font.label=list(size=16),widths=c(3,2)),
gg_all),nrow=2,labels=c("","C"), font.label=list(size=16))


#ggsave("figures/deltafreq_nvt.pdf", height=10,width=14)
