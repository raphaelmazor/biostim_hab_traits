library(tidyverse)

# assemblage_df<-read_csv("Data/TITAN/assemblage.responses.csv")
titan_df<-read_csv("Data/TITAN/titan_derived_traits_biostim.csv")
taxa_df<-read_csv("Data/Taxa/toxin producing taxa.csv")


taxa_df$InTitan<-
  sapply(taxa_df$Genus_orig, function(x){x %in% titan_df$Taxon})


taxa_df %>% clipr::write_clip()

titan_df %>% filter(Assemblage=="SBA") %>%
  select(Taxon) %>% distinct() %>% clipr::write_clip()

titan_df %>%
  filter(Taxon %in% taxa_df$Genus_orig) %>%
  write_csv("Output/titan_hab.csv")


plot_dat<- titan_df %>%
  filter(Taxon %in% taxa_df$Genus_orig)

ggplot(data=plot_dat,
       aes(x=Taxon, y=zenv.cp))+
  geom_point()+
  facet_wrap(~Stressor, scales="free")+
  coord_flip()


plot_dat_afdm<-plot_dat %>% 
  filter(Stressor=="AFDM") %>%
  arrange(zenv.cp) %>%
  mutate(Taxon=factor(Taxon, levels=Taxon))
afdm_hab_plot<-
ggplot(data=plot_dat_afdm,
       aes(x=Taxon, y=zenv.cp))+
  # geom_hline(data=plot.tn.sumdat, aes(yintercept=cp), linetype="dashed", color="white")+
  geom_linerange(aes(ymin=`5%`, ymax=`95%`), color="gray")+
  geom_point(aes( fill=Group, shape=Group), size=1.5)+
  # geom_hline(data=tn_cp_df, aes(yintercept=zenv.cp, color=Assemblage, linetype=Response), size=0.5)+scale_linetype_manual(values=c("dotted","dashed"), name="Mean response")+
  scale_shape_manual(values=c(24,25,21), name="Response")+
  scale_fill_manual(values=c("skyblue","tomato1","darkgreen"))+
  coord_flip()+
  theme_bw(base_size=10)+
  xlab("")+ylab("Changepoint")+
  ggtitle("Ash-Free Dry Mass (mg/cm2)")+
  scale_size_continuous("Indicator value\n(z-score)")+
  theme(axis.text.y = element_text(face = "italic", size=6), panel.grid = element_blank())+
  guides(fill = guide_legend(override.aes=list(shape=22, color=NA, size=5)),
         shape = guide_legend(override.aes=list(size=2)))
ggsave(afdm_hab_plot, filename="Output/afdm_hab_plot.png", height=6, width=6)


plot_dat_tn<-plot_dat %>% 
  filter(Stressor=="TN") %>%
  arrange(zenv.cp) %>%
  mutate(Taxon=factor(Taxon, levels=Taxon))
tn_hab_plot<-
  ggplot(data=plot_dat_tn,
         aes(x=Taxon, y=zenv.cp))+
  # geom_hline(data=plot.tn.sumdat, aes(yintercept=cp), linetype="dashed", color="white")+
  geom_linerange(aes(ymin=`5%`, ymax=`95%`), color="gray")+
  geom_point(aes( fill=Group, shape=Group), size=1.5)+
  # geom_hline(data=tn_cp_df, aes(yintercept=zenv.cp, color=Assemblage, linetype=Response), size=0.5)+scale_linetype_manual(values=c("dotted","dashed"), name="Mean response")+
  scale_shape_manual(values=c(24,25,21), name="Response")+
  scale_fill_manual(values=c("skyblue","tomato1","darkgreen"))+
  coord_flip()+
  theme_bw(base_size=10)+
  xlab("")+ylab("Changepoint")+
  ggtitle("Total N (mg/L)")+
  scale_size_continuous("Indicator value\n(z-score)")+
  theme(axis.text.y = element_text(face = "italic", size=6), panel.grid = element_blank())+
  guides(fill = guide_legend(override.aes=list(shape=22, color=NA, size=5)),
         shape = guide_legend(override.aes=list(size=2)))
ggsave(tn_hab_plot, filename="Output/tn_hab_plot.png", height=6, width=6)


plot_dat_tp<-plot_dat %>% 
  filter(Stressor=="TP") %>%
  arrange(zenv.cp) %>%
  mutate(Taxon=factor(Taxon, levels=Taxon))
tp_hab_plot<-
  ggplot(data=plot_dat_tp,
         aes(x=Taxon, y=zenv.cp))+
  # geom_hline(data=plot.tn.sumdat, aes(yintercept=cp), linetype="dashed", color="white")+
  geom_linerange(aes(ymin=`5%`, ymax=`95%`), color="gray")+
  geom_point(aes( fill=Group, shape=Group), size=1.5)+
  # geom_hline(data=tn_cp_df, aes(yintercept=zenv.cp, color=Assemblage, linetype=Response), size=0.5)+scale_linetype_manual(values=c("dotted","dashed"), name="Mean response")+
  scale_shape_manual(values=c(24,25,21), name="Response")+
  scale_fill_manual(values=c("skyblue","tomato1","darkgreen"))+
  coord_flip()+
  theme_bw(base_size=10)+
  xlab("")+ylab("Changepoint")+
  ggtitle("Total P (mg/L)")+
  scale_size_continuous("Indicator value\n(z-score)")+
  theme(axis.text.y = element_text(face = "italic", size=6), panel.grid = element_blank())+
  guides(fill = guide_legend(override.aes=list(shape=22, color=NA, size=5)),
         shape = guide_legend(override.aes=list(size=2)))
ggsave(tp_hab_plot, filename="Output/tp_hab_plot.png", height=6, width=6)

plot_dat_chla<-plot_dat %>% 
  filter(Stressor=="Chl-a") %>%
  arrange(zenv.cp) %>%
  mutate(Taxon=factor(Taxon, levels=Taxon))
chla_hab_plot<-
  ggplot(data=plot_dat_chla,
         aes(x=Taxon, y=zenv.cp))+
  # geom_hline(data=plot.tn.sumdat, aes(yintercept=cp), linetype="dashed", color="white")+
  geom_linerange(aes(ymin=`5%`, ymax=`95%`), color="gray")+
  geom_point(aes( fill=Group, shape=Group), size=1.5)+
  # geom_hline(data=tn_cp_df, aes(yintercept=zenv.cp, color=Assemblage, linetype=Response), size=0.5)+scale_linetype_manual(values=c("dotted","dashed"), name="Mean response")+
  scale_shape_manual(values=c(24,25,21), name="Response")+
  scale_fill_manual(values=c("skyblue","tomato1","darkgreen"))+
  coord_flip()+
  theme_bw(base_size=10)+
  xlab("")+ylab("Changepoint")+
  ggtitle("Chlorophyll-a (mg/m2)")+
  scale_size_continuous("Indicator value\n(z-score)")+
  theme(axis.text.y = element_text(face = "italic", size=6), panel.grid = element_blank())+
  guides(fill = guide_legend(override.aes=list(shape=22, color=NA, size=5)),
         shape = guide_legend(override.aes=list(size=2)))
ggsave(chla_hab_plot, filename="Output/chla_hab_plot.png", height=6, width=6)



plot_dat_map<-plot_dat %>% 
  filter(Stressor=="PCT_MAP") %>%
  arrange(zenv.cp) %>%
  mutate(Taxon=factor(Taxon, levels=Taxon))
map_hab_plot<-
  ggplot(data=plot_dat_map,
         aes(x=Taxon, y=zenv.cp))+
  # geom_hline(data=plot.tn.sumdat, aes(yintercept=cp), linetype="dashed", color="white")+
  geom_linerange(aes(ymin=`5%`, ymax=`95%`), color="gray")+
  geom_point(aes( fill=Group, shape=Group), size=1.5)+
  # geom_hline(data=tn_cp_df, aes(yintercept=zenv.cp, color=Assemblage, linetype=Response), size=0.5)+scale_linetype_manual(values=c("dotted","dashed"), name="Mean response")+
  scale_shape_manual(values=c(24,25,21), name="Response")+
  scale_fill_manual(values=c("skyblue","tomato1","darkgreen"))+
  coord_flip()+
  theme_bw(base_size=10)+
  xlab("")+ylab("Changepoint")+
  ggtitle("% macroalgal cover")+
  scale_size_continuous("Indicator value\n(z-score)")+
  theme(axis.text.y = element_text(face = "italic", size=6), panel.grid = element_blank())+
  guides(fill = guide_legend(override.aes=list(shape=22, color=NA, size=5)),
         shape = guide_legend(override.aes=list(size=2)))
ggsave(map_hab_plot, filename="Output/map_hab_plot.png", height=6, width=6)
