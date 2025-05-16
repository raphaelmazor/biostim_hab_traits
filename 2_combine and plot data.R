library(tidyverse)

algae_df2<-read_csv("NotForGit/Step1/algae_df2.csv") %>%
  rename(biovolume=result)
agg_bs_data_df<-read_csv("NotForGit/Step1/agg_bs_data_df.csv")


agg_df<-algae_df2 %>%
  filter(Toxic) %>%
  inner_join(agg_bs_data_df)


taxa_vs_indicators_plot<-
ggplot(data=agg_df,
       aes(x=result, y=biovolume))+
  geom_point()+
  facet_grid(finalid~analytename, scales="free")+
  geom_smooth()+
  scale_y_log10()

ggsave(taxa_vs_indicators_plot, filename="Output/taxa_vs_indicators_plot.png",
       
       height=9, width=6.5)
