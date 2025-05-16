library(tidyverse)
#Evaluation of sites needing channel engineering summary
library(DBI) # needed to connect to data.dfbase
library(dbplyr) # needed to connect to data.dfbase
library(RPostgreSQL) # needed to connect to our data.dfbase
library(rstudioapi) # just so we can type the password as we run the script, so it is not written in the clear
library(tidyverse)
library(lubridate)

#########Import data#########
# taxa_df<-read_csv("Data/Taxa/toxin producing taxa.csv")

con <- dbConnect(
  PostgreSQL(),
  host = "geobiology.cottkh4djef2.us-west-2.rds.amazonaws.com",
  dbname = 'smc',
  user='smcread',
  password='1969$Harbor'
)

# #Admin
# con <- dbConnect(
#   PostgreSQL(),
#   host = "geobiology.cottkh4djef2.us-west-2.rds.amazonaws.com",
#   dbname = 'smc',
#   user='sde',
#   password='Q6DC2pgGWUxtKrmz'
# )

######import######
lustations_df<-dbGetQuery(con, sql(
  "SELECT *
FROM
sde.lu_stations"
)) %>% as_tibble()

algae_df<-dbGetQuery(con, sql(
  "SELECT *
FROM
sde.unified_algae"
)) %>% as_tibble()

nutrients_df<-dbGetQuery(con, sql(
  "SELECT *
FROM
sde.analysis_chem_nutrients_0"
)) %>% as_tibble()

chem_df<-dbGetQuery(con, sql(
"SELECT *
FROM
sde.unified_chemistry
WHERE
sde.unified_chemistry.analytename LIKE ('Chlorophyll a%') OR sde.unified_chemistry.analytename LIKE ('Ash%')"
)) %>% as_tibble()

phab_df<-dbGetQuery(con, sql(
  "SELECT *
FROM
sde.analysis_phabmetrics"
)) %>% as_tibble()

write_csv(algae_df,"NotForGit/Step0/algae_df.csv")
write_csv(chem_df,"NotForGit/Step0/chem_df.csv")
write_csv(lustations_df,"NotForGit/Step0/lustations_df.csv")
write_csv(nutrients_df,"NotForGit/Step0/nutrients_df.csv")
write_csv(phab_df,"NotForGit/Step0/phab_df.csv")
###########################################


#####Clean and manipulate data####
algae_df %>%
  print(width=Inf)
algae_df2<-algae_df %>%
  transmute(
    stationcode,
    sampledate=as_date(sampledate),
    collectiondevicename,
    grabsize,
    sampletypecode,
    finalid, 
    baresult, 
    result = as.numeric(result), 
    unitname
    ) %>%
  filter(result>0) %>%
  inner_join(lustations_df %>%
               select(stationcode=stationid, masterid))

algae_df2$unitname %>% unique()

tox_alg<-algae_df2 %>% 
  filter(
    finalid %in% taxa_df$Genus_orig,
    # sampletypecode!="Integrated",
  unitname!="count"
)

nutrients_df %>%
  filter(str_detect(masterid, "\t"))
nutrients_df2<-nutrients_df %>%
  mutate(
    masterid=case_when(masterid=="\t\r\n204COY425"~"204COY425", T~masterid),
    sampledate=as_date(sampledate)) %>%
  # filter(sampletypecode %in% c("grab")) %>%
  filter(total_n_mgl_method %in% c("reported","calculated")) %>%
  group_by(masterid, sampledate, analyte="TN_mgL") %>%
  summarize(result = mean(total_n_mgl, na.rm=T)) %>%
  bind_rows(
    nutrients_df %>%
      mutate(
        masterid=case_when(masterid=="\t\r\n204COY425"~"204COY425", T~masterid),
        sampledate=as_date(sampledate)) %>%
      # filter(sampletypecode %in% c("grab")) %>%
      group_by(masterid, sampledate, analyte="TP_mgL") %>%
      summarize(result = mean(total_p_mgl, na.rm=T))
  )
chem_df %>% group_by(analytename, unit) %>% tally()

# chem_df2<-
  chem_df %>%
  filter(stationcode=="SMC00345") %>%
  filter(stationcode!="000NONPJ") %>%
    filter(unit %in%c("mg/m2","ug/cm2","mg/cm2","mg/cm²","ug/cm²","g/m2"),
           sampletypecode %in% c("Integrated","Grab","Composite","integrated","IntegratedSplit")) %>%
  mutate( sampledate=as_date(sampledate),
          result=as.numeric(result),
          result=case_when(
            resqualcode %in% c("ND")~0,
            analytename=="Ash Free Dry Mass" & unit %in% c("g/m2")~result * 1, 
            analytename=="Ash Free Dry Mass" & unit %in% c("mg/cm2","mg/cm²")~result*10,
            str_detect(analytename, "Chlorophyll a") & unit %in% c("mg/m2")~result*1,
            str_detect(analytename, "Chlorophyll a") & unit %in% c("ug/cm2","ug/cm²")~result*10,
            T~NA_real_)
          ) %>%
    
  inner_join(lustations_df %>% select(stationcode=stationid, masterid)) %>%
  group_by(masterid, sampledate, analytename) %>%
  summarize(result=mean(result, na.rm=T)) %>%
  ungroup()
  # filter(sampletypecode %in% c("grab","integrated"))

chem_df2 %>% filter(is.na(result))
chem_df %>% filter(stationcode=="SMC00345") %>% as.data.frame()

plot_dat<-algae_df2 %>%
  rename(biovolume=result) %>%
  inner_join(nutrients_df2) %>%
  mutate(ToxGen = finalid %in% taxa_df$Genus_orig,
         ToxGenName = case_When(ToxGen~finalid,T~"Other"))

ggplot(plot_dat, aes(x=result, y=biovolume))+
  geom_point()+
  facet_grid(analyte~., scales="free")

tox_taxa_present = taxa_df$Genus_orig[taxa_df$Genus_orig %in% algae_df2$finalid] %>% unique()

alg_tax_summary<-algae_df2 %>% 
  group_by(masterid, sampledate) %>% 
  summarize(tot_biovol=sum(result),
            tot_biovol_tox=sum(result[finalid %in% tox_taxa_present]),
            tot_taxa=length(finalid %>% unique()),
            tot_taxa_tox=length(intersect(finalid, tox_taxa_present))) %>%
  ungroup() %>%
  crossing(finalid=tox_taxa_present) %>%
  left_join(algae_df2 %>%
              select(masterid, sampledate, finalid, biovol=result)) %>%
  mutate(biovol=case_when(is.na(biovol)~0,T~biovol)) %>%
  inner_join(
    nutrients_df2
    )


ggplot(alg_tax_summary %>%
         select(-finalid, -biovol) %>%
         distinct(),
       aes(x=result, y=tot_biovol_tox)
         )+
  geom_point(aes(color=tot_taxa_tox))+
  facet_wrap(~analyte, scales="free")+
  scale_y_sqrt()+
  scale_x_sqrt()

total_biovol<-
ggplot(alg_tax_summary %>%
         select(-finalid, -biovol) %>%
         distinct() %>%
         filter(tot_biovol_tox>0),
       aes(x=result, y=tot_biovol_tox)
)+
  geom_point(aes(color=tot_taxa_tox %>% as.factor()))+
  scale_color_viridis_d( name="# toxic taxa")+
  facet_wrap(~analyte, scales="free")+
  scale_y_log10("Total biovolume of toxic genera")+
  scale_x_log10("Concentration")

ggsave(total_biovol, filename="Output/Scatterplots/total_biovol.png",
       height=4, width=6.5)

ggplot(alg_tax_summary %>%
         filter(biovol>0,
                analyte=="TN_mgL"),
       aes(x=result, y=biovol)
)+
  geom_point(aes(color=tot_biovol  ))+
  scale_color_viridis_c(name="Total biovolume")+
  # scale_color_viridis_d( name="# toxic taxa")+
  facet_wrap(~finalid, scales="free")+
  scale_y_log10("Biovolume")+
  scale_x_sqrt("Total N (mg/L)")
