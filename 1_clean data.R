library(tidyverse)

#Read in files
lustations_df<-read_csv("NotForGit/Step0/lustations_df.csv")
algae_df<-read_csv("NotForGit/Step0/algae_df.csv")
chem_df<-read_csv("NotForGit/Step0/chem_df.csv")
nutrients_df<-read_csv("NotForGit/Step0/nutrients_df.csv")
phab_df<-read_csv("NotForGit/Step0/phab_df.csv")
taxa_df<-read_csv("Data/Taxa/toxin producing taxa.csv")


#add masterid
#change datetime to date
#go through samples

algae_df2<-algae_df %>%
  filter(!sampletypecode %in% c("Integrated", "Qualitative"),
         unitname=="um3/cm2",
         !is.na(result)
  ) %>%
  inner_join(
    lustations_df %>%
      select(stationcode=stationid, masterid)
  ) %>%
  transmute(
    masterid,
    alg_replicate=replicate,
    sampledate=as_date(sampledate),
    sampletypecode,
    finalid, 
    result, 
    # baresult, 
    unitname
  ) %>%
  mutate(Toxic = finalid %in% taxa_df$Genus_orig)

chem_df2<-chem_df %>%
  inner_join(
    lustations_df %>%
      select(stationcode=stationid, masterid)
  ) %>%
  filter(stationcode!="000NONPJ") %>%
  filter(unit %in%c("mg/m2","ug/cm2","mg/cm2","mg/cm²","ug/cm²","g/m2"),
         sampletypecode %in% c("Integrated","Grab","Composite","integrated","IntegratedSplit")) %>%
  transmute(
    masterid,
    chem_fieldreplicate=fieldreplicate, chem_labreplicate=labreplicate,
    sampledate=as_date(sampledate),
    analytename = case_when(analytename=="Ash Free Dry Mass"~"Ash Free Dry Mass",
                            str_detect(analytename,"Chlorophyll")~"Chlorophyll a",
                            T~"OTHER"),
    result=as.numeric(result),
    result=case_when(
      resqualcode %in% c("ND")~0,
      analytename=="Ash Free Dry Mass" & unit %in% c("g/m2")~result * 1, 
      analytename=="Ash Free Dry Mass" & unit %in% c("mg/cm2","mg/cm²")~result*10,
      str_detect(analytename, "Chlorophyll a") & unit %in% c("mg/m2")~result*1,
      str_detect(analytename, "Chlorophyll a") & unit %in% c("ug/cm2","ug/cm²")~result*10,
      T~NA_real_),
    unit=case_when(analytename=="Ash Free Dry Mass"~"g/m2",
                   analytename=="Chlorophyll a"~"mg/m2",
                   T~"OTHER")
  ) %>%
  group_by(masterid,
           sampledate,
           chem_fieldreplicate,
           analytename, unit) %>%
  summarize(result=mean(result, na.rm=T)) %>%
  group_by(masterid,
           sampledate,
           analytename, unit) %>%
  summarize(result=mean(result, na.rm=T)) %>%
  ungroup() 

nutrients_df2<-nutrients_df %>%
  mutate(masterid = case_when(masterid=="\n204COY425"~"204COY425",
                              T~masterid)) %>%
  #Nitrogen
  filter(total_n_mgl_method %in% c("reported","calculated"),
         !is.na(total_n_mgl)) %>%
  transmute(masterid, 
            sampledate= as_date(sampledate),
            fieldreplicate, labreplicate,
            analytename="TN",
            unit="mg/L",
            result=total_n_mgl
            ) %>%
  #Add Phosphorus
  bind_rows(
    nutrients_df %>%
      mutate(masterid = case_when(masterid=="\n204COY425"~"204COY425",
                                  T~masterid)) %>%
      filter(!is.na(total_p_mgl)) %>%
      transmute(masterid, 
                sampledate= as_date(sampledate),
                fieldreplicate, labreplicate,
                analytename="TP",
                unit="mg/L",
                result=total_p_mgl
      )
  ) %>%
  group_by(masterid, sampledate, fieldreplicate, analytename, unit  ) %>%
  summarise(result=mean(result, na.rm=T)) %>%
  ungroup( ) %>%
  group_by(masterid, sampledate,analytename, unit  ) %>%
  summarise(result=mean(result, na.rm=T)) %>%
  ungroup( )

phab_df2<-phab_df %>%
  filter(variable=="PCT_MAP",
         !is.na(result)) %>%
  inner_join(
    lustations_df %>% 
      select(stationcode=stationid, masterid)
  ) %>%
  transmute(
    masterid, 
    sampledate=as_date(sampledate),
    analytename="PCT_MAP",
    unit="Percent",
    result=result
  )

######
agg_bs_data_df<-
  bind_rows(chem_df2,
            nutrients_df2,
            phab_df2)
  
write_csv(agg_bs_data_df,"NotForGit/Step1/agg_bs_data_df.csv")
write_csv(algae_df2,"NotForGit/Step1/algae_df2.csv")
