#Prep data
# setwd("D:/Documents/SCCWRP/Biointegrity Biostim policy_RM/Data/RawData/Raw data 031617")
setwd("C:/Users/Raphaelm/Documents/SCCWRP/Biointegrity Biostim policy_RM/Data/RawData/Raw data 031617")
library(TITAN2)
library(ggplot2)
library(quantregForest)
library(dplyr)
library(plyr)
library(mgcv)
library(stringr)
library(sf)
library(lubridate)
library(tidyverse)
library(reshape2)
library(MASS)
library(lubridate)
# source("parse.R")

#Import the data
mydf<-read.csv("SMR/smc.swamp_mastertable_SMR_CRAM.csv", stringsAsFactors = F)
asci.df<-read.csv("SMR/asci.sites.forRafi_v3.csv", stringsAsFactors = F)
alg.df<-read.csv("SMR/algae.taxonomic.data.07252018.csv", stringsAsFactors = F)

thresholds.bs.df<-read.csv("R_Outputs_082318/tab.threshold.rr.summary.csv", stringsAsFactors = F) %>%
  dplyr::select(BSPretty, BIgoal, Prob, Response, Est, l95, u95, RR.l95.cal, RR.l95.val) %>%
  filter(RR.l95.cal>1 & RR.l95.val>1 & Prob=="p90")

thresholds.bs.df_sum<-ddply(thresholds.bs.df, .(BSPretty, BIgoal), summarize, MostConservative=min(Est))







# bugs.df<-ddply(bugs.df, .(StationCode,OTU), summarize, BAResult=sum(Iteration1)) #Work on this...

#MMI, if you prefer
asci.df$ASCI_D<-asci.df$MMI.diatoms 
asci.df$ASCI_S<-asci.df$MMI.sba
asci.df$ASCI_H<-asci.df$MMI.hybrid 

#Biostiulatory variables
chem.varz<-c("Nitrogen_Total_mgPerL","Phosphorus_as_P_mgPerL",
             # "Nitrate_as_N_mgPerL","Nitrite_as_N_mgPerL","Nitrate_Nitrite_as_N_mgPerL","OrthoPhosphate_as_P_mgPerL","Ammonia_as_N_mgPerL",
             NULL)
om.varz<-c(  "Chlorophyll_a_mgPerm2","Ash_Free_Dry_Mass_mgPercm2",
             "PCT_MAP",
             NULL)

bs.varz<-c(chem.varz, om.varz)
bs.varz_pretty<-c("Total N", "Total P", "Chl-a", "AFDM", "% cover")
# bs.varz_pretty<-c("Total N", "Total P","NO3-N", "NO2-N","NO3NO2-N","OrthoP","NH4-N", "Chl-a", "AFDM", "% cover")
bs.varz.df<-data.frame(BiostimVar=bs.varz, BSPretty=bs.varz_pretty, stringsAsFactors = F)

#Format dates and designate sampleIDs
mydf$Date2<-lubridate::mdy(mydf$SampleDate)
mydf$MasterDate<-paste(mydf$MasterID, mydf$Date2, sep="___")
mydf$MasterDate_rep<-paste(mydf$MasterID, mydf$Date2, mydf$Replicate ,sep="___")
#Define complete cases
mydf$Complete<-ifelse(!is.na(mydf$CSCI) &
                        !is.na(mydf$Nitrogen_Total_mgPerL) & !is.na(mydf$Phosphorus_as_P_mgPerL) & 
                        !is.na(mydf$Chlorophyll_a_mgPerm2) & !is.na(mydf$Ash_Free_Dry_Mass_mgPercm2),
                      TRUE, FALSE)
mydf.c<-mydf[which(mydf$Complete==T),]


#Select one sample per site, at random
sites.samps<-unique(mydf.c[,c("MasterID","MasterDate", "New_Lat","New_Long")])
samps.t<- sites.samps %>% group_by(MasterID)
set.seed(501)
samps.sel<-sample_n(samps.t, 1)
mydf.c$SelectedSample<-ifelse(mydf.c$MasterDate %in% samps.sel$MasterDate, "Selected","NotSelected")


#Repeat for ASCI
#Format dates and designate sampleIDs
asci.df$Date2<-lubridate::mdy(asci.df$SampleDate)
asci.df$MasterDate<-paste(asci.df$MasterID, asci.df$Date2, sep="___")
asci.df$MasterDate_rep<-paste(asci.df$MasterID, asci.df$Date2, asci.df$Replicate ,sep="___")
asci.df$PSA6c[which(asci.df$PSA6c=="NV")]<-"SN" #Correct input

alg.df$Date2<-lubridate::mdy(alg.df$SampleDate)
alg.df$MasterDate<-paste(alg.df$MasterID, alg.df$Date2, sep="___")
alg.df$MasterDate_rep<-paste(alg.df$MasterID, alg.df$Date2, alg.df$Replicate ,sep="___")


#Define complete cases
asci.df$Complete<-ifelse(!is.na(asci.df$ASCI_D) & !is.na(asci.df$ASCI_S) & !is.na(asci.df$ASCI_H) &
                           !is.na(asci.df$Nitrogen_Total_mgPerL) & !is.na(asci.df$Phosphorus_as_P_mgPerL) & 
                           !is.na(asci.df$Chlorophyll_a_mgPerm2) & !is.na(asci.df$Ash_Free_Dry_Mass_mgPercm2),
                         TRUE, FALSE)
asci.df.c<-asci.df[which(asci.df$Complete==T),]

#Select one sample per site, at random
sites.samps<-unique(asci.df.c[,c("MasterID","MasterDate", "New_Lat","New_Long")])
samps.t<- sites.samps %>% group_by(MasterID)
set.seed(601)
samps.sel<-sample_n(samps.t, 1)
asci.df.c$SelectedSample<-ifelse(asci.df.c$MasterDate %in% samps.sel$MasterDate, "Selected","NotSelected")




#Identify all possible samples and region
samples.df<-unique(data.frame(
  rbind(mydf.c[,c("MasterID","Date2","MasterDate","PSA6c", "New_Lat", "New_Long")],
        asci.df.c[,c("MasterID","Date2","MasterDate", "PSA6c", "New_Lat", "New_Long")])))

#Get in some site data from mydf, going to asci.df where necessary
#ugh
samples.df<-unique(data.frame(  rbind(mydf.c[,c("MasterID","Date2","MasterDate")],        asci.df.c[,c("MasterID","Date2","MasterDate")])))
samples.df<-join(samples.df, unique(mydf[,c("MasterID","PSA6c", "New_Lat", "New_Long")]))
samples.df$PSA6c<-  sapply(1:nrow(samples.df), function(i){
  psa.i<-samples.df$PSA6c[i]
  if(is.na(psa.i))
  {
    site.i<-samples.df$MasterID[i]
    unique(asci.df$PSA6c[which(asci.df$MasterID==site.i)])
  }
  else psa.i
})
samples.df$New_Lat<-  sapply(1:nrow(samples.df), function(i){
  psa.i<-samples.df$New_Lat[i]
  if(is.na(psa.i))
  {
    site.i<-samples.df$MasterID[i]
    unique(asci.df$New_Lat[which(asci.df$MasterID==site.i)])
  }
  else psa.i
})
samples.df$New_Long<-  sapply(1:nrow(samples.df), function(i){
  psa.i<-samples.df$New_Long[i]
  if(is.na(psa.i))
  {
    site.i<-samples.df$MasterID[i]
    unique(asci.df$New_Long[which(asci.df$MasterID==site.i)])
  }
  else psa.i
})


#Designate type of biointegrity data available
samples.df$BugSamples<-samples.df$MasterDate %in% mydf$MasterDate
samples.df$AlgaeSamples<-samples.df$MasterDate %in% asci.df$MasterDate
samples.df$BugAlgaeSamples<-samples.df$MasterDate %in% asci.df$MasterDate & samples.df$MasterDate %in% mydf$MasterDate
samples.df$SampleType<-ifelse(samples.df$BugSamples & samples.df$AlgaeSamples, "Both",
                              ifelse(samples.df$BugSamples & !samples.df$AlgaeSamples, "Bugs only",
                                     ifelse(!samples.df$BugSamples & samples.df$AlgaeSamples, "Algae only","Neither")))
table(samples.df$SampleType)
length(unique(samples.df$MasterID))
#Select 75% of each PSA region for calibration
#Stratify by data available.
#This way, the same sites will be selected for both algae and bug models, where possible

sites.psa<-ddply(samples.df, .(MasterID, PSA6c, New_Lat, New_Long), summarise, 
                 BugSamples=any(BugSamples), 
                 AlgaeSamples=any(AlgaeSamples),
                 BugAlgaeSamples=any(BugAlgaeSamples))
sites.psa$SampleType<-ifelse(sites.psa$BugSamples & sites.psa$AlgaeSamples, "Both",
                             ifelse(sites.psa$BugSamples & !sites.psa$AlgaeSamples, "Bugs only",
                                    ifelse(!sites.psa$BugSamples & sites.psa$AlgaeSamples, "Algae only","Neither")))
sites.t<-sites.psa %>% group_by(PSA6c, SampleType)
set.seed(502)
sites.sel<-sample_frac(sites.t, 0.75, replace=F)
sites.psa$DevSet<-ifelse(sites.psa$MasterID %in% sites.sel$MasterID,"Cal","Val")
mydf.c$DevSet <- ifelse(mydf.c$MasterID %in% sites.sel$MasterID,"Cal","Val")
asci.df.c$DevSet <- ifelse(asci.df.c$MasterID %in% sites.sel$MasterID,"Cal","Val")

addmargins(table(sites.psa$DevSet,sites.psa$SampleType))
addmargins(table(sites.psa$PSA6c,sites.psa$SampleType))
addmargins(table(mydf.c$PSA6c,mydf.c$DevSet))
addmargins(table(asci.df.c$PSA6c,asci.df.c$DevSet))
length(unique(mydf.c$MasterID))



# mydf2<-join(mydf, samples.df[,c("MasterDate","DevSet","SampleType")])

basemap.simple<-
  ggplot(data=  subset(map_data('state'), region=="california"), aes(x=long,y=lat))+geom_polygon(color=NA, fill="gray80")+
  theme_minimal()+theme(axis.text=element_blank(), axis.ticks=element_blank(), axis.title=element_blank())+
  coord_map()

dev.set.map<-basemap.simple+
  geom_point(data=sites.psa, aes(x=New_Long, y=New_Lat, fill=DevSet), shape=21, size=1)+
  geom_point(data=sites.psa[which(sites.psa$DevSet=="Val"),], aes(x=New_Long, y=New_Lat, fill=DevSet), shape=21, size=1)+
  scale_fill_manual(values=c("gray25","white"), name="Site type", labels=c("Calibration","Validation"))
# ggsave(dev.set.map, filename="R_Outputs_082318/map.jpg", dpi=300,height=2.5, width=3)

#Designate samples by biointegrity condition class

#MMI+OE thresholds:
# asci_d_thresh<-  qnorm(mean=1, sd=0.14, p=c(0.3,0.1,.01))
# asci_s_thresh<-  qnorm(mean=1, sd=0.16, p=c(0.3,0.1,.01))
# asci_h_thresh<-  qnorm(mean=1, sd=0.21, p=c(0.3,0.1,.01))

#MMI thresholds:
asci_d_thresh<-  qnorm(mean=1, sd=0.16, p=c(0.3,0.1,.01))
asci_s_thresh<-  qnorm(mean=1, sd=0.14, p=c(0.3,0.1,.01))
asci_h_thresh<-  qnorm(mean=1, sd=0.13, p=c(0.3,0.1,.01))

thresholds.df<-data.frame(Index=c("CSCI","ASCI_D","ASCI_S","ASCI_H"),
                          Ref30=c(0.92, asci_d_thresh[1], asci_s_thresh[1], asci_h_thresh[1]),
                          Ref10=c(0.79, asci_d_thresh[2], asci_s_thresh[2], asci_h_thresh[2]),
                          Ref01=c(0.63, asci_d_thresh[3], asci_s_thresh[3], asci_h_thresh[3]),
                          BCG2=c(1.025, 1.31, 1.36, 1.23),
                          BCG3=c(0.825, 0.95, 0.86, 0.97),
                          BCG4=c(0.625, 0.54, 0.36, 0.67),
                          BCG5=c(0.325, 0, 0, 0.30))

#TITAN
library(TITAN2)

alg.df2<-merge(alg.df, asci.df.c[,c("MasterID", "SampleID_old","DevSet","SelectedSample")])
setdiff(asci.df.c$SampleID_old,alg.df2$SampleID_old)
setdiff(alg.df2$SampleID_old,asci.df.c$SampleID_old)

#Make a site by gradient dataframe
alg_env.cal<-asci.df.c[which(asci.df.c$DevSet=="Cal" & asci.df.c$SelectedSample=="Selected"), c(chem.varz, om.varz)]
alg_map.cal<-na.omit(asci.df.c[which(asci.df.c$DevSet=="Cal" & asci.df.c$SelectedSample=="Selected"), c(chem.varz, om.varz)])
alg_map.cal.samps<-na.omit(asci.df.c[which(asci.df.c$DevSet=="Cal" & asci.df.c$SelectedSample=="Selected"), c(chem.varz, om.varz,"SampleID_old")])
str(alg_env.cal)




#Make a site by taxon dataframe with raw abundance
#Do diatoms alone based on abundance, SBA alone based on biovolume
alg.cal<-alg.df2[which(alg.df2$DevSet=="Cal" & alg.df2$SelectedSample=="Selected"),]
alg.diatom<-ddply(alg.cal[!is.na(alg.cal$Genus) & alg.cal$SampleTypeCode=="Integrated",], .(SampleID_old, Genus), summarize, BAResult=sum(BAResult, na.rm=T))
alg.sba<-ddply(alg.cal[!is.na(alg.cal$Genus) & alg.cal$SampleTypeCode!="Integrated",], .(SampleID_old, Genus), summarize, Result=sum(Result, na.rm=T))


alg.diatom.tot<-ddply(alg.cal[!is.na(alg.cal$Genus) & alg.cal$SampleTypeCode=="Integrated",], .(SampleID_old), summarize, Total_BAResult=sum(BAResult, na.rm=T))
alg.diatom<-join(alg.diatom, alg.diatom.tot)
alg.diatom$RelCount<-alg.diatom$BAResult/alg.diatom$Total_BAResult

alg.sba.tot<-ddply(alg.cal[!is.na(alg.cal$Genus) & alg.cal$SampleTypeCode!="Integrated",], .(SampleID_old), summarize, Total_Result=sum(Result, na.rm=T))
alg.sba<-join(alg.sba, alg.sba.tot)
alg.sba$RelBiovol<-alg.sba$Result/alg.sba$Total_Result


diatom.freq.cal<-ddply(alg.diatom, .(Genus), summarize, freq=sum(BAResult>0))
sba.freq.cal<-ddply(alg.sba, .(Genus), summarize, freq=sum(Result>0))


min.freq<-length(unique(alg.cal$SampleID_old))*.05
diatom.freq.cal.min<-diatom.freq.cal$Genus[diatom.freq.cal$freq>=min.freq]
sba.freq.cal.min<-sba.freq.cal$Genus[sba.freq.cal$freq>=min.freq]




diatom_tax.cal<-dcast(alg.diatom[which(alg.diatom$Genus %in% diatom.freq.cal.min) ,], 
                      SampleID_old~Genus, 
                      value.var="RelCount", fill = 0)
diatom_tax.cal<-diatom_tax.cal[,diatom.freq.cal.min]
diatom.map.cal<-dcast(alg.diatom[which(alg.diatom$Genus %in% diatom.freq.cal.min & alg.diatom$SampleID_old %in% alg_map.cal.samps$SampleID_old ) ,], 
                      SampleID_old~Genus, 
                      value.var="RelCount", fill = 0)
diatom.map.cal<-diatom.map.cal[,diatom.freq.cal.min]

sba_tax.cal<-dcast(alg.sba[which(alg.sba$Genus %in% sba.freq.cal.min) ,], 
                   SampleID_old~Genus, 
                   value.var="RelBiovol", fill = 0)
sba_tax.cal<-sba_tax.cal[,sba.freq.cal.min]
sba.map.cal<-dcast(alg.sba[which(alg.sba$Genus %in% sba.freq.cal.min & alg.sba$SampleID_old %in% alg_map.cal.samps$SampleID_old ) ,], 
                   SampleID_old~Genus, 
                   value.var="RelBiovol", fill = 0)
sba.map.cal<-sba.map.cal[,sba.freq.cal.min]



#Default settings: nBoot = 500, numPerm = 250. Set both to 10 for exploration.
nb.i<-1000
perm.i<-500

set.seed(2)
diatom_tn_cal<-titan(alg_env.cal$Nitrogen_Total_mgPerL, diatom_tax.cal, nBoot = nb.i, numPerm = perm.i)
set.seed(3)
diatom_tp_cal<-titan(alg_env.cal$Phosphorus_as_P_mgPerL, diatom_tax.cal, nBoot = nb.i, numPerm = perm.i)
set.seed(4)
diatom_chla_cal<-titan(alg_env.cal$Chlorophyll_a_mgPerm2, diatom_tax.cal, nBoot = nb.i, numPerm = perm.i)
set.seed(5)
diatom_afdm_cal<-titan(alg_env.cal$Ash_Free_Dry_Mass_mgPercm2, diatom_tax.cal, nBoot = nb.i, numPerm = perm.i)
set.seed(6)
diatom_map_cal<-titan(alg_map.cal$PCT_MAP, diatom.map.cal, nBoot = nb.i, numPerm = perm.i)

set.seed(12)
sba_tn_cal<-titan(alg_env.cal$Nitrogen_Total_mgPerL, sba_tax.cal, nBoot = nb.i, numPerm = perm.i)
save.image("R_Outputs_082318/titan_analyses_082918.Rdata")
set.seed(13)
sba_tp_cal<-titan(alg_env.cal$Phosphorus_as_P_mgPerL, sba_tax.cal, nBoot = nb.i, numPerm = perm.i)
save.image("R_Outputs_082318/titan_analyses_082918.Rdata")
set.seed(14)
sba_chla_cal<-titan(alg_env.cal$Chlorophyll_a_mgPerm2, sba_tax.cal, nBoot = nb.i, numPerm = perm.i)
save.image("R_Outputs_082318/titan_analyses_082918.Rdata")
set.seed(15)
sba_afdm_cal<-titan(alg_env.cal$Ash_Free_Dry_Mass_mgPercm2, sba_tax.cal, nBoot = nb.i, numPerm = perm.i)
save.image("R_Outputs_082318/titan_analyses_082918.Rdata")
#START HERE
set.seed(16)
sba_map_cal<-titan(alg_map.cal$PCT_MAP, sba.map.cal, nBoot = nb.i, numPerm = perm.i)
save.image("R_Outputs_082318/titan_analyses_082918.Rdata")

####
#BUGS
# mydf.c$CollectionMethodCode<-mydf.c$CollectionMethodCode_csci
mydf.cA<-mydf.c
mydf.cA$SampleIDx<-
  paste(
    substr(mydf.cA$SampleID,1, nchar(mydf.cA$SampleID)-2),
    mydf.cA$CollectionMethodCode_csci,
    mydf.cA$Replicate,
    sep="_"
  )

# bugs.df<-join(read.csv("SMR/Bug_OE.csv", stringsAsFactors = F),
#               read.csv("SMR/bug.samples_parsed.csv", stringsAsFactors = F))
bugs.df<-read.csv("SMR/Bug_OE.csv", stringsAsFactors = F)
bugs.df<-bugs.df[which(!bugs.df$OTU %in% c("OoverE")),]
bugs.df$BAResult<-bugs.df$Iteration1
bugs.df$SampleIDx<-bugs.df$SampleID
bugs.df$SampleIDx<-gsub("801S02167","SMC02167",bugs.df$SampleIDx)
bugs.df$SampleIDx<-gsub("801S02573","SMC02573",bugs.df$SampleIDx)
bugs.df$SampleIDx<-gsub("801S03533","SMC03533",bugs.df$SampleIDx)
bugs.df$SampleIDx<-gsub("801S09591","SMC09591",bugs.df$SampleIDx)
bugs.df$SampleIDx<-gsub("801S00135","SMC00135",bugs.df$SampleIDx)
bugs.df$SampleIDx<-gsub("801S00375","SMC00375",bugs.df$SampleIDx)
bugs.df$SampleIDx<-gsub("801S01383","SMC01383",bugs.df$SampleIDx)
bugs.df$SampleIDx<-gsub("801S02059","SMC02059",bugs.df$SampleIDx)
bugs.df$SampleIDx<-gsub("801S03133","SMC03133",bugs.df$SampleIDx)
bugs.df$SampleIDx<-gsub("801S03687","SMC03687",bugs.df$SampleIDx)
bugs.df$SampleIDx<-gsub("802S09698","SMC09698",bugs.df$SampleIDx)
bugs.df$SampleIDx<-gsub("SMC09698_04262016_BMI_RWB_1","SMC09698_04212015_BMI_RWB_1",bugs.df$SampleIDx) #Date error?
bugs.df$SampleIDx<-gsub("802S27709","SMC27709",bugs.df$SampleIDx)

bugs.df2<-merge(bugs.df, mydf.cA[,c("SampleIDx","DevSet","SelectedSample")])
# intersect(bugs.df$SampleIDx, mydf.cA$SampleIDx)
setdiff(mydf.cA$SampleIDx, bugs.df2$SampleIDx)

#Make a site by gradient dataframe


bug_env.cal<-mydf.cA[which(mydf.cA$DevSet=="Cal" & mydf.cA$SelectedSample=="Selected"), c(chem.varz, om.varz)] #There should be no supes among selected SampleIDx
bug_map.cal<-na.omit(mydf.cA[which(mydf.c$DevSet=="Cal" & mydf.cA$SelectedSample=="Selected"), c(chem.varz, om.varz)])
bug_map.cal.samps<-na.omit(mydf.cA[which(mydf.cA$DevSet=="Cal" & mydf.cA$SelectedSample=="Selected"), c("SampleIDx", chem.varz, om.varz)])
summary(bug_map.cal.samps)

#Make a site by taxon dataframe with raw abundance
head(bugs.df)
bug.cal<-bugs.df2[which(bugs.df2$DevSet=="Cal" & bugs.df2$SelectedSample=="Selected"),]
# bug.cal<-ddply(bug.cal, .(SampleIDx, OTU), summarise, BAResult=sum(BAResult, na.rm=T))

bugs.tot<-ddply(bug.cal[!is.na(bug.cal$OTU),], .(SampleIDx), summarize, Total_BAResult=sum(BAResult, na.rm=T))
bug.cal<-join(bug.cal, bugs.tot)
bug.cal$RelCount<-bug.cal$BAResult/bug.cal$Total_BAResult

bugs.freq.cal<-ddply(bug.cal, .(OTU), summarize, freq=sum(BAResult>0))
min.freq_bugs<-length(unique(bug.cal$StationCode))*.05
bugs.freq.cal.min<-bugs.freq.cal$OTU[bugs.freq.cal$freq>=min.freq_bugs]

bug_tax.cal<-dcast(bug.cal[which(bug.cal$OTU %in% bugs.freq.cal.min) ,], 
                   SampleIDx~OTU, #fun.aggregate = mean, #there are some collection method replicates that are duplicating records
                   value.var="RelCount", fill = 0)

bug_tax.cal<-bug_tax.cal[,bugs.freq.cal.min]
bug.map.cal<-dcast(bug.cal[which(bug.cal$OTU %in% bugs.freq.cal.min & bug.cal$SampleIDx %in% bug_map.cal.samps$SampleIDx ) ,], 
                   SampleIDx~OTU, #fun.aggregate = mean, #there are some collection method replicates that are duplicating records
                   value.var="RelCount", fill = 0)
bug.map.cal<-bug.map.cal[,intersect(names(bug.map.cal), bugs.freq.cal.min)]



str(bug_tax.cal)
###

set.seed(120)
bugs_tn_cal<-titan(bug_env.cal$Nitrogen_Total_mgPerL, bug_tax.cal)
save.image("R_Outputs_082318/titan_analyses_082918.Rdata")
set.seed(130)
bugs_tp_cal<-titan(bug_env.cal$Phosphorus_as_P_mgPerL, bug_tax.cal)
save.image("R_Outputs_082318/titan_analyses_082918.Rdata")
set.seed(140)
bugs_chla_cal<-titan(bug_env.cal$Chlorophyll_a_mgPerm2, bug_tax.cal, nBoot=100, numPerm=100)
save.image("R_Outputs_082318/titan_analyses_082918.Rdata")
#Redo from here (set to 100/100 instead of defaults)
set.seed(150)
bugs_afdm_cal<-titan(bug_env.cal$Ash_Free_Dry_Mass_mgPercm2, bug_tax.cal, nBoot=100, numPerm=100)
save.image("R_Outputs_082318/titan_analyses_082918.Rdata")
set.seed(160)
bugs_map_cal<-titan(bug_map.cal$PCT_MAP, bug.map.cal, nBoot=100, numPerm=100)
save.image("R_Outputs_082318/titan_analyses_082918.Rdata")

# save.image("R_Outputs_082318/titan_analyses_031319.Rdata")
# load("D:/Documents/SCCWRP/Biointegrity Biostim policy_RM/Data/RawData/Raw data 031617/R_Outputs_082318/titan_analyses_031319.Rdata")
load("C:/Users/Raphaelm/Documents/SCCWRP/Biointegrity Biostim policy_RM/Data/RawData/Raw data 031617/R_Outputs_082318/titan_analyses_031319.Rdata")

head(diatom_tn_cal$sppmax)
head(diatom_tn_cal$sumz.cp)
head(diatom_tn_cal$maxSumz)


sba_chla_cal$maxSumz
sba_tn_cal$sumz.cp
plotSumz(diatom_tn_cal, filter = T, xlabel="Total N (mg/L)")
plotTaxa(diatom_tn_cal, filter=T, xlabel="Total N (mg/L)")

head(diatom_tn_cal$sppmax)
diatom_tn_cal$sppmax["Denticula",]

#Make plot data
#tn
dia.plot.tn.dat<-as.data.frame(diatom_tn_cal$sppmax)
dia.plot.tn.dat$Taxon<-row.names(dia.plot.tn.dat)
dia.plot.tn.dat$Assemblage<-"Diatom"
sba.plot.tn.dat<-as.data.frame(sba_tn_cal$sppmax)
sba.plot.tn.dat$Taxon<-row.names(sba.plot.tn.dat)
sba.plot.tn.dat$Assemblage<-"SBA"
bug.plot.tn.dat<-as.data.frame(bugs_tn_cal$sppmax)
bug.plot.tn.dat$Taxon<-row.names(bug.plot.tn.dat)
bug.plot.tn.dat$Assemblage<-"BMI"

plot.tn.dat<-rbind(dia.plot.tn.dat, sba.plot.tn.dat,bug.plot.tn.dat)
plot.tn.dat$Group<-factor(plot.tn.dat$filter, levels=c(0,2,1), labels = c("Neither","Increasers","Decreasers"))
plot.tn.dat$Taxon<-factor(plot.tn.dat$Taxon, levels=plot.tn.dat$Taxon[order(plot.tn.dat$zenv.cp, decreasing=T)])
plot.tn.sumdat<-data.frame(Assemblage=c("Diatom","Diatom","SBA","SBA","BMI","BMI"),
                           Group=c("Increasers","Decreasers","Increasers","Decreasers","Increasers","Decreasers"), 
                           cp=c(diatom_tn_cal$sumz.cp[4,1], diatom_tn_cal$sumz.cp[3,1],
                                sba_tn_cal$sumz.cp[4,1], sba_tn_cal$sumz.cp[3,1],
                                bugs_tn_cal$sumz.cp[4,1], bugs_tn_cal$sumz.cp[3,1]),
                           cp.05=c(diatom_tn_cal$sumz.cp[4,2], diatom_tn_cal$sumz.cp[3,2],
                                   sba_tn_cal$sumz.cp[4,2], sba_tn_cal$sumz.cp[3,2],
                                   bugs_tn_cal$sumz.cp[4,2], bugs_tn_cal$sumz.cp[3,2]),
                           cp.95=c(diatom_tn_cal$sumz.cp[4,6], diatom_tn_cal$sumz.cp[3,6],
                                   sba_tn_cal$sumz.cp[4,6], sba_tn_cal$sumz.cp[3,6],
                                   bugs_tn_cal$sumz.cp[4,6], bugs_tn_cal$sumz.cp[3,6]),
                           stringsAsFactors = F)
plot.tn.dat$z2<-ifelse(plot.tn.dat$Group=="Decreasers", -plot.tn.dat$zscore,plot.tn.dat$zscore)





#tp
dia.plot.tp.dat<-as.data.frame(diatom_tp_cal$sppmax)
dia.plot.tp.dat$Taxon<-row.names(dia.plot.tp.dat)
dia.plot.tp.dat$Assemblage<-"Diatom"
sba.plot.tp.dat<-as.data.frame(sba_tp_cal$sppmax)
sba.plot.tp.dat$Taxon<-row.names(sba.plot.tp.dat)
sba.plot.tp.dat$Assemblage<-"SBA"
bug.plot.tp.dat<-as.data.frame(bugs_tp_cal$sppmax)
bug.plot.tp.dat$Taxon<-row.names(bug.plot.tp.dat)
bug.plot.tp.dat$Assemblage<-"BMI"

plot.tp.dat<-rbind(dia.plot.tp.dat, sba.plot.tp.dat,bug.plot.tp.dat)
plot.tp.dat$Group<-factor(plot.tp.dat$filter, levels=c(0,2,1), labels = c("Neither","Increasers","Decreasers"))
plot.tp.dat$Taxon<-factor(plot.tp.dat$Taxon, levels=plot.tp.dat$Taxon[order(plot.tp.dat$zenv.cp, decreasing=T)])
plot.tp.sumdat<-data.frame(Assemblage=c("Diatom","Diatom","SBA","SBA","BMI","BMI"),
                           Group=c("Increasers","Decreasers","Increasers","Decreasers","Increasers","Decreasers"), 
                           cp=c(diatom_tp_cal$sumz.cp[4,1], diatom_tp_cal$sumz.cp[3,1],
                                sba_tp_cal$sumz.cp[4,1], sba_tp_cal$sumz.cp[3,1],
                                bugs_tp_cal$sumz.cp[4,1], bugs_tp_cal$sumz.cp[3,1]),
                           cp.05=c(diatom_tp_cal$sumz.cp[4,2], diatom_tp_cal$sumz.cp[3,2],
                                   sba_tp_cal$sumz.cp[4,2], sba_tp_cal$sumz.cp[3,2],
                                   bugs_tp_cal$sumz.cp[4,2], bugs_tp_cal$sumz.cp[3,2]),
                           cp.95=c(diatom_tp_cal$sumz.cp[4,6], diatom_tp_cal$sumz.cp[3,6],
                                   sba_tp_cal$sumz.cp[4,6], sba_tp_cal$sumz.cp[3,6],
                                   bugs_tp_cal$sumz.cp[4,6], bugs_tp_cal$sumz.cp[3,6]),
                           stringsAsFactors = F)
plot.tp.dat$z2<-ifelse(plot.tp.dat$Group=="Decreasers", -plot.tp.dat$zscore,plot.tp.dat$zscore)


plot.tp.dat %>%
  group_by(Assemblage) %>%
  filter(Group!="Neither") %>%
  filter(IndVal == max(IndVal) ) %>%
  data.frame()


#chla
dia.plot.chla.dat<-as.data.frame(diatom_chla_cal$sppmax)
dia.plot.chla.dat$Taxon<-row.names(dia.plot.chla.dat)
dia.plot.chla.dat$Assemblage<-"Diatom"
sba.plot.chla.dat<-as.data.frame(sba_chla_cal$sppmax)
sba.plot.chla.dat$Taxon<-row.names(sba.plot.chla.dat)
sba.plot.chla.dat$Assemblage<-"SBA"
bug.plot.chla.dat<-as.data.frame(bugs_chla_cal$sppmax)
bug.plot.chla.dat$Taxon<-row.names(bug.plot.chla.dat)
bug.plot.chla.dat$Assemblage<-"BMI"

plot.chla.dat<-rbind(dia.plot.chla.dat, sba.plot.chla.dat,bug.plot.chla.dat)
plot.chla.dat$Group<-factor(plot.chla.dat$filter, levels=c(0,2,1), labels = c("Neither","Increasers","Decreasers"))
plot.chla.dat$Taxon<-factor(plot.chla.dat$Taxon, levels=plot.chla.dat$Taxon[order(plot.chla.dat$zenv.cp, decreasing=T)])
plot.chla.sumdat<-data.frame(Assemblage=c("Diatom","Diatom","SBA","SBA","BMI","BMI"),
                             Group=c("Increasers","Decreasers","Increasers","Decreasers","Increasers","Decreasers"), 
                             cp=c(diatom_chla_cal$sumz.cp[4,1], diatom_chla_cal$sumz.cp[3,1],
                                  sba_chla_cal$sumz.cp[4,1], sba_chla_cal$sumz.cp[3,1],
                                  bugs_chla_cal$sumz.cp[4,1], bugs_chla_cal$sumz.cp[3,1]),
                             cp.05=c(diatom_chla_cal$sumz.cp[4,2], diatom_chla_cal$sumz.cp[3,2],
                                     sba_chla_cal$sumz.cp[4,2], sba_chla_cal$sumz.cp[3,2],
                                     bugs_chla_cal$sumz.cp[4,2], bugs_chla_cal$sumz.cp[3,2]),
                             cp.95=c(diatom_chla_cal$sumz.cp[4,6], diatom_chla_cal$sumz.cp[3,6],
                                     sba_chla_cal$sumz.cp[4,6], sba_chla_cal$sumz.cp[3,6],
                                     bugs_chla_cal$sumz.cp[4,6], bugs_chla_cal$sumz.cp[3,6]),
                             stringsAsFactors = F)
plot.chla.dat$z2<-ifelse(plot.chla.dat$Group=="Decreasers", -plot.chla.dat$zscore,plot.chla.dat$zscore)

#afdm
dia.plot.afdm.dat<-as.data.frame(diatom_afdm_cal$sppmax)
dia.plot.afdm.dat$Taxon<-row.names(dia.plot.afdm.dat)
dia.plot.afdm.dat$Assemblage<-"Diatom"
sba.plot.afdm.dat<-as.data.frame(sba_afdm_cal$sppmax)
sba.plot.afdm.dat$Taxon<-row.names(sba.plot.afdm.dat)
sba.plot.afdm.dat$Assemblage<-"SBA"
bug.plot.afdm.dat<-as.data.frame(bugs_afdm_cal$sppmax)
bug.plot.afdm.dat$Taxon<-row.names(bug.plot.afdm.dat)
bug.plot.afdm.dat$Assemblage<-"BMI"

plot.afdm.dat<-rbind(dia.plot.afdm.dat, sba.plot.afdm.dat,bug.plot.afdm.dat)
plot.afdm.dat$Group<-factor(plot.afdm.dat$filter, levels=c(0,2,1), labels = c("Neither","Increasers","Decreasers"))
plot.afdm.dat$Taxon<-factor(plot.afdm.dat$Taxon, levels=plot.afdm.dat$Taxon[order(plot.afdm.dat$zenv.cp, decreasing=T)])
plot.afdm.sumdat<-data.frame(Assemblage=c("Diatom","Diatom","SBA","SBA","BMI","BMI"),
                             Group=c("Increasers","Decreasers","Increasers","Decreasers","Increasers","Decreasers"), 
                             cp=c(diatom_afdm_cal$sumz.cp[4,1], diatom_afdm_cal$sumz.cp[3,1],
                                  sba_afdm_cal$sumz.cp[4,1], sba_afdm_cal$sumz.cp[3,1],
                                  bugs_afdm_cal$sumz.cp[4,1], bugs_afdm_cal$sumz.cp[3,1]),
                             cp.05=c(diatom_afdm_cal$sumz.cp[4,2], diatom_afdm_cal$sumz.cp[3,2],
                                     sba_afdm_cal$sumz.cp[4,2], sba_afdm_cal$sumz.cp[3,2],
                                     bugs_afdm_cal$sumz.cp[4,2], bugs_afdm_cal$sumz.cp[3,2]),
                             cp.95=c(diatom_afdm_cal$sumz.cp[4,6], diatom_afdm_cal$sumz.cp[3,6],
                                     sba_afdm_cal$sumz.cp[4,6], sba_afdm_cal$sumz.cp[3,6],
                                     bugs_afdm_cal$sumz.cp[4,6], bugs_afdm_cal$sumz.cp[3,6]),
                             stringsAsFactors = F)
plot.afdm.dat$z2<-ifelse(plot.afdm.dat$Group=="Decreasers", -plot.afdm.dat$zscore,plot.afdm.dat$zscore)

#map
dia.plot.map.dat<-as.data.frame(diatom_map_cal$sppmax)
dia.plot.map.dat$Taxon<-row.names(dia.plot.map.dat)
dia.plot.map.dat$Assemblage<-"Diatom"
sba.plot.map.dat<-as.data.frame(sba_map_cal$sppmax)
sba.plot.map.dat$Taxon<-row.names(sba.plot.map.dat)
sba.plot.map.dat$Assemblage<-"SBA"
bug.plot.map.dat<-as.data.frame(bugs_map_cal$sppmax)
bug.plot.map.dat$Taxon<-row.names(bug.plot.map.dat)
bug.plot.map.dat$Assemblage<-"BMI"

plot.map.dat<-rbind(dia.plot.map.dat, sba.plot.map.dat,bug.plot.map.dat)
plot.map.dat$Group<-factor(plot.map.dat$filter, levels=c(0,2,1), labels = c("Neither","Increasers","Decreasers"))
plot.map.dat$Taxon<-factor(plot.map.dat$Taxon, levels=plot.map.dat$Taxon[order(plot.map.dat$zenv.cp, decreasing=T)])
plot.map.sumdat<-data.frame(Assemblage=c("Diatom","Diatom","SBA","SBA","BMI","BMI"),
                            Group=c("Increasers","Decreasers","Increasers","Decreasers","Increasers","Decreasers"), 
                            cp=c(diatom_map_cal$sumz.cp[4,1], diatom_map_cal$sumz.cp[3,1],
                                 sba_map_cal$sumz.cp[4,1], sba_map_cal$sumz.cp[3,1],
                                 bugs_map_cal$sumz.cp[4,1], bugs_map_cal$sumz.cp[3,1]),
                            cp.05=c(diatom_map_cal$sumz.cp[4,2], diatom_map_cal$sumz.cp[3,2],
                                    sba_map_cal$sumz.cp[4,2], sba_map_cal$sumz.cp[3,2],
                                    bugs_map_cal$sumz.cp[4,2], bugs_map_cal$sumz.cp[3,2]),
                            cp.95=c(diatom_map_cal$sumz.cp[4,6], diatom_map_cal$sumz.cp[3,6],
                                    sba_map_cal$sumz.cp[4,6], sba_map_cal$sumz.cp[3,6],
                                    bugs_map_cal$sumz.cp[4,6], bugs_map_cal$sumz.cp[3,6]),
                            stringsAsFactors = F)
plot.map.dat$z2<-ifelse(plot.map.dat$Group=="Decreasers", -plot.map.dat$zscore,plot.map.dat$zscore)



min.z<-6
head(plot.tp.dat)
tn.taxa*.07
tn.taxa<-length(unique(droplevels(plot.tn.dat[which(plot.tn.dat$filter>0 & plot.tn.dat$zscore>min.z),"Taxon"])))
tp.taxa<-length(unique(droplevels(plot.tp.dat[which(plot.tp.dat$filter>0 & plot.tp.dat$zscore>min.z),"Taxon"])))
chla.taxa<-length(unique(droplevels(plot.chla.dat[which(plot.chla.dat$filter>0 & plot.chla.dat$zscore>min.z),"Taxon"])))
afdm.taxa<-length(unique(droplevels(plot.afdm.dat[which(plot.afdm.dat$filter>0 & plot.afdm.dat$zscore>min.z),"Taxon"])))
map.taxa<-length(unique(droplevels(plot.map.dat[which(plot.map.dat$filter>0 & plot.map.dat$zscore>min.z),"Taxon"])))

tn_cp_df<-data.frame(Assemblage=c("Diatom","Diatom","SBA","SBA","BMI","BMI"),
                     Response=c("Increasers","Decreasers"),
                     zenv.cp=c(.44,.38,.58,.17,.65,.65))

tn.titan.plot<-
  ggplot(data=droplevels(plot.tn.dat[which(plot.tn.dat$filter>0 & plot.tn.dat$zscore>min.z),]), aes(x=Taxon, y=zenv.cp))+
  # geom_rect(data=plot.tn.sumdat, aes(ymin=cp.05, ymax=cp.95, xmin=-Inf, xmax=Inf, x=NULL, y=NULL), fill="gray",color=NA, alpha=0.35)+
  geom_hline(data=plot.tn.sumdat, aes(yintercept=cp), linetype="dashed", color="white")+
  geom_linerange(aes(ymin=`5%`, ymax=`95%`), color="gray")+
  # geom_point(aes(size=zscore, fill=Assemblage, shape=Group))+
  geom_point(aes(fill=Assemblage, shape=Group), size=1.5)+
  # geom_hline(data=thresholds.bs.df_sum[which(thresholds.bs.df_sum$BSPretty=="Total N" & thresholds.bs.df_sum$BIgoal %in% c("BCG3","BCG4")),], aes(yintercept=MostConservative), linetype="dashed", color="#fed9a6", size=1)+
  # geom_hline(data=thresholds.bs.df_sum[which(thresholds.bs.df_sum$BSPretty=="Total N" & thresholds.bs.df_sum$BIgoal %in% c("Ref30","Ref10","Ref01")),], aes(yintercept=MostConservative), linetype="dotted", color="black", size=0.5)+
  geom_hline(data=tn_cp_df, aes(yintercept=zenv.cp, color=Assemblage, linetype=Response), size=0.5)+scale_linetype_manual(values=c("dotted","dashed"), name="Mean response")+
  scale_shape_manual(values=c(24,25), name="Response")+
  # scale_fill_manual(values=c("#1f78b4","#b2df8a","#33a02c"))+
  scale_fill_manual(values=c("skyblue","tomato1","darkgreen"))+
  # facet_wrap(~Group, scales="free")+
  scale_y_continuous(trans="log10", breaks=c(.1,.5,1,5,10))+
  coord_flip()+
  theme_bw(base_size=10)+
  xlab("")+ylab("Changepoint")+
  ggtitle("Total N (mg/L)")+
  scale_size_continuous("Indicator value\n(z-score)")+
  theme(axis.text.y = element_text(face = "italic", size=6), panel.grid = element_blank())+
  guides(fill = guide_legend(override.aes=list(shape=22, color=NA, size=5)),
         shape = guide_legend(override.aes=list(size=2)))
tn.titan.plot
ggsave(tn.titan.plot, filename="Figures_stratified/TITAN/tn.titan.plot_assemblage_mean.jpg", dpi=300, width=6, height=tn.taxa*.1)

tn.titan.plot2<-
  ggplot(data=droplevels(plot.tn.dat[which(plot.tn.dat$filter>0 & plot.tn.dat$zscore>min.z),]), aes(x=Taxon, y=zenv.cp))+
  # geom_rect(data=plot.tn.sumdat, aes(ymin=cp.05, ymax=cp.95, xmin=-Inf, xmax=Inf, x=NULL, y=NULL), fill="gray",color=NA, alpha=0.35)+
  geom_hline(data=plot.tn.sumdat, aes(yintercept=cp), linetype="dashed", color="white")+
  geom_linerange(aes(ymin=`5%`, ymax=`95%`), color="gray")+
  # geom_point(aes(size=zscore, fill=Assemblage, shape=Group))+
  geom_point(aes(fill=Assemblage, shape=Group), size=1.5)+
  # geom_hline(data=thresholds.bs.df_sum[which(thresholds.bs.df_sum$BSPretty=="Total N" & thresholds.bs.df_sum$BIgoal %in% c("BCG3","BCG4")),], aes(yintercept=MostConservative), linetype="dashed", color="#fed9a6", size=1)+
  # geom_hline(data=thresholds.bs.df_sum[which(thresholds.bs.df_sum$BSPretty=="Total N" & thresholds.bs.df_sum$BIgoal %in% c("Ref30","Ref10","Ref01")),], aes(yintercept=MostConservative), linetype="dotted", color="black", size=0.5)+
  # geom_hline(data=tn_cp_df, aes(yintercept=zenv.cp, color=Assemblage, linetype=Response), size=0.5)+scale_linetype_manual(values=c("dotted","dashed"), name="Mean response")+
  geom_hline(yintercept=c(0.5375, 0.21,0.129,0.069), linetype="dashed")+
  scale_shape_manual(values=c(24,25), name="Response")+
  # scale_fill_manual(values=c("#1f78b4","#b2df8a","#33a02c"))+
  scale_fill_manual(values=c("skyblue","tomato1","darkgreen"))+
  # facet_wrap(~Group, scales="free")+
  scale_y_continuous(trans="log10", breaks=c(.1,.5,1,5,10))+
  coord_flip()+
  theme_bw(base_size=10)+
  xlab("")+ylab("Changepoint")+
  ggtitle("Total N (mg/L)")+
  scale_size_continuous("Indicator value\n(z-score)")+
  theme(axis.text.y = element_text(face = "italic", size=6), panel.grid = element_blank())+
  guides(fill = guide_legend(override.aes=list(shape=22, color=NA, size=5)),
         shape = guide_legend(override.aes=list(size=2)))
tn.titan.plot2
ggsave(tn.titan.plot2, filename="Figures_stratified/TITAN/tn.titan.plot_ref10_probs.jpg", dpi=300, width=6, height=tn.taxa*.1)




tn.ref10<-thresholds.bs.df_sum[which(thresholds.bs.df_sum$BSPretty=="Total N" &
                                       thresholds.bs.df_sum$BIgoal %in% c("Ref01")),"MostConservative"]
plot.tn.dat %>%
  group_by(Assemblage) %>%
  filter(Group!="Neither" & zenv.cp>=tn.ref10) %>%
  ddply(.(Assemblage),summarise, ntaxa=length(Taxon))

tp.ref10<-thresholds.bs.df_sum[which(thresholds.bs.df_sum$BSPretty=="Total P" &
                                       thresholds.bs.df_sum$BIgoal %in% c("Ref01")),"MostConservative"]
plot.tp.dat %>%
  group_by(Assemblage) %>%
  filter(Group!="Neither" & zenv.cp>=tp.ref10) %>%
  ddply(.(Assemblage),summarise, ntaxa=length(Taxon))

chla.ref10<-thresholds.bs.df_sum[which(thresholds.bs.df_sum$BSPretty=="Chl-a" &
                                         thresholds.bs.df_sum$BIgoal %in% c("Ref01")),"MostConservative"]
plot.chla.dat %>%
  group_by(Assemblage) %>%
  filter(Group!="Neither" & zenv.cp>=chla.ref10) %>%
  ddply(.(Assemblage),summarise, ntaxa=length(Taxon))

afdm.ref10<-thresholds.bs.df_sum[which(thresholds.bs.df_sum$BSPretty=="AFDM" &
                                         thresholds.bs.df_sum$BIgoal %in% c("Ref01")),"MostConservative"]
plot.afdm.dat %>%
  group_by(Assemblage) %>%
  filter(Group!="Neither" & zenv.cp>=afdm.ref10) %>%
  ddply(.(Assemblage),summarise, ntaxa=length(Taxon))

map.ref10<-thresholds.bs.df_sum[which(thresholds.bs.df_sum$BSPretty=="% cover" &
                                        thresholds.bs.df_sum$BIgoal %in% c("Ref01")),"MostConservative"]
plot.map.dat %>%
  group_by(Assemblage) %>%
  filter(Group!="Neither" & zenv.cp>=map.ref10) %>%
  ddply(.(Assemblage),summarise, ntaxa=length(Taxon))


tp_cp_df<-data.frame(Assemblage=c("Diatom","Diatom","SBA","SBA","BMI","BMI"),
                     Response=c("Increasers","Decreasers"),
                     zenv.cp=c(.082,.048,.075,.034,.091,.080))

tp.titan.plot<-
  ggplot(data=droplevels(plot.tp.dat[which(plot.tp.dat$filter>0 & plot.tp.dat$zscore>min.z),]), aes(x=Taxon, y=zenv.cp))+
  # geom_rect(data=plot.tp.sumdat, aes(ymin=cp.05, ymax=cp.95, xmin=-Inf, xmax=Inf, x=NULL, y=NULL), fill="gray",color=NA, alpha=0.35)+
  geom_hline(data=plot.tp.sumdat, aes(yintercept=cp), linetype="dashed", color="white")+
  geom_linerange(aes(ymin=`5%`, ymax=`95%`), color="gray")+
  # geom_point(aes(size=zscore, fill=Assemblage, shape=Group))+
  geom_point(aes(fill=Assemblage, shape=Group), size=1.5)+
  # geom_hline(data=thresholds.bs.df_sum[which(thresholds.bs.df_sum$BSPretty=="Total P" & thresholds.bs.df_sum$BIgoal %in% c("BCG3","BCG4")),], aes(yintercept=MostConservative), linetype="dashed", color="#fed9a6", size=1)+
  # geom_hline(data=thresholds.bs.df_sum[which(thresholds.bs.df_sum$BSPretty=="Total P" & thresholds.bs.df_sum$BIgoal %in% c("Ref30","Ref10","Ref01")),], aes(yintercept=MostConservative), linetype="dotted", color="black", size=0.5)+
  geom_hline(data=tp_cp_df, aes(yintercept=zenv.cp, color=Assemblage, linetype=Response), size=0.5)+scale_linetype_manual(values=c("dotted","dashed"), name="Mean response")+
  scale_shape_manual(values=c(24,25), name="Response")+
  scale_fill_manual(values=c("skyblue","tomato1","darkgreen"))+
  # facet_wrap(~Group, scales="free")+
  scale_y_continuous(trans="log10", breaks=c(0.01,0.05,0.1,0.5,1))+
  coord_flip()+
  theme_bw()+
  xlab("")+ylab("Changepoint")+
  ggtitle("Total P (mg/L)")+
  scale_size_continuous("Indicator value\n(z-score)")+
  theme(axis.text.y = element_text(face = "italic", size=6), panel.grid = element_blank())+
  guides(fill = guide_legend(override.aes=list(shape=22, color=NA, size=5)),
         shape = guide_legend(override.aes=list(size=2)))
tp.titan.plot
ggsave(tp.titan.plot, filename="Figures_stratified/TITAN/tp.titan.plot_assemblage_mean.jpg", dpi=300, width=6, height=tp.taxa*.1)

tp.titan.plot2<-
  ggplot(data=droplevels(plot.tp.dat[which(plot.tp.dat$filter>0 & plot.tp.dat$zscore>min.z),]), aes(x=Taxon, y=zenv.cp))+
  # geom_rect(data=plot.tp.sumdat, aes(ymin=cp.05, ymax=cp.95, xmin=-Inf, xmax=Inf, x=NULL, y=NULL), fill="gray",color=NA, alpha=0.35)+
  geom_hline(data=plot.tp.sumdat, aes(yintercept=cp), linetype="dashed", color="white")+
  geom_linerange(aes(ymin=`5%`, ymax=`95%`), color="gray")+
  # geom_point(aes(size=zscore, fill=Assemblage, shape=Group))+
  geom_point(aes(fill=Assemblage, shape=Group), size=1.5)+
  # geom_hline(data=thresholds.bs.df_sum[which(thresholds.bs.df_sum$BSPretty=="Total P" & thresholds.bs.df_sum$BIgoal %in% c("BCG3","BCG4")),], aes(yintercept=MostConservative), linetype="dashed", color="#fed9a6", size=1)+
  # geom_hline(data=thresholds.bs.df_sum[which(thresholds.bs.df_sum$BSPretty=="Total P" & thresholds.bs.df_sum$BIgoal %in% c("Ref30","Ref10","Ref01")),], aes(yintercept=MostConservative), linetype="dotted", color="black", size=0.5)+
  # geom_hline(data=tp_cp_df, aes(yintercept=zenv.cp, color=Assemblage, linetype=Response), size=0.5)+scale_linetype_manual(values=c("dotted","dashed"), name="Mean response")+
  geom_hline(linetype="dashed", yintercept=c(.1276,0.054,0.027,0.1351))+
  scale_shape_manual(values=c(24,25), name="Response")+
  scale_fill_manual(values=c("skyblue","tomato1","darkgreen"))+
  # facet_wrap(~Group, scales="free")+
  scale_y_continuous(trans="log10", breaks=c(0.01,0.05,0.1,0.5,1))+
  coord_flip()+
  theme_bw()+
  xlab("")+ylab("Changepoint")+
  ggtitle("Total P (mg/L)")+
  scale_size_continuous("Indicator value\n(z-score)")+
  theme(axis.text.y = element_text(face = "italic", size=6), panel.grid = element_blank())+
  guides(fill = guide_legend(override.aes=list(shape=22, color=NA, size=5)),
         shape = guide_legend(override.aes=list(size=2)))
tp.titan.plot2
ggsave(tp.titan.plot2, filename="Figures_stratified/TITAN/tp.titan.plot_ref10_probs.jpg", dpi=300, width=6, height=tp.taxa*.1)



chla_cp_df<-data.frame(Assemblage=c("Diatom","Diatom","SBA","SBA","BMI","BMI"),
                       Response=c("Increasers","Decreasers"),
                       zenv.cp=c(47,11,26,36,71,31))

chla.titan.plot<-
  ggplot(data=droplevels(plot.chla.dat[which(plot.chla.dat$filter>0 & plot.chla.dat$zscore>min.z),]), aes(x=Taxon, y=zenv.cp))+
  # geom_rect(data=plot.chla.sumdat, aes(ymin=cp.05, ymax=cp.95, xmin=-Inf, xmax=Inf, x=NULL, y=NULL), fill="gray",color=NA, alpha=0.35)+
  geom_hline(data=plot.chla.sumdat, aes(yintercept=cp), linetype="dashed", color="white")+
  geom_linerange(aes(ymin=`5%`, ymax=`95%`), color="gray")+
  # geom_point(aes(size=zscore, fill=Assemblage, shape=Group))+
  geom_point(aes(fill=Assemblage, shape=Group), size=1.5)+
  # geom_hline(data=thresholds.bs.df_sum[which(thresholds.bs.df_sum$BSPretty=="Chl-a" & thresholds.bs.df_sum$BIgoal %in% c("BCG3","BCG4")),], aes(yintercept=MostConservative), linetype="dashed", color="#fed9a6", size=1)+
  # geom_hline(data=thresholds.bs.df_sum[which(thresholds.bs.df_sum$BSPretty=="Chl-a" & thresholds.bs.df_sum$BIgoal %in% c("Ref30","Ref10","Ref01")),], aes(yintercept=MostConservative), linetype="dotted", color="black", size=0.5)+
  geom_hline(data=chla_cp_df, aes(yintercept=zenv.cp, color=Assemblage, linetype=Response), size=0.5)+scale_linetype_manual(values=c("dotted","dashed"), name="Mean response")+
  scale_shape_manual(values=c(24,25), name="Response")+
  scale_fill_manual(values=c("skyblue","tomato1","darkgreen"))+
  # facet_wrap(~Group, scales="free")+
  scale_y_continuous(trans="log10", breaks=c(5,10,50,100,500,1000))+
  coord_flip()+
  theme_bw()+
  xlab("")+ylab("Changepoint")+
  ggtitle("Chlorophyll-a (mg/m2)")+
  scale_size_continuous("Indicator value\n(z-score)")+
  theme(axis.text.y = element_text(face = "italic", size=6), panel.grid = element_blank())+
  guides(fill = guide_legend(override.aes=list(shape=22, color=NA, size=5)),
         shape = guide_legend(override.aes=list(size=2)))
chla.titan.plot
ggsave(chla.titan.plot, filename="Figures_stratified//TITAN/chla.titan.plot_assemblage_mean.jpg", dpi=300, width=6, height=chla.taxa*.1)


chla.titan.plot2<-
  ggplot(data=droplevels(plot.chla.dat[which(plot.chla.dat$filter>0 & plot.chla.dat$zscore>min.z),]), aes(x=Taxon, y=zenv.cp))+
  # geom_rect(data=plot.chla.sumdat, aes(ymin=cp.05, ymax=cp.95, xmin=-Inf, xmax=Inf, x=NULL, y=NULL), fill="gray",color=NA, alpha=0.35)+
  geom_hline(data=plot.chla.sumdat, aes(yintercept=cp), linetype="dashed", color="white")+
  geom_linerange(aes(ymin=`5%`, ymax=`95%`), color="gray")+
  # geom_point(aes(size=zscore, fill=Assemblage, shape=Group))+
  geom_point(aes(fill=Assemblage, shape=Group), size=1.5)+
  # geom_hline(data=thresholds.bs.df_sum[which(thresholds.bs.df_sum$BSPretty=="Chl-a" & thresholds.bs.df_sum$BIgoal %in% c("BCG3","BCG4")),], aes(yintercept=MostConservative), linetype="dashed", color="#fed9a6", size=1)+
  # geom_hline(data=thresholds.bs.df_sum[which(thresholds.bs.df_sum$BSPretty=="Chl-a" & thresholds.bs.df_sum$BIgoal %in% c("Ref30","Ref10","Ref01")),], aes(yintercept=MostConservative), linetype="dotted", color="black", size=0.5)+
  # geom_hline(data=chla_cp_df, aes(yintercept=zenv.cp, color=Assemblage, linetype=Response), size=0.5)+scale_linetype_manual(values=c("dotted","dashed"), name="Mean response")+
  geom_hline(linetype="dashed", yintercept=c(109.3,44.1,9.009,11.411))+
  scale_shape_manual(values=c(24,25), name="Response")+
  scale_fill_manual(values=c("skyblue","tomato1","darkgreen"))+
  # facet_wrap(~Group, scales="free")+
  scale_y_continuous(trans="log10", breaks=c(5,10,50,100,500,1000))+
  coord_flip()+
  theme_bw()+
  xlab("")+ylab("Changepoint")+
  ggtitle("Chlorophyll-a (mg/m2)")+
  scale_size_continuous("Indicator value\n(z-score)")+
  theme(axis.text.y = element_text(face = "italic", size=6), panel.grid = element_blank())+
  guides(fill = guide_legend(override.aes=list(shape=22, color=NA, size=5)),
         shape = guide_legend(override.aes=list(size=2)))
chla.titan.plot2
ggsave(chla.titan.plot2, filename="Figures_stratified//TITAN/chla.titan.plot_ref10_probs.jpg", dpi=300, width=6, height=chla.taxa*.1)


afdm_cp_df<-data.frame(Assemblage=c("Diatom","Diatom","SBA","SBA","BMI","BMI"),
                       Response=c("Increasers","Decreasers"),
                       zenv.cp=c(1.8,1.1,1.9,1.5,3.1,2))

afdm.titan.plot<-
  ggplot(data=droplevels(plot.afdm.dat[which(plot.afdm.dat$filter>0 & plot.afdm.dat$zscore>min.z),]), aes(x=Taxon, y=zenv.cp))+
  # geom_rect(data=plot.afdm.sumdat, aes(ymin=cp.05, ymax=cp.95, xmin=-Inf, xmax=Inf, x=NULL, y=NULL), fill="gray",color=NA, alpha=0.35)+
  geom_hline(data=plot.afdm.sumdat, aes(yintercept=cp), linetype="dashed", color="white")+
  geom_linerange(aes(ymin=`5%`, ymax=`95%`), color="gray")+
  # geom_point(aes(size=zscore, fill=Assemblage, shape=Group))+
  geom_point(aes(fill=Assemblage, shape=Group), size=1.5)+
  # geom_hline(data=thresholds.bs.df_sum[which(thresholds.bs.df_sum$BSPretty=="AFDM" & thresholds.bs.df_sum$BIgoal %in% c("BCG3","BCG4")),], aes(yintercept=MostConservative), linetype="dashed", color="#fed9a6", size=1)+
  # geom_hline(data=thresholds.bs.df_sum[which(thresholds.bs.df_sum$BSPretty=="AFDM" & thresholds.bs.df_sum$BIgoal %in% c("Ref30","Ref10","Ref01")),], aes(yintercept=MostConservative), linetype="dotted", color="black", size=0.5)+
  geom_hline(data=afdm_cp_df, aes(yintercept=zenv.cp, color=Assemblage, linetype=Response), size=0.5)+scale_linetype_manual(values=c("dotted","dashed"), name="Mean response")+  
  scale_shape_manual(values=c(24,25), name="Response")+
  scale_fill_manual(values=c("skyblue","tomato1","darkgreen"))+
  # facet_wrap(~Group, scales="free")+
  scale_y_continuous(trans="log10", breaks=c(0.5,1,5,10))+
  coord_flip()+
  theme_bw()+
  xlab("")+ylab("Changepoint")+
  ggtitle("Ash-Free Dry Mass (mg/cm2)")+
  scale_size_continuous("Indicator value\n(z-score)")+
  theme(axis.text.y = element_text(face = "italic", size=6), panel.grid = element_blank())+
  guides(fill = guide_legend(override.aes=list(shape=22, color=NA, size=5)),
         shape = guide_legend(override.aes=list(size=2)))
afdm.titan.plot
ggsave(afdm.titan.plot, filename="Figures_stratified/TITAN/afdm.titan.plot_assemblage_mean.jpg", dpi=300, width=6, height=afdm.taxa*.1)


afdm.titan.plot2<-
  ggplot(data=droplevels(plot.afdm.dat[which(plot.afdm.dat$filter>0 & plot.afdm.dat$zscore>min.z),]), aes(x=Taxon, y=zenv.cp))+
  # geom_rect(data=plot.afdm.sumdat, aes(ymin=cp.05, ymax=cp.95, xmin=-Inf, xmax=Inf, x=NULL, y=NULL), fill="gray",color=NA, alpha=0.35)+
  geom_hline(data=plot.afdm.sumdat, aes(yintercept=cp), linetype="dashed", color="white")+
  geom_linerange(aes(ymin=`5%`, ymax=`95%`), color="gray")+
  # geom_point(aes(size=zscore, fill=Assemblage, shape=Group))+
  geom_point(aes(fill=Assemblage, shape=Group), size=1.5)+
  # geom_hline(data=thresholds.bs.df_sum[which(thresholds.bs.df_sum$BSPretty=="AFDM" & thresholds.bs.df_sum$BIgoal %in% c("BCG3","BCG4")),], aes(yintercept=MostConservative), linetype="dashed", color="#fed9a6", size=1)+
  # geom_hline(data=thresholds.bs.df_sum[which(thresholds.bs.df_sum$BSPretty=="AFDM" & thresholds.bs.df_sum$BIgoal %in% c("Ref30","Ref10","Ref01")),], aes(yintercept=MostConservative), linetype="dotted", color="black", size=0.5)+
  # geom_hline(data=afdm_cp_df, aes(yintercept=zenv.cp, color=Assemblage, linetype=Response), size=0.5)+scale_linetype_manual(values=c("dotted","dashed"), name="Mean response")+
  geom_hline(linetype="dashed", yintercept=c(6.08,2.477,1.276276,0.600601))+
  scale_shape_manual(values=c(24,25), name="Response")+
  scale_fill_manual(values=c("skyblue","tomato1","darkgreen"))+
  # facet_wrap(~Group, scales="free")+
  scale_y_continuous(trans="log10", breaks=c(0.5,1,5,10))+
  coord_flip()+
  theme_bw()+
  xlab("")+ylab("Changepoint")+
  ggtitle("Ash-Free Dry Mass (mg/cm2)")+
  scale_size_continuous("Indicator value\n(z-score)")+
  theme(axis.text.y = element_text(face = "italic", size=6), panel.grid = element_blank())+
  guides(fill = guide_legend(override.aes=list(shape=22, color=NA, size=5)),
         shape = guide_legend(override.aes=list(size=2)))
afdm.titan.plot2
ggsave(afdm.titan.plot2, filename="Figures_stratified/TITAN/afdm.titan.plot_ref10_probs.jpg", dpi=300, width=6, height=afdm.taxa*.1)


map_cp_df<-data.frame(Assemblage=c("Diatom","Diatom","SBA","SBA","BMI","BMI"),
                      Response=c("Increasers","Decreasers"),
                      zenv.cp=c(17,18,16,23,68,28))

map.titan.plot<-
  ggplot(data=droplevels(plot.map.dat[which(plot.map.dat$filter>0 & plot.map.dat$zscore>min.z),]), aes(x=Taxon, y=zenv.cp))+
  # geom_rect(data=plot.map.sumdat, aes(ymin=cp.05, ymax=cp.95, xmin=-Inf, xmax=Inf, x=NULL, y=NULL), fill="gray",color=NA, alpha=0.35)+
  geom_hline(data=plot.map.sumdat, aes(yintercept=cp), linetype="dashed", color="white")+
  geom_linerange(aes(ymin=`5%`, ymax=`95%`), color="gray")+
  # geom_point(aes(size=zscore, fill=Assemblage, shape=Group))+
  geom_point(aes(fill=Assemblage, shape=Group), size=1.5)+
  # geom_hline(data=thresholds.bs.df_sum[which(thresholds.bs.df_sum$BSPretty=="% cover" & thresholds.bs.df_sum$BIgoal %in% c("BCG3","BCG4")),], aes(yintercept=MostConservative), linetype="dashed", color="#fed9a6", size=1)+
  # geom_hline(data=thresholds.bs.df_sum[which(thresholds.bs.df_sum$BSPretty=="% cover" & thresholds.bs.df_sum$BIgoal %in% c("Ref30","Ref10","Ref01")),], aes(yintercept=MostConservative), linetype="dotted", color="black", size=0.5)+
  geom_hline(data=map_cp_df, aes(yintercept=zenv.cp, color=Assemblage, linetype=Response), size=0.5)+scale_linetype_manual(values=c("dotted","dashed"), name="Mean response")+  
  scale_shape_manual(values=c(24,25), name="Response")+
  scale_fill_manual(values=c("skyblue","tomato1","darkgreen"))+
  # facet_wrap(~Group, scales="free")+
  scale_y_continuous()+
  coord_flip()+
  theme_bw()+
  xlab("")+ylab("Changepoint")+
  ggtitle("% macroalgal cover")+
  scale_size_continuous("Indicator value\n(z-score)")+
  theme(axis.text.y = element_text(face = "italic", size=6), panel.grid = element_blank())+
  guides(fill = guide_legend(override.aes=list(shape=22, color=NA, size=5)),
         shape = guide_legend(override.aes=list(size=2)))
map.titan.plot
ggsave(map.titan.plot, filename="Figures_stratified/TITAN/map.titan.plot_assemblage_mean.jpg", dpi=300, width=6, height=map.taxa*.11)

map.titan.plot2<-
  ggplot(data=droplevels(plot.map.dat[which(plot.map.dat$filter>0 & plot.map.dat$zscore>min.z),]), aes(x=Taxon, y=zenv.cp))+
  # geom_rect(data=plot.map.sumdat, aes(ymin=cp.05, ymax=cp.95, xmin=-Inf, xmax=Inf, x=NULL, y=NULL), fill="gray",color=NA, alpha=0.35)+
  geom_hline(data=plot.map.sumdat, aes(yintercept=cp), linetype="dashed", color="white")+
  geom_linerange(aes(ymin=`5%`, ymax=`95%`), color="gray")+
  # geom_point(aes(size=zscore, fill=Assemblage, shape=Group))+
  geom_point(aes(fill=Assemblage, shape=Group), size=1.5)+
  # geom_hline(data=thresholds.bs.df_sum[which(thresholds.bs.df_sum$BSPretty=="% cover" & thresholds.bs.df_sum$BIgoal %in% c("BCG3","BCG4")),], aes(yintercept=MostConservative), linetype="dashed", color="#fed9a6", size=1)+
  # geom_hline(data=thresholds.bs.df_sum[which(thresholds.bs.df_sum$BSPretty=="% cover" & thresholds.bs.df_sum$BIgoal %in% c("Ref30","Ref10","Ref01")),], aes(yintercept=MostConservative), linetype="dotted", color="black", size=0.5)+
  # geom_hline(data=map_cp_df, aes(yintercept=zenv.cp, color=Assemblage, linetype=Response), size=0.5)+scale_linetype_manual(values=c("dotted","dashed"), name="Mean response")+
  geom_hline(linetype="dashed", yintercept=c(67.76777,26.72673,13.61361,6.806807))+
  scale_shape_manual(values=c(24,25), name="Response")+
  scale_fill_manual(values=c("skyblue","tomato1","darkgreen"))+
  # facet_wrap(~Group, scales="free")+
  scale_y_continuous()+
  coord_flip()+
  theme_bw()+
  xlab("")+ylab("Changepoint")+
  ggtitle("% macroalgal cover")+
  scale_size_continuous("Indicator value\n(z-score)")+
  theme(axis.text.y = element_text(face = "italic", size=6), panel.grid = element_blank())+
  guides(fill = guide_legend(override.aes=list(shape=22, color=NA, size=5)),
         shape = guide_legend(override.aes=list(size=2)))
map.titan.plot2
ggsave(map.titan.plot2, filename="Figures_stratified/TITAN/map.titan.plot_ref10_probs.jpg", dpi=300, width=6, height=map.taxa*.11)



####Other responses
plot.tn.sumdat$BSPretty<-"Total N"
plot.tp.sumdat$BSPretty<-"Total P"
plot.chla.sumdat$BSPretty<-"Chl-a"
plot.afdm.sumdat$BSPretty<-"AFDM"
plot.map.sumdat$BSPretty<-"% cover"

assemblage.responses<-rbind.fill(plot.tn.sumdat,
                                 plot.tp.sumdat,
                                 plot.chla.sumdat,
                                 plot.afdm.sumdat,
                                 plot.map.sumdat)

write.csv(assemblage.responses, "R_Outputs_082318/TITAN/assemblage.responses.csv", row.names=F)

ref.dists.x<-read.csv("R_Outputs_082318/state_summary.ref.table2.csv") 

ref.dists<-read.csv("R_Outputs_082318/state_summary.ref.table2.csv") %>%
  filter(BSPretty %in% bs.varz.df$BSPretty & Dist==0.9) %>%
  mutate(#Group="Not applicable",
    Type="Ref-90th",
    l95=NA, u95=NA) %>%
  dplyr::rename(Est=Ref) %>%
  dplyr::select(BSPretty,Type,Est, l95, u95)
titan.ass.dists<-assemblage.responses %>% 
  mutate(Type=paste(Assemblage, Group)) %>%
  dplyr::rename(Est=cp,
                l95=cp.05, u95=cp.95) %>%
  dplyr::select(BSPretty,Type,Est, l95, u95)

other.dists<-rbind.fill(ref.dists, titan.ass.dists)  
other.dists$Type<-factor(other.dists$Type, levels=c("BMI Increasers", "BMI Decreasers", "Diatom Increasers", "SBA Increasers", "Diatom Decreasers", "SBA Decreasers", "Ref-90th"))

myvar.i<-"Total N"
other.tn.plot<-  ggplot(data=other.dists[which(other.dists$BSPretty==myvar.i),], aes(x=Type, y=Est))+
  geom_linerange(aes(ymin=l95, ymax=u95), color="gray")+
  geom_point()+
  # geom_hline(data=thresholds.bs.df_sum[which(thresholds.bs.df_sum$BSPretty==myvar.i & thresholds.bs.df_sum$BIgoal %in% c("BCG3","BCG4")),], aes(yintercept=MostConservative), linetype="dashed", color="#fed9a6", size=1)+
  geom_hline(data=thresholds.bs.df_sum[which(thresholds.bs.df_sum$BSPretty==myvar.i & thresholds.bs.df_sum$BIgoal %in% c("Ref30","Ref10","Ref01")),], aes(yintercept=MostConservative), linetype="dotted", color="black", size=0.5)+
  scale_y_continuous()+
  coord_flip()+
  theme_bw(base_size=10)+
  xlab("")+ylab("Total N (mg/L)")+
  theme(panel.grid = element_blank())
other.tn.plot
ggsave(other.tn.plot,filename="R_Outputs_082318/TITAN/other.tn.plot.jpg", height=3.5, width=3.6)

myvar.i<-"Total P"
other.tp.plot<-  ggplot(data=other.dists[which(other.dists$BSPretty==myvar.i),], aes(x=Type, y=Est))+
  geom_linerange(aes(ymin=l95, ymax=u95), color="gray")+
  geom_point()+
  # geom_hline(data=thresholds.bs.df_sum[which(thresholds.bs.df_sum$BSPretty==myvar.i & thresholds.bs.df_sum$BIgoal %in% c("BCG3","BCG4")),], aes(yintercept=MostConservative), linetype="dashed", color="#fed9a6", size=1)+
  geom_hline(data=thresholds.bs.df_sum[which(thresholds.bs.df_sum$BSPretty==myvar.i & thresholds.bs.df_sum$BIgoal %in% c("Ref30","Ref10","Ref01")),], aes(yintercept=MostConservative), linetype="dotted", color="black", size=0.5)+
  scale_y_continuous()+
  coord_flip()+
  theme_bw(base_size=10)+
  xlab("")+ylab("Total P (mg/L)")+
  theme(panel.grid = element_blank())
other.tp.plot
ggsave(other.tp.plot,filename="R_Outputs_082318/TITAN/other.tp.plot.jpg", height=3.5, width=3.6)

myvar.i<-"Chl-a"
other.chla.plot<-  ggplot(data=other.dists[which(other.dists$BSPretty==myvar.i),], aes(x=Type, y=Est))+
  geom_linerange(aes(ymin=l95, ymax=u95), color="gray")+
  geom_point()+
  # geom_hline(data=thresholds.bs.df_sum[which(thresholds.bs.df_sum$BSPretty==myvar.i & thresholds.bs.df_sum$BIgoal %in% c("BCG3","BCG4")),], aes(yintercept=MostConservative), linetype="dashed", color="#fed9a6", size=1)+
  geom_hline(data=thresholds.bs.df_sum[which(thresholds.bs.df_sum$BSPretty==myvar.i & thresholds.bs.df_sum$BIgoal %in% c("Ref30","Ref10","Ref01")),], aes(yintercept=MostConservative), linetype="dotted", color="black", size=0.5)+
  scale_y_continuous(breaks=c(0,25,50,75,100))+
  coord_flip()+
  theme_bw(base_size=10)+
  xlab("")+ylab("Chlorophyll-a (mg/m2)")+
  theme(panel.grid = element_blank())
other.chla.plot
ggsave(other.chla.plot,filename="R_Outputs_082318/TITAN/other.chla.plot.jpg", height=3.5, width=3.6)

myvar.i<-"AFDM"
other.afdm.plot<-  ggplot(data=other.dists[which(other.dists$BSPretty==myvar.i),], aes(x=Type, y=Est))+
  geom_linerange(aes(ymin=l95, ymax=u95), color="gray")+
  geom_point()+
  # geom_hline(data=thresholds.bs.df_sum[which(thresholds.bs.df_sum$BSPretty==myvar.i & thresholds.bs.df_sum$BIgoal %in% c("BCG3","BCG4")),], aes(yintercept=MostConservative), linetype="dashed", color="#fed9a6", size=1)+
  geom_hline(data=thresholds.bs.df_sum[which(thresholds.bs.df_sum$BSPretty==myvar.i & thresholds.bs.df_sum$BIgoal %in% c("Ref30","Ref10","Ref01")),], aes(yintercept=MostConservative), linetype="dotted", color="black", size=0.5)+
  scale_y_continuous(breaks=c(0,1,2,3,4))+
  coord_flip()+
  theme_bw(base_size=8)+
  xlab("")+ylab("Ash-Free Dry Mass (mg/cm2)")+
  theme(panel.grid = element_blank())
other.afdm.plot
ggsave(other.afdm.plot,filename="R_Outputs_082318/TITAN/other.afdm.plot.jpg", height=3.5, width=3.6)

myvar.i<-"% cover"
other.map.plot<-  ggplot(data=other.dists[which(other.dists$BSPretty==myvar.i),], aes(x=Type, y=Est))+
  geom_linerange(aes(ymin=l95, ymax=u95), color="gray")+
  geom_point()+
  # geom_hline(data=thresholds.bs.df_sum[which(thresholds.bs.df_sum$BSPretty==myvar.i & thresholds.bs.df_sum$BIgoal %in% c("BCG3","BCG4")),], aes(yintercept=MostConservative), linetype="dashed", color="#fed9a6", size=1)+
  geom_hline(data=thresholds.bs.df_sum[which(thresholds.bs.df_sum$BSPretty==myvar.i & thresholds.bs.df_sum$BIgoal %in% c("Ref30","Ref10","Ref01")),], aes(yintercept=MostConservative), linetype="dotted", color="black", size=0.5)+
  scale_y_continuous()+
  coord_flip()+
  theme_bw(base_size=8)+
  xlab("")+ylab("% macroalgal cover")+
  theme(panel.grid = element_blank())
other.map.plot
ggsave(other.map.plot,filename="R_Outputs_082318/TITAN/other.map.plot.jpg", height=3.5, width=3.6)

write.csv(other.dists,"R_Outputs_082318/TITAN/other.dists.csv", row.names=F)

# 
# ggplot(data=droplevels(plot.tn.dat[which(plot.tn.dat$filter>0),]), aes(x=Taxon, y=zenv.cp, color=Group))+
#   geom_hline(data=plot.tn.sumdat, aes(yintercept=cp, color=Group), linetype="dashed")+
#   geom_linerange(aes(ymin=`5%`, ymax=`95%`))+
#   geom_point(aes(size=zscore))+
#   facet_wrap(~Assemblage, scales="free")+
#   scale_y_continuous(trans="log10")+
#   coord_flip()+
#   theme_bw()+
#   xlab("")+ylab("Changepoint\nTotal N (mg/L)")+
#   theme(axis.text = element_text(face = "italic"))
# 
# #
# 

#Ntaxa summary
ntax.summary<-expand.grid(Assemblage=c("BMI","Diatom","SBA"),
                          BSPretty=bs.varz.df$BSPretty,
                          BIgoal=c("Ref30","Ref10","Ref01"),
                          # Protected=c("Yes","No"),
                          stringsAsFactors = F)
ntax.summary$ntax_eval<-sapply(1:nrow(ntax.summary), function(i){
  ass.i<-ntax.summary$Assemblage[i]
  if(ass.i=="BMI")
    dim(bug_tax.cal)[2]
  else
    if(ass.i=="Diatom")
      dim(diatom_tax.cal)[2]
  else
    if(ass.i=="SBA")
      dim(sba_tax.cal)[2]
  else
    NA
})

plot.tn.dat %>%
  group_by(Assemblage) %>%
  filter(Group!="Neither" & zenv.cp<tn.ref10) %>%
  ddply(.(Assemblage),summarise, ntaxa=length(Taxon))

i<-3

thresholds.bs.df_sum<-read_csv("Outputs_stratified/tab.threshold.rr.summary.csv") %>%
  filter(Stratum=="California" & 
           Prob=="p80" & 
           RR.l95.cal>1 &
           RR.l95.val>1) %>%
  group_by(BSPretty, BIgoal) %>%
  summarize(MostConservative=min(Est)) %>%
  ungroup()

ntax.summary$ntax_LowerCP<-sapply(1:nrow(ntax.summary), function(i){
  ass.i<-as.character(ntax.summary$Assemblage[i])
  goal.i<-as.character(ntax.summary$BIgoal[i])
  var.i<-as.character(ntax.summary$BSPretty[i])
  targ.i<-thresholds.bs.df_sum$MostConservative[which(thresholds.bs.df_sum$BSPretty==var.i & thresholds.bs.df_sum$BIgoal==goal.i)]
  if(var.i=="Total N"){
    x.df<-plot.tn.dat[which(plot.tn.dat$Assemblage==ass.i & plot.tn.dat$Group!="Neither" & plot.tn.dat$zenv.cp<targ.i),]
    return(length(x.df$Taxon)  )}
  if(var.i=="Total P"){
    x.df<-plot.tp.dat[which(plot.tp.dat$Assemblage==ass.i & plot.tp.dat$Group!="Neither" & plot.tp.dat$zenv.cp<targ.i),]
    return(length(x.df$Taxon))    }
  if(var.i=="Chl-a"){
    x.df<-plot.chla.dat[which(plot.chla.dat$Assemblage==ass.i & plot.chla.dat$Group!="Neither" & plot.chla.dat$zenv.cp<targ.i),]
    return(length(x.df$Taxon))  }
  if(var.i=="AFDM"){
    x.df<-plot.afdm.dat[which(plot.afdm.dat$Assemblage==ass.i & plot.afdm.dat$Group!="Neither" & plot.afdm.dat$zenv.cp<targ.i),]
    return(length(x.df$Taxon))    }
  if(var.i=="% cover"){
    x.df<-plot.map.dat[which(plot.map.dat$Assemblage==ass.i & plot.map.dat$Group!="Neither" & plot.map.dat$zenv.cp<targ.i),]
    return(length(x.df$Taxon))    }
})

ntax.summary$ntax_HigherCP<-sapply(1:nrow(ntax.summary), function(i){
  ass.i<-as.character(ntax.summary$Assemblage[i])
  goal.i<-as.character(ntax.summary$BIgoal[i])
  var.i<-as.character(ntax.summary$BSPretty[i])
  targ.i<-thresholds.bs.df_sum$MostConservative[which(thresholds.bs.df_sum$BSPretty==var.i & thresholds.bs.df_sum$BIgoal==goal.i)]
  if(var.i=="Total N"){
    x.df<-plot.tn.dat[which(plot.tn.dat$Assemblage==ass.i & plot.tn.dat$Group!="Neither" & plot.tn.dat$zenv.cp>=targ.i),]
    return(length(x.df$Taxon)  )}
  if(var.i=="Total P"){
    x.df<-plot.tp.dat[which(plot.tp.dat$Assemblage==ass.i & plot.tp.dat$Group!="Neither" & plot.tp.dat$zenv.cp>=targ.i),]
    return(length(x.df$Taxon))    }
  if(var.i=="Chl-a"){
    x.df<-plot.chla.dat[which(plot.chla.dat$Assemblage==ass.i & plot.chla.dat$Group!="Neither" & plot.chla.dat$zenv.cp>=targ.i),]
    return(length(x.df$Taxon))  }
  if(var.i=="AFDM"){
    x.df<-plot.afdm.dat[which(plot.afdm.dat$Assemblage==ass.i & plot.afdm.dat$Group!="Neither" & plot.afdm.dat$zenv.cp>=targ.i),]
    return(length(x.df$Taxon))    }
  if(var.i=="% cover"){
    x.df<-plot.map.dat[which(plot.map.dat$Assemblage==ass.i & plot.map.dat$Group!="Neither" & plot.map.dat$zenv.cp>=targ.i),]
    return(length(x.df$Taxon))    }
})

head(ntax.summary)
ntax.summary$Indifferent<-ntax.summary$ntax_eval-(ntax.summary$ntax_LowerCP+ntax.summary$ntax_HigherCP)

ntax.summary.m<-#melt(ntax.summary, id.var=c("Assemblage","BSPretty","BIgoal","ntax_eval"))
  ntax.summary %>%
  pivot_longer(cols=c(ntax_LowerCP, ntax_HigherCP, Indifferent))
ntax.summary.m$Pct<-100*ntax.summary.m$value/ntax.summary.m$ntax_eval
ntax.summary.m$name<-factor(ntax.summary.m$name, levels=c("ntax_LowerCP","ntax_HigherCP","Indifferent"),
                            # labels=c("No","Yes (sensitive)","Yes (indifferent)"))
                            labels=c("No","Yes","Yes (no cp detected)"))
# head(ntax.summary.m)
# protected.taxa.plot<-
#   ggplot(data=ntax.summary.m, aes(x=BIgoal, y=Pct))+
#   geom_bar(aes(fill=variable), position=position_stack(), stat="identity", color=NA)+
#   scale_fill_manual(values=c("tomato1","#1f78b4","#a6cee3"), name="Is changepoint above threshold?")+
#   facet_grid(Assemblage~BSPretty,scales="free")+
#   scale_y_continuous(name="% of taxa",breaks=c(25,50,75,100), expand=c(0,0))+
#   xlab("Biointegrity basis for eutrophication threshold")+
#   theme_bw(base_size = 10)+
#   theme(legend.position = "bottom", panel.grid = element_blank(), panel.border = element_blank(),
#         strip.background = element_blank())
# ggsave(protected.taxa.plot, filename="R_Outputs_082318/TITAN/protected.taxa.plot.jpg", dpi=300, height=5, width=6.5)

protected.taxa.plot<-
  ggplot(data=ntax.summary.m, aes(x=BIgoal, y=value))+
  geom_bar(aes(fill=name), position=position_stack(), stat="identity", color=NA)+
  scale_fill_manual(values=c("tomato1","#1f78b4","#a6cee3"), name="Is changepoint above threshold?")+
  facet_grid(Assemblage~BSPretty,scales="free")+
  scale_y_continuous(name="# of taxa", expand=c(0,0))+
  xlab("Biointegrity basis for eutrophication threshold")+
  theme_bw(base_size = 10)+
  theme(legend.position = "bottom", panel.grid = element_blank(), panel.border = element_blank(),
        strip.background = element_blank())
ggsave(protected.taxa.plot, filename="Figures_stratified/protected.taxa.plot_p80.jpg", dpi=300, height=5, width=6.5)

head(other.dists)
# other.dists2<-other.dists
# other.dists2$BSPretty<-factor(other.dists2$BSPretty, levels=rev(bs.varz.df$BSPretty))
other.all.plot<-  ggplot(data=other.dists, aes(x=Type, y=Est))+
  geom_linerange(aes(ymin=l95, ymax=u95), color="gray")+
  geom_point()+
  # geom_hline(data=thresholds.bs.df_sum[which(thresholds.bs.df_sum$BSPretty==myvar.i & thresholds.bs.df_sum$BIgoal %in% c("BCG3","BCG4")),], aes(yintercept=MostConservative), linetype="dashed", color="#fed9a6", size=1)+
  geom_hline(data=thresholds.bs.df_sum[which(thresholds.bs.df_sum$BIgoal %in% c("Ref30","Ref10","Ref01")),], aes(yintercept=MostConservative), linetype="dotted", color="black", size=0.5)+
  scale_y_continuous()+
  coord_flip()+
  theme_bw(base_size=8)+
  facet_wrap(~BSPretty, scales="free_x", strip.position = "bottom")+
  xlab("")+ylab("")+
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.placement = "outside")
other.all.plot
ggsave(other.all.plot,filename="R_Outputs_082318/TITAN/other.all.plot.jpg", height=5, width=6.5)

alldata<-rbind.fill(plot.tn.dat, plot.tp.dat, plot.chla.dat, plot.afdm.dat, plot.map.dat)
alldata.2<-alldata[which(alldata$Group!="Neither"),]
ddply(alldata, .(Assemblage), summarize, ntax=length(unique(Taxon[which(Group!="Neither")]))/length(unique(Taxon)))

ddply(ntax.summary, .(BIgoal,BSPretty), summarise, nVul=sum(ntax_LowerCP))


#################
#New figures 05312022

dummy_line_plot_dat<-tibble(
  Series = c("BMI CP","Diatom CP","SBA CP","Ref10 response model","CP 95% confidence interval"),
  Series2=c("Changepoint","Changepoint","Changepoint","Response model threshold","CI"),
  v1=c(1,2,3,4,5),
  v2=c(1,2,3,4,5)+10)
dummy_line_plot_dat$Series<-factor(dummy_line_plot_dat$Series, levels=dummy_line_plot_dat$Series)

ggplot(data=dummy_line_plot_dat)+
  geom_line(aes(x=v1, y=v2,
                color=Series, linetype=Series))+
  scale_color_manual(values=c("skyblue","tomato1","darkgreen","red","gray"), name="")+
  scale_linetype_manual(values=c("dotted","dotted","dotted","dashed", "solid"), name="")+
  theme_bw()+
  guides(linetype = guide_legend(override.aes=list(size=1, angle=0)))

levels(plot.tn.dat$Group)<-c("Neither", "Increaser taxa","Decreaser taxa")
plot.tn.sumdat<-plot.tn.sumdat %>%  mutate(Group=case_when(Group=="Increasers"~"Increaser taxa",T~"Decreaser taxa"))

tn.titan.plot3<-
  ggplot(data=plot.tn.dat %>%
           filter(filter>0 & zscore>min.z) %>%
           droplevels(),
         aes(x=Taxon, y=zenv.cp))+
  geom_linerange(aes(ymin=`5%`, ymax=`95%`, color="gray"),
                 key_glyph="path")+ 
  scale_color_manual(values='gray', name="", labels="95% confidence interval")+
  geom_point(aes(fill=Assemblage), size=1.5, shape=21)+
  scale_fill_manual(values=c("skyblue","tomato1","darkgreen"), labels=c("Benthic macroinvertebrates","Diatoms","Soft-bodied algae"))+
  facet_wrap(~Group, scales="free")+
  theme_bw(base_size=10)+
  xlab("")+ylab("Changepoint")+
  ggtitle("Total N (mg/L)")+
  scale_size_continuous("Indicator value\n(z-score)")+
  theme(axis.text.y = element_text(face = "italic", size=6), panel.grid = element_blank())+
  guides(fill = guide_legend(override.aes=list(shape=21, size=3) , order=1),
         color = guide_legend(order=2))+
  coord_flip()
tn.titan.plot3
ggsave(tn.titan.plot3, filename="Figures_stratified/TITAN/tn.titan.plot_faceted.jpg", dpi=300, width=6, height=tn.taxa*.1)


tn.titan.plot3_reflines<-
  ggplot(data=plot.tn.dat %>%
           filter(filter>0 & zscore>min.z) %>%
           droplevels(),
         aes(x=Taxon, y=zenv.cp))+
  geom_linerange(aes(ymin=`5%`, ymax=`95%`), color="gray")+
  geom_hline(data=plot.tn.sumdat, aes(yintercept=cp, color=Assemblage), linetype="dashed", key_glyph="vline")+
  scale_color_manual(values=c("skyblue","tomato1","darkgreen"), name="Mean changepoint", labels=c("Benthic macroinvertebrates","Diatoms","Soft-bodied algae"))+
  
  geom_point(aes(fill=Assemblage), size=1.5, shape=21)+
  scale_fill_manual(values=c("skyblue","tomato1","darkgreen"), labels=c("Benthic macroinvertebrates","Diatoms","Soft-bodied algae"))+
  facet_wrap(~Group, scales="free")+
  theme_bw(base_size=10)+
  xlab("")+ylab("Changepoint")+
  ggtitle("Total N (mg/L)")+
  scale_size_continuous("Indicator value\n(z-score)")+
  
  
  theme(axis.text.y = element_text(face = "italic", size=6), panel.grid = element_blank())+
  guides(fill = guide_legend(override.aes=list(shape=21, size=3) , order=1),
         color = guide_legend(order=2))+
  coord_flip()
tn.titan.plot3_reflines
ggsave(tn.titan.plot3_reflines, filename="Figures_stratified/TITAN/tn.titan.plot3_reflines.jpg", dpi=300, width=6, height=tn.taxa*.1)


levels(plot.tp.dat$Group)<-c("Neither", "Increaser taxa","Decreaser taxa")
plot.tp.sumdat<-plot.tp.sumdat %>%  mutate(Group=case_when(Group=="Increasers"~"Increaser taxa",T~"Decreaser taxa"))

tp.titan.plot3<-
  ggplot(data=plot.tp.dat %>%
           filter(filter>0 & zscore>min.z) %>%
           droplevels(),
         aes(x=Taxon, y=zenv.cp))+
  geom_linerange(aes(ymin=`5%`, ymax=`95%`, color="gray"),
                 key_glyph="path")+ 
  scale_color_manual(values='gray', name="", labels="95% confidence interval")+
  geom_point(aes(fill=Assemblage), size=1.5, shape=21)+
  scale_fill_manual(values=c("skyblue","tomato1","darkgreen"), labels=c("Benthic macroinvertebrates","Diatoms","Soft-bodied algae"))+
  facet_wrap(~Group, scales="free")+
  theme_bw(base_size=10)+
  xlab("")+ylab("Changepoint")+
  ggtitle("Total P (mg/L)")+
  scale_size_continuous("Indicator value\n(z-score)")+
  theme(axis.text.y = element_text(face = "italic", size=6), panel.grid = element_blank())+
  guides(fill = guide_legend(override.aes=list(shape=21, size=3) , order=1),
         color = guide_legend(order=2))+
  coord_flip()
tp.titan.plot3
ggsave(tp.titan.plot3, filename="Figures_stratified/TITAN/tp.titan.plot_faceted.jpg", dpi=300, width=6, height=tn.taxa*.1)



tp.titan.plot3_reflines<-
  ggplot(data=plot.tp.dat %>%
           filter(filter>0 & zscore>min.z) %>%
           droplevels(),
         aes(x=Taxon, y=zenv.cp))+
  geom_linerange(aes(ymin=`5%`, ymax=`95%`), color="gray")+
  geom_hline(data=plot.tp.sumdat, aes(yintercept=cp, color=Assemblage), linetype="dashed", key_glyph="vline")+
  scale_color_manual(values=c("skyblue","tomato1","darkgreen"), name="Mean changepoint", labels=c("Benthic macroinvertebrates","Diatoms","Soft-bodied algae"))+
  
  geom_point(aes(fill=Assemblage), size=1.5, shape=21)+
  scale_fill_manual(values=c("skyblue","tomato1","darkgreen"), labels=c("Benthic macroinvertebrates","Diatoms","Soft-bodied algae"))+
  facet_wrap(~Group, scales="free")+
  theme_bw(base_size=10)+
  xlab("")+ylab("Changepoint")+
  ggtitle("Total P (mg/L)")+
  scale_size_continuous("Indicator value\n(z-score)")+
  theme(axis.text.y = element_text(face = "italic", size=6), panel.grid = element_blank())+
  guides(fill = guide_legend(override.aes=list(shape=21, size=3) , order=1),
         color = guide_legend(order=2))+
  coord_flip()
tp.titan.plot3_reflines
ggsave(tp.titan.plot3_reflines, filename="Figures_stratified/TITAN/tp.titan.plot3_reflines.jpg", dpi=300, width=6, height=tn.taxa*.1)

levels(plot.chla.dat$Group)<-c("Neither", "Increaser taxa","Decreaser taxa")
plot.chla.sumdat<-plot.chla.sumdat %>%  mutate(Group=case_when(Group=="Increasers"~"Increaser taxa",T~"Decreaser taxa"))

chla.titan.plot3<-
  ggplot(data=plot.chla.dat %>%
           filter(filter>0 & zscore>min.z) %>%
           droplevels(),
         aes(x=Taxon, y=zenv.cp))+
  geom_linerange(aes(ymin=`5%`, ymax=`95%`, color="gray"),
                 key_glyph="path")+ 
  scale_color_manual(values='gray', name="", labels="95% confidence interval")+
  geom_point(aes(fill=Assemblage), size=1.5, shape=21)+
  scale_fill_manual(values=c("skyblue","tomato1","darkgreen"), labels=c("Benthic macroinvertebrates","Diatoms","Soft-bodied algae"))+
  facet_wrap(~Group, scales="free")+
  theme_bw(base_size=10)+
  xlab("")+ylab("Changepoint")+
  ggtitle("Chlorophyll a (mg/m2)")+
  scale_size_continuous("Indicator value\n(z-score)")+
  theme(axis.text.y = element_text(face = "italic", size=6), panel.grid = element_blank())+
  guides(fill = guide_legend(override.aes=list(shape=21, size=3) , order=1),
         color = guide_legend(order=2))+
  coord_flip()
chla.titan.plot3
ggsave(chla.titan.plot3, filename="Figures_stratified/TITAN/chla.titan.plot_faceted.jpg", dpi=300, width=6, height=tn.taxa*.1)


chla.titan.plot3_reflines<-
  ggplot(data=plot.chla.dat %>%
           filter(filter>0 & zscore>min.z) %>%
           droplevels(),
         aes(x=Taxon, y=zenv.cp))+
  geom_linerange(aes(ymin=`5%`, ymax=`95%`), color="gray")+
  geom_hline(data=plot.chla.sumdat, aes(yintercept=cp, color=Assemblage), linetype="dashed", key_glyph="vline")+
  scale_color_manual(values=c("skyblue","tomato1","darkgreen"), name="Mean changepoint", labels=c("Benthic macroinvertebrates","Diatoms","Soft-bodied algae"))+
  geom_point(aes(fill=Assemblage), size=1.5, shape=21)+
  scale_fill_manual(values=c("skyblue","tomato1","darkgreen"), labels=c("Benthic macroinvertebrates","Diatoms","Soft-bodied algae"))+
  facet_wrap(~Group, scales="free")+
  theme_bw(base_size=10)+
  xlab("")+ylab("Changepoint")+
  ggtitle("Chlorophyll a (mg/m2)")+
  scale_size_continuous("Indicator value\n(z-score)")+
  theme(axis.text.y = element_text(face = "italic", size=6), panel.grid = element_blank())+
  guides(fill = guide_legend(override.aes=list(shape=21, size=3) , order=1),
         color = guide_legend(order=2))+
  coord_flip()
chla.titan.plot3_reflines
ggsave(chla.titan.plot3_reflines, filename="Figures_stratified/TITAN/chla.titan.plot3_reflines.jpg", dpi=300, width=6, height=tn.taxa*.1)

levels(plot.afdm.dat$Group)<-c("Neither", "Increaser taxa","Decreaser taxa")
plot.afdm.sumdat<-plot.afdm.sumdat %>%  mutate(Group=case_when(Group=="Increasers"~"Increaser taxa",T~"Decreaser taxa"))

afdm.titan.plot3<-
  ggplot(data=plot.afdm.dat %>%
           filter(filter>0 & zscore>min.z) %>%
           droplevels(),
         aes(x=Taxon, y=zenv.cp))+
  geom_linerange(aes(ymin=`5%`, ymax=`95%`, color="gray"),
                 key_glyph="path")+ 
  scale_color_manual(values='gray', name="", labels="95% confidence interval")+
  geom_point(aes(fill=Assemblage), size=1.5, shape=21)+
  scale_fill_manual(values=c("skyblue","tomato1","darkgreen"), labels=c("Benthic macroinvertebrates","Diatoms","Soft-bodied algae"))+
  facet_wrap(~Group, scales="free")+
  theme_bw(base_size=10)+
  xlab("")+ylab("Changepoint")+
  ggtitle("Ash-free dry mass (mg/cm2)")+
  scale_size_continuous("Indicator value\n(z-score)")+
  theme(axis.text.y = element_text(face = "italic", size=6), panel.grid = element_blank())+
  guides(fill = guide_legend(override.aes=list(shape=21, size=3) , order=1),
         color = guide_legend(order=2))+
  coord_flip()
afdm.titan.plot3
ggsave(afdm.titan.plot3, filename="Figures_stratified/TITAN/afdm.titan.plot_faceted.jpg", dpi=300, width=6, height=tn.taxa*.1)

afdm.titan.plot3_reflines<-
  ggplot(data=plot.afdm.dat %>%
           filter(filter>0 & zscore>min.z) %>%
           droplevels(),
         aes(x=Taxon, y=zenv.cp))+
  geom_linerange(aes(ymin=`5%`, ymax=`95%`), color="gray")+
  geom_hline(data=plot.afdm.sumdat, aes(yintercept=cp, color=Assemblage), linetype="dashed", key_glyph="vline")+
  scale_color_manual(values=c("skyblue","tomato1","darkgreen"), name="Mean changepoint", labels=c("Benthic macroinvertebrates","Diatoms","Soft-bodied algae"))+
  geom_point(aes(fill=Assemblage), size=1.5, shape=21)+
  scale_fill_manual(values=c("skyblue","tomato1","darkgreen"), labels=c("Benthic macroinvertebrates","Diatoms","Soft-bodied algae"))+
  facet_wrap(~Group, scales="free")+
  theme_bw(base_size=10)+
  xlab("")+ylab("Changepoint")+
  ggtitle("Ash-free dry mass (mg/cm2)")+
  scale_size_continuous("Indicator value\n(z-score)")+
  theme(axis.text.y = element_text(face = "italic", size=6), panel.grid = element_blank())+
  guides(fill = guide_legend(override.aes=list(shape=21, size=3) , order=1),
         color = guide_legend(order=2))+
  coord_flip()
afdm.titan.plot3_reflines
ggsave(afdm.titan.plot3_reflines, filename="Figures_stratified/TITAN/afdm.titan.plot3_reflines.jpg", dpi=300, width=6, height=tn.taxa*.1)


levels(plot.map.dat$Group)<-c("Neither", "Increaser taxa","Decreaser taxa")
plot.map.sumdat<-plot.map.sumdat %>%  mutate(Group=case_when(Group=="Increasers"~"Increaser taxa",T~"Decreaser taxa"))

map.titan.plot3<-
  ggplot(data=plot.map.dat %>%
           filter(filter>0 & zscore>min.z) %>%
           droplevels(),
         aes(x=Taxon, y=zenv.cp))+
  geom_linerange(aes(ymin=`5%`, ymax=`95%`, color="gray"),
                 key_glyph="path")+ 
  scale_color_manual(values='gray', name="", labels="95% confidence interval")+
  geom_point(aes(fill=Assemblage), size=1.5, shape=21)+
  scale_fill_manual(values=c("skyblue","tomato1","darkgreen"), labels=c("Benthic macroinvertebrates","Diatoms","Soft-bodied algae"))+
  facet_wrap(~Group, scales="free")+
  theme_bw(base_size=10)+
  xlab("")+ylab("Changepoint")+
  ggtitle("Percent macroalgal cover")+
  scale_size_continuous("Indicator value\n(z-score)")+
  theme(axis.text.y = element_text(face = "italic", size=6), panel.grid = element_blank())+
  guides(fill = guide_legend(override.aes=list(shape=21, size=3) , order=1),
         color = guide_legend(order=2))+
  coord_flip()
map.titan.plot3
ggsave(map.titan.plot3, filename="Figures_stratified/TITAN/map.titan.plot_faceted.jpg", dpi=300, width=6, height=tn.taxa*.1)


map.titan.plot3_reflines<-
  ggplot(data=plot.map.dat %>%
           filter(filter>0 & zscore>min.z) %>%
           droplevels(),
         aes(x=Taxon, y=zenv.cp))+
  geom_linerange(aes(ymin=`5%`, ymax=`95%`), color="gray")+
  geom_hline(data=plot.map.sumdat, aes(yintercept=cp, color=Assemblage), linetype="dashed", key_glyph="vline")+
  scale_color_manual(values=c("skyblue","tomato1","darkgreen"), name="Mean changepoint", labels=c("Benthic macroinvertebrates","Diatoms","Soft-bodied algae"))+
  geom_point(aes(fill=Assemblage), size=1.5, shape=21)+
  scale_fill_manual(values=c("skyblue","tomato1","darkgreen"), labels=c("Benthic macroinvertebrates","Diatoms","Soft-bodied algae"))+
  facet_wrap(~Group, scales="free")+
  theme_bw(base_size=10)+
  xlab("")+ylab("Changepoint")+
  ggtitle("Percent macroalgal cover")+
  scale_size_continuous("Indicator value\n(z-score)")+
  theme(axis.text.y = element_text(face = "italic", size=6), panel.grid = element_blank())+
  guides(fill = guide_legend(override.aes=list(shape=21, size=3) , order=1),
         color = guide_legend(order=2))+
  coord_flip()
map.titan.plot3_reflines
ggsave(map.titan.plot3_reflines, filename="Figures_stratified/TITAN/map.titan.plot3_reflines.jpg", dpi=300, width=6, height=tn.taxa*.1)



##############
mydata_sba<-bind_cols(sba_tax.cal %>%as.data.frame(),
                      alg_env.cal %>% as.data.frame() )

mydata_bmi<-bind_cols(bug_tax.cal %>%as.data.frame(),
                      bug_env.cal %>% as.data.frame() )

mytaxon<-"Cladophora"
mytaxon_tn_cp<-plot.tn.dat[mytaxon,"zenv.cp"]

increaser_example_plot<-ggplot(data=mydata_sba %>%
                                 mutate(Pres = case_when(Cladophora==0~0,T~1)),
                               aes(x=Nitrogen_Total_mgPerL, y=Pres))+
  geom_point(aes(y=Cladophora))+
  geom_vline(xintercept=mytaxon_tn_cp, linetype="dashed", color="red", size=1)+
  ylab("Relative abundance (points) or\n probability of occurrence (blue line)")+
  xlab("Total N (mg/L)")+
  geom_smooth(method = "glm", method.args = list(family = "binomial"), se = FALSE)+
  theme_bw(base_size = 14)
ggsave(increaser_example_plot+
         ggtitle(expression(italic(Cladophora)), subtitle=("A soft-bodied algal increaser taxon")), filename="Figures_stratified/TITAN/increaser_example_plot.jpg", dpi=300, width=6, height=5)  

plot.chla.dat %>% filter(filter==1 & Assemblage=="BMI")

mytaxon<-"Ephemerella"
mytaxon_chl_cp<-plot.chla.dat[mytaxon,"zenv.cp"]

decreaser_example_plot<-ggplot(data=mydata_bmi %>%
                                 mutate(Pres = case_when(Ephemerella==0~0,T~1)),
                               aes(x=Chlorophyll_a_mgPerm2, y=Pres))+
  geom_point(aes(y=Ephemerella))+
  geom_vline(xintercept=mytaxon_chl_cp, linetype="dashed", color="red", size=1)+
  ylab("Relative abundance (points) or\nprobability of occurrence (blue line)")+
  xlab("Chlorophyll a (mg/m2)")+
  geom_smooth(method = "glm", method.args = list(family = "binomial"), se = FALSE)+
  theme_bw(base_size = 14)
decreaser_example_plot
ggsave(decreaser_example_plot+
         ggtitle(expression(italic(Ephemerella)), 
                 subtitle=("A benthic macroinvertebrate decreaser taxon")), 
       filename="Figures_stratified/TITAN/decreaser_example_plot.jpg", dpi=300, width=6, height=5)
