#------------------------------------------------------------------
# Title: Rate inputs
# Author: Eleanor Hayes-Larson
# Purpose: This program loads in data to be used in calibrating the 
#     multistate model.
#------------------------------------------------------------------

#------------------------------------------------------------------
# Loading packages, options and initializations #
#------------------------------------------------------------------
if (!require("pacman")) 
  install.packages("pacman", repos='http://cran.us.r-project.org')

p_load("tidyverse", "magrittr", "foreign",
       "deSolve","numDeriv", "ggplot2", "here")

#------------------------------------------------------------------
# Defining mortality rates #
#------------------------------------------------------------------
      #Make overall mortality rate a function of time
      #--- Load survival data from lifetables ---
      Lifetables<-read.csv(file=here("Calibration data","Lifetables.csv"),
                           header=T, sep=",", stringsAsFactors=F)
      
      Lifetables<-Lifetables[,c("Age", "Allrace_all", "Allrace_M", "Allrace_F")]

      #--- Calculate conditional mortality rate and cumulative survival ---
      #Mort rate is # deaths/person-time (5 years for survivors and 2.5 for deaths in interval)
      Lifetables$Mortality_obs_all<- 
        (lag(Lifetables$Allrace_all)-Lifetables$Allrace_all)/
        (5*Lifetables$Allrace_all+2.5*(lag(Lifetables$Allrace_all)-Lifetables$Allrace_all))
      
      Lifetables$Cum_surv_obs_all<-(Lifetables$Allrace_all/Lifetables$Allrace_all[1])
      
      #Same calculation for men only
      Lifetables$Mortality_obs_male<-
        (lag(Lifetables$Allrace_M)-Lifetables$Allrace_M)/
        (5*Lifetables$Allrace_M+2.5*(lag(Lifetables$Allrace_M)-Lifetables$Allrace_M))
      
      Lifetables$Cum_surv_obs_male<-(Lifetables$Allrace_M/Lifetables$Allrace_M[1])
      
      #Same calculation for women only
      Lifetables$Mortality_obs_female<-
        (lag(Lifetables$Allrace_F)-Lifetables$Allrace_F)/
        (5*Lifetables$Allrace_F+2.5*(lag(Lifetables$Allrace_F)-Lifetables$Allrace_F))
      
      Lifetables$Cum_surv_obs_female<-(Lifetables$Allrace_F/Lifetables$Allrace_F[1])
      
      Lifetables$time<-seq(from=0, to=40, by=5)
      Lifetables
      

#------------------------------------------------------------------
# Defining cancer incidence rates #
#------------------------------------------------------------------
    #make cancer incidence rate a function of time (study time = age, starting at age 65) 
      ages<-c("65-69","70-74","75-79","80-84","85+")
      SEER_inc_obs<-data.frame(Age = ages, time = seq(from=5, to=25, by=5), 
                                Rate_all = NA, 
                                Rate_lung = NA, 
                                Rate_breast = NA, 
                                Rate_prostate=NA) 
      
      types<-c("all", "lung", "breast", "prostate")
      sexes=c("Both Sexes","Both Sexes", "Female", "Male")
      
      for (i in 1:4){
        SEER_inc <- read.csv(here("Calibration data", 
                                paste("SEER ",types[i], " cancer incidence.csv", sep="")),
                                skip=3, stringsAsFactors = F, header=T)
        SEER_inc_obs[,i+2]<- #Pull rates from SEER, and convert from rate per 100k.
          if (i!=2) {as.numeric(SEER_inc$Rate.per.100.000[
            ((SEER_inc$Age=="65-69" | 
                SEER_inc$Age=="70-74" | 
                SEER_inc$Age=="75-79" | 
                SEER_inc$Age=="80-84" |
                SEER_inc$Age=="85+") &
               SEER_inc$Sex==sexes[i] & 
               SEER_inc$Race.Ethnicity=="All Races (includes Hispanic)" &
               SEER_inc$Rate.Type=="Delay-adjusted Rates")])/100000 
          } else {
            (as.numeric(SEER_inc$Rate.per.100.000[
              ((SEER_inc$Age=="65-69" | 
                  SEER_inc$Age=="70-74" | 
                  SEER_inc$Age=="75-79" | 
                  SEER_inc$Age=="80-84" |
                  SEER_inc$Age=="85+") &
                 SEER_inc$Sex==sexes[i] & 
                 SEER_inc$Race.Ethnicity=="All Races (includes Hispanic)")])/100000)[6:10]
          } 
        #Lung cancer SEER data doesn't
        # have a label for observed vs. 
        # delay-adjusted, but larger rates (rows 6-10) 
        # are delay-adjusted, so I take these.
      }
      SEER_inc_obs

    
#------------------------------------------------------------------
# Defining cancer relative survival #
#------------------------------------------------------------------
    
      #make cancer relative survival a function of time (study time = age, starting at age 65) 
      ages2<-c("65-74","75+")
      SEER_rel.surv<-data.frame(Age = ages2,  
                               Rel.surv_all = NA, 
                               Rel.surv_lung = NA, 
                               Rel.surv_breast = NA, 
                               Rel.surv_prostate=NA) 
      
      for (i in 1:4){
        SEER_surv <- read.csv(here("Calibration data", 
                                  paste("SEER ",types[i], " cancer survival.csv", sep="")),
                             skip=3, stringsAsFactors = F, header=T)
        SEER_rel.surv[,i+1]<-as.numeric(SEER_surv$Relative.Survival.Rate....[
            ((SEER_surv$Age=="Ages 65-74" | 
                SEER_surv$Age=="Ages 75+") &
               SEER_surv$Sex==sexes[i] & 
               SEER_surv$Race.Ethnicity=="All Races (includes Hispanic)" &
               SEER_surv$Stage.at.Diagnosis=="All Stages")])/100 #convert to %
      }
      SEER_rel.surv
      
#------------------------------------------------------------------
# Defining dementia incidence rates #
#------------------------------------------------------------------
    #make dementia incidence rate a function of time (study time = age, starting at age 65)
    #Load ACT study data--use all-cause dementia, not just AD.
    AD_inc_obs<-read.csv(file=here("Calibration data","AD_inc_Tom_etal.csv"),
                         header=T, sep=",", stringsAsFactors=F)
    names(AD_inc_obs)[names(AD_inc_obs) == "ï..Age"] <- "Age"
    AD_inc_obs$Rate_all<-AD_inc_obs$Rate_all/1000 #Convert from rate/1000
    AD_inc_obs$Rate_M<-AD_inc_obs$Rate_M/1000
    AD_inc_obs$Rate_F<-AD_inc_obs$Rate_F/1000
    AD_inc_obs$time<-seq(from=0, to=40, by=5)
    AD_inc_obs
    
    
    
#------------------------------------------------------------------
# Defining dementia mortality hazard ratios  #
#------------------------------------------------------------------
    #Load Mayeda 2017 study data--use results for whites
    AD_HR_obs<-read.csv(file=here("Calibration data","Dem_mort_Mayeda_etal.csv"),
                         header=T, sep=",", stringsAsFactors=F)
   
    AD_HR_obs
    
#------------------------------------------------------------------
# Dropping unneeded objects
#------------------------------------------------------------------
    
    rm(list=setdiff(ls(), c("Lifetables","SEER_inc_obs","SEER_rel.surv","AD_inc_obs","AD_HR_obs",lsf.str())))
    
    