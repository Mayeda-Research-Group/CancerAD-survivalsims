#------------------------------------------------------------------
# Title: Analysis of simulation output
# Author: Eleanor Hayes-Larson
#------------------------------------------------------------------
#------------------------------------------------------------------
# Loading packages, options and initializations #
#------------------------------------------------------------------
if (!require("pacman")) 
  install.packages("pacman", repos='http://cran.us.r-project.org')

p_load("tidyverse", "magrittr", "foreign",
       "deSolve","numDeriv", "ggplot2", "here", "rlang")

#------------------------------------------------------------------
# Analyze output from simulation 
#------------------------------------------------------------------

calib_plots<-function(immort, dcanmort, ddemmort, SdemIRR, SmortIRR, p_S, cancertype){

    #Run compartmodel function for input scenario
      modresults<-compartmodel(immortal=immort, 
                       diff.cancer.mort = dcanmort, 
                       diff.dem.mort = ddemmort, 
                       S1.dem.rateratio = SdemIRR, 
                       S1.mort.rateratio = SmortIRR, 
                       Sprev = p_S,
                       times=c(seq(from=0, to=40, by=1/12),10000),
                       type=cancertype,
                       pars=c())
      
    #Format data for plotting
    to.plot <- gather(modresults$out.data %>% as.data.frame,
                      "state","proportion",-time)
    to.plot <- to.plot[to.plot$time!=10000,] #drop longest time point
    
    out.data<-modresults$out.data

  
    ##########################################################    
    # Calculate cumulative survival calibration
    ##########################################################    
    
    
    #Plot cumulative survival (1-dead) in model and observed data
      #Take rows 2:7 because row 1 is NaN -- no mortality at start fo follow-up.
      if (cancertype=="all" | cancertype=="lung") {
        Life<-data.frame(time=Lifetables$time[2:7], Surv=Lifetables$Cum_surv_obs_all[2:7])
      } else if (cancertype=="breast") {
        Life<-data.frame(time=Lifetables$time[2:7], Surv=Lifetables$Cum_surv_obs_female[2:7])
      }else if (cancertype=="prostate") {
        Life<-data.frame(time=Lifetables$time[2:7], Surv=Lifetables$Cum_surv_obs_male[2:7])}
      
      cumsurv.calib.plot<-ggplot(to.plot[to.plot$state %in% c('DEAD'),],aes(time,1-proportion))+
                            geom_line(color="blue", size=1.25)+xlim(0,40)+ylab("Cumulative survival")+
                            geom_point(data=Life,aes(Life$time,Life$Surv), color="red", size=2)
    
    
    ##########################################################    
    # Calculate mortality rate ratios for dementia and cancer--should match inputs.
    ##########################################################    
    
      rate.results<-data.frame(modresults$rate.results)
      
      testratesum<-rate.results[,c("dC0D0S0","dC0D1S0","dC1D0S0","dC1D1S0","dC0D0S1","dC0D1S1","dC1D0S1","dC1D1S1", "dDEAD")]
      testratesum$sum<-rowSums(testratesum)
      sumrates<-round(mean(testratesum$sum), 10)
      
      #Get time- specific mortality rate ratio dementia from diffEQ model output
        #Used no cancer groups because cancer also affects mortality rate
        #stratify by S because S also affects mortality rate
        rate.nodem.S0.mort <- #weighted average rate of death from no dementia states
              rate.results$dDEADC0D0S0
            rate.dem.S0.mort<-  #weighted average rate of death from dementia states
              rate.results$dDEADC0D1S0
            rate.ratio.dem.S0.mort<- #take ratio, adjusted for proportion in each state
              (rate.dem.S0.mort/rate.nodem.S0.mort)*((out.data$C0D0S0)/(out.data$C0D1S0))  
            
        rate.nodem.S1.mort <- #weighted average rate of death from no dementia states
              (rate.results$dDEADC0D0S1*out.data$C0D0S1)/(out.data$C0D0S1)
            rate.dem.S1.mort<-  #weighted average rate of death from dementia states
              (rate.results$dDEADC0D1S1*out.data$C0D1S1)/(out.data$C0D1S1)
            rate.ratio.dem.S1.mort<- #take ratio, adjusted for proportion in each state
              (rate.dem.S1.mort/rate.nodem.S1.mort)*((out.data$C0D0S1)/(out.data$C0D1S1))  
            
        #Compare to observed HRs input to parameterize model
          rate.ratio.dem.mort.check<-data.frame(age=AD_HR_obs$Age, 
                                                Sim.S0.HR=c(mean(rate.ratio.dem.S0.mort[2:60]), #Start at 2 because rate ratio undefine in first month 
                                                            mean(rate.ratio.dem.S0.mort[61:120]), 
                                                            mean(rate.ratio.dem.S0.mort[121:180]), 
                                                            mean(rate.ratio.dem.S0.mort[181:240]),
                                                            mean(rate.ratio.dem.S0.mort[241:300]),
                                                            mean(rate.ratio.dem.S0.mort[301:480])),
                                                Sim.S1.HR=c(mean(rate.ratio.dem.S1.mort[2:60]), #Start at 2 because rate ratio undefine in first month 
                                                            mean(rate.ratio.dem.S1.mort[61:120]), 
                                                            mean(rate.ratio.dem.S1.mort[121:180]), 
                                                            mean(rate.ratio.dem.S1.mort[181:240]),
                                                            mean(rate.ratio.dem.S1.mort[241:300]),
                                                            mean(rate.ratio.dem.S1.mort[301:480])),                                                  
                                                Obs.HR=AD_HR_obs$Dementia.HR)
          rate.ratio.dem.mort.check           
            
    
    
      #Get time-specific S mortality rate ratio from diffEQ model output
          #Use no dementia groups because dementia affects mortality
          rate.cancer.S0.mort<-rate.results$dDEADC1D0S0 #rate in cancer+ without S
          rate.cancer.S1.mort <- rate.results$dDEADC1D0S1 #rate in cancer+ with S
        
          rate.ratio.S1.mort<- #take ratio, adjust for proprotion in each state.
            (rate.cancer.S1.mort/rate.cancer.S0.mort)*((out.data$C1D0S0)/(out.data$C1D0S1))  
          
          #Compare to IRRs input to parameterize model
          rate.ratio.S1.mort.check<-data.frame(Sim_SmortIRR=mean(rate.ratio.S1.mort[2:480]), 
                                               Input_SmortIRR=SmortIRR)
        
          rate.ratio.S1.mort.check           
        
      
      #Check relative survival over 5 years to check cancer mortality rate ratio
        #Get month-specific survival risk ratio for cancer
        #Use no dementia group because dementia affects mortality
          num.dead.nocancer<-#diff in #dead without cancer at t and t-1 month
            out.data$DEADC0D0S0-lag(out.data$DEADC0D0S0) 
          risk.surv.nocancer<- #surv risk = 1-risk(dead), where risk(dead)= #dead/#at risk
            1-(num.dead.nocancer/(lag(out.data$C0D0S0)))
          
          num.dead.cancer<-#diff in #dead with cancer at t and t-1 month
            out.data$DEADC1D0S0-lag(out.data$DEADC1D0S0)
          risk.surv.cancer<- #surv risk = 1-risk(dead), where risk(dead)= #dead/#at risk
            1-(num.dead.cancer/(lag(out.data$C1D0S0)))
          
          mort.riskratio.cancer<-risk.surv.cancer/risk.surv.nocancer #monthly surv risk ratio
          
      #Convert to 5 year relative survival by multiplying monthly rel surv over 5 months.
          rel.surv.6574<-prod(mort.riskratio.cancer[61:120]) 
                #I picked years 5-10 for the check because years 0-5 the risk ratio is undefined in months 1-2. 
          rel.surv.75up<-prod(mort.riskratio.cancer[(12*15+1):(12*20)]) 
                #I picked years 15-20 to check this relative survival--could have picked a different window.
          
          #Compare to observed 5-year relative survivl from SEER input to parameterize model
          cancer.rel.surv.check<-cbind(data.frame(Age=SEER_rel.surv$Age, Cancertype=cancertype,
                                            Sim_relsurv=c(rel.surv.6574,rel.surv.75up)),
                                       SEER_rel.surv[,2:5]) 
       
          cancer.rel.surv.check
          
        #Check dem rate ratio with/without S
          
          rate.S0.dem<-rate.results$dDEMNOSNOC #rate in cancer+ without S
          rate.S1.dem <- rate.results$dDEMWITHSNOC #rate in cancer+ with S
          
          rate.ratio.S1.dem<- #take ratio, adjust for proprotion in each state.
            (rate.S1.dem/rate.S0.dem)*((out.data$C0D0S0)/(out.data$C0D0S1))  
         
          #Compare to IRRs input to parameterize model
          rate.ratio.S1.dem.check<-data.frame(Sim_SdemIRR=mean(rate.ratio.S1.dem[2:480]), 
                                               Input_SdemIRR=SdemIRR)
          
          rate.ratio.S1.dem.check      
    ##########################################################    
    # Calculate dementia incidence rates to check calibration
    ##########################################################    
    
      #Calculate dementia incidence rates in simulation
      dem_inc_mod<-rep(NA,6)
            
        
      
      for (i in seq(from=5, to=30, by=5)){
        #This is a loop that goes in 5 year intervals to add up incident dementia and person time
        if (i<30){ #Need a longer interval for 30+ years, because ACT data has 5 year age bands and then 90+
          
          ncases.dem.interval<- #diff in number of ever-dementia at t and t-5 years
            out.data$DEMENTIA[out.data$time==i]-out.data$DEMENTIA[out.data$time==(i-5)]
          
          ncases.dem.months<-rep(NA,60)
          ndeathsnodem.months<-rep(NA,60)
          PT.dem.death<-0
           
          #Add up person time in each month. 
          #This is a nested loop that goes in monthly (1/12 yearly) intervals to add up person time
    
          for (j in seq(from=1/12, to=5, by=1/12)){
            ncases.dem.months[j*12] <- #tally new dementia cases in month
              out.data$DEMENTIA[round(out.data$time,2)==round(i-(5-j),2)]-out.data$DEMENTIA[round(out.data$time,2)==round((i-(5-j)-1/12),2)]
            ndeathsnodem.months[j*12] <- #tally number of deaths from people at risk for dementia in month
                out.data$DEADC0D0S0[round(out.data$time,2)==round(i-(5-j),2)]-out.data$DEADC0D0S0[round(out.data$time,2)==round((i-(5-j)-1/12),2)] +
                out.data$DEADC1D0S0[round(out.data$time,2)==round(i-(5-j),2)]-out.data$DEADC1D0S0[round(out.data$time,2)==round((i-(5-j)-1/12),2)] +
                out.data$DEADC0D0S1[round(out.data$time,2)==round(i-(5-j),2)]-out.data$DEADC0D0S1[round(out.data$time,2)==round((i-(5-j)-1/12),2)] +
                out.data$DEADC1D0S1[round(out.data$time,2)==round(i-(5-j),2)]-out.data$DEADC1D0S1[round(out.data$time,2)==round((i-(5-j)-1/12),2)]
            
            PT.dem.death<-PT.dem.death+((j-(0.5/12))*(ndeathsnodem.months[j*12]+ncases.dem.months[j*12]))
            #People who get dementia or die in a month get 1/2*1/12 year of person time
          }
          
          #Rate is cases/persontime. People surviving to next interval get 5 years, + PT for dem/death cases
          dem_inc_mod[i/5]<-ncases.dem.interval/(5*(out.data$C0D0S0[round(out.data$time,2)==round(i,2)] +
                                                  out.data$C1D0S0[round(out.data$time,2)==round(i,2)] +
                                                  out.data$C0D0S1[round(out.data$time,2)==round(i,2)] +
                                                  out.data$C1D0S1[round(out.data$time,2)==round(i,2)] ) +
                                               PT.dem.death)
        }
        else { #For 30+ years, because ACT data has all 90+, go until year 40. 
                  #Same logic as above but for 10 year interval
          ncases.dem.interval<-out.data$DEMENTIA[out.data$time==(i+10)]-out.data$DEMENTIA[out.data$time==(i-5)]
          ncases.dem.months<-rep(NA,15*12)
          ndeathsnodem.months<-rep(NA,15*12)
          PT.dem.death<-0
          
          #Add up person time in each year
          for (j in seq(from=1/12, to=15, by=1/12)){
            ncases.dem.months[j*12]<-out.data$DEMENTIA[round(out.data$time,2)==round((i+10)-(15-j),2)]-out.data$DEMENTIA[round(out.data$time,2)==round(((i+10)-(15-j)-1/12),2)]
            ndeathsnodem.months[j*12]<-  out.data$DEADC0D0S0[round(out.data$time,2)==round((i+10)-(15-j),2)] - out.data$DEADC0D0S0[round(out.data$time,2)==round(((i+10)-(15-j)-1/12),2)] +
                                         out.data$DEADC1D0S0[round(out.data$time,2)==round((i+10)-(15-j),2)] - out.data$DEADC1D0S0[round(out.data$time,2)==round(((i+10)-(15-j)-1/12),2)] + 
                                         out.data$DEADC0D0S1[round(out.data$time,2)==round((i+10)-(15-j),2)] - out.data$DEADC0D0S1[round(out.data$time,2)==round(((i+10)-(15-j)-1/12),2)] +
                                         out.data$DEADC1D0S1[round(out.data$time,2)==round((i+10)-(15-j),2)] - out.data$DEADC1D0S1[round(out.data$time,2)==round(((i+10)-(15-j)-1/12),2)] 
                                    
            PT.dem.death<-PT.dem.death+((j-1/24)*(ndeathsnodem.months[j*12]+ncases.dem.months[j*12]))
          }
          
          #Rate is cases/persontime 
          dem_inc_mod[i/5]<-ncases.dem.interval/(15*(out.data$C0D0S0[out.data$time==(i+10)] +
                                                     out.data$C1D0S0[out.data$time==(i+10)] +
                                                     out.data$C0D0S1[out.data$time==(i+10)] +
                                                     out.data$C1D0S1[out.data$time==(i+10)] ) +
                                                 PT.dem.death)    
        }}
      
      
      #Store incident rates over time
      dem_inc_mod<-data.frame(time=seq(from=5, to=30, by=5), dem_inc_rate=dem_inc_mod)
      
      #Plot dementia incidence rate in model and observed data, pulling in correct sex AD rates
      if (cancertype=="all" | cancertype=="lung") {AD_obs<-data.frame(time=AD_inc_obs$time[2:7], Rate=AD_inc_obs$Rate_all[2:7])
          } else if (cancertype=="breast") {AD_obs<-data.frame(time=AD_inc_obs$time[2:7], Rate=AD_inc_obs$Rate_F[2:7])
              }else if (cancertype=="prostate") {AD_obs<-data.frame(time=AD_inc_obs$time[2:7], Rate=AD_inc_obs$Rate_M[2:7])}
      
      dem.calib.plot<-ggplot(data=AD_obs,aes(time,Rate))+
        geom_point(color="red", size=2)+xlim(0,35)+ylim(0,0.125)+ylab("Dementia incidence rate")+
        geom_point(data=dem_inc_mod[1:6,],aes(dem_inc_mod$time[1:6],dem_inc_mod$dem_inc_rate[1:6]), color="blue", size=2)
    
    ##########################################################    
    # Calculate cancer incidence rates to check calibration
    ##########################################################    
        #Calculate cancer incidence rate in simulation
        #All of this works the same ways as the dementia calculations.  
        cancer_inc_mod<-rep(NA,5)
        
        for (i in seq(from=5, to=25, by=5)){
          if (i<25){
            ncases.cancer.interval<-out.data$CANCER[out.data$time==i]-out.data$CANCER[out.data$time==(i-5)]
            
            ncases.cancer.months<-rep(NA,60)
            ndeathsnocancer.months<-rep(NA,60)
            PT.cancer.death<-0
            
            #Add up person time in each month
            for (j in seq(from=1/12, to=5, by=1/12)){
              ncases.cancer.months[j*12] <- out.data$CANCER[round(out.data$time,2)==round(i-(5-j),2)]-out.data$CANCER[round(out.data$time,2)==round((i-(5-j)-1/12),2)]
              ndeathsnocancer.months[j*12] <- out.data$DEADC0D0S0[round(out.data$time,2)==round(i-(5-j),2)]-out.data$DEADC0D0S0[round(out.data$time,2)==round((i-(5-j)-1/12),2)] +
                                              out.data$DEADC0D1S0[round(out.data$time,2)==round(i-(5-j),2)]-out.data$DEADC0D1S0[round(out.data$time,2)==round((i-(5-j)-1/12),2)] +
                                              out.data$DEADC0D0S1[round(out.data$time,2)==round(i-(5-j),2)]-out.data$DEADC0D0S1[round(out.data$time,2)==round((i-(5-j)-1/12),2)] +
                                              out.data$DEADC0D1S1[round(out.data$time,2)==round(i-(5-j),2)]-out.data$DEADC0D1S1[round(out.data$time,2)==round((i-(5-j)-1/12),2)]
              PT.cancer.death<-PT.cancer.death+((j-(0.5/12))*(ndeathsnocancer.months[j*12]+ncases.cancer.months[j*12]))
            }
            
            #Rate is cases/persontime 
            cancer_inc_mod[i/5]<-ncases.cancer.interval/(5*(out.data$C0D0S0[round(out.data$time,2)==round(i,2)] +
                                                        out.data$C0D1S0[round(out.data$time,2)==round(i,2)] +
                                                        out.data$C0D0S1[round(out.data$time,2)==round(i,2)] +
                                                        out.data$C0D1S1[round(out.data$time,2)==round(i,2)] ) +
                                                     PT.cancer.death)
          }
          else {
            ncases.cancer.interval<-out.data$CANCER[out.data$time==(i+15)]-out.data$CANCER[out.data$time==(i-5)]
            ncases.cancer.months<-rep(NA,20*12)
            ndeathsnocancer.months<-rep(NA,20*12)
            PT.cancer.death<-0
        
            #Add up person time in each year
            for (j in seq(from=1/12, to=20, by=1/12)){
              ncases.cancer.months[j*12]<-out.data$CANCER[round(out.data$time,2)==round((i+15)-(20-j),2)]-out.data$CANCER[round(out.data$time,2)==round(((i+15)-(20-j)-1/12),2)]
              ndeathsnocancer.months[j*12]<-  out.data$DEADC0D0S0[round(out.data$time,2)==round((i+15)-(20-j),2)] - out.data$DEADC0D0S0[round(out.data$time,2)==round(((i+15)-(20-j)-1/12),2)] +
                                              out.data$DEADC0D1S0[round(out.data$time,2)==round((i+15)-(20-j),2)] - out.data$DEADC0D1S0[round(out.data$time,2)==round(((i+15)-(20-j)-1/12),2)] + 
                                              out.data$DEADC0D0S1[round(out.data$time,2)==round((i+15)-(20-j),2)] - out.data$DEADC0D0S1[round(out.data$time,2)==round(((i+15)-(20-j)-1/12),2)] +
                                              out.data$DEADC0D1S1[round(out.data$time,2)==round((i+15)-(20-j),2)] - out.data$DEADC0D1S1[round(out.data$time,2)==round(((i+15)-(20-j)-1/12),2)] 
              PT.cancer.death<-PT.cancer.death+((j-1/24)*(ndeathsnocancer.months[j*12]+ncases.cancer.months[j*12]))
            }
            
            #Rate is cases/persontime 
            cancer_inc_mod[i/5]<-ncases.cancer.interval/(15*(out.data$C0D0S0[out.data$time==(i+15)] +
                                                         out.data$C0D1S0[out.data$time==(i+15)] +
                                                         out.data$C0D0S1[out.data$time==(i+15)] +
                                                         out.data$C0D1S1[out.data$time==(i+15)] ) +
                                                     PT.cancer.death)    
          }}
        
        
        cancer_inc_mod<-data.frame(time=seq(from=5, to=25, by=5), cancer_inc_rate=cancer_inc_mod)
        
        
        #Plot cancer incidence rate in model and observed data
        cancer.calib.plot<-ggplot(data=SEER_inc_obs,aes(time, eval(parse_expr(paste("Rate_",cancertype, sep="")))))+
                            geom_point(color="red", size=2)+xlim(0,30)+ylab(paste("Rate of",cancertype, "cancer",sep=" "))+
                            geom_point(data=cancer_inc_mod[1:5,],aes(cancer_inc_mod$time[1:5],cancer_inc_mod$cancer_inc_rate[1:5]), color="blue", size=2)
    
    ##########################################################    
    # Calib_plots function output
    ##########################################################    
        
    return(list(sumrates=sumrates,
                rate.ratio.dem.mort.check=rate.ratio.dem.mort.check,
                rate.ratio.S1.mort.check=rate.ratio.S1.mort.check,
                cancer.rel.surv.check=cancer.rel.surv.check, 
                rate.ratio.S1.dem.check=rate.ratio.S1.dem.check,
                cumsurv.calib.plot=cumsurv.calib.plot, 
                dem.calib.plot=dem.calib.plot, 
                cancer.calib.plot=cancer.calib.plot))

}

#calib_plots function returns the rate ratios and 3 plots to confirm that the 
    # model is correctly calibrated to total dementia, cancer, and mortality incidence
#this function needs to be run for each set of the parameter values for which the model is calibrated, and checked.




# This is a user-interactive loop to manually inspect the calibration plots 
  # for each scenario. The loop produces 3 plots per scenario, and I broke the 
  # loop when a scenario did not look calibrated. 
# For code review purposes, all the scenarios are already calibrated. 
  # I have commented out this loop and all the calibration plots can be inspected 
  # together in another script.


# library(here)
# source(here("Code", "Multistate model","Rate inputs.R"))
# # 
# S1demvals<-c(1,0.9,0.7, 0.5,0.3)
# S1mortvals<-c(1,0.9,0.7, 0.5,0.3)
# Sprevvals<-c(0.5,0.4,0.3,0.2,0.1)
# typevals<-c("all", "lung", "breast", "prostate")
# 
# args<-expand.grid(S1demvals=S1demvals, S1mortvals=S1mortvals, Sprevvals=Sprevvals, typevals=typevals)
# args_short<-args[args$Sprevvals!=0 | (args$Sprevvals==0 & args$S1demvals==1 & args$S1mortvals==1),]
# row.names(args_short) <- NULL
# 
#source(here("Code", "Multistate model","Model.R"))
# 
# 
# for (i in 1:nrow(args_short)){
#     checks<-calib_plots(immort=0, dcanmort=1, ddemmort=1, SdemIRR=args_short$S1demvals[i],
#                 SmortIRR=args_short$S1mortvals[i], p_S=args_short$Sprevvals[i],
#                 cancertype=args_short$typevals[i])
# 
#           print(checks$cumsurv.calib.plot)
#           print(checks$dem.calib.plot)
#           print(checks$cancer.calib.plot)
# 
#        print(paste("Plots for ",args_short$typevals[i], " cancer, S mortality IRR = ",
#                    args_short$S1mortvals[i],", S dementia IRR = ", args_short$S1demvals[i],
#                    " and S prevalence = ", args_short$Sprevvals[i], ".",sep=""))
#         question <- readline(prompt="Proceed to calibration plots for next scenario?
#                                       Type 'yes' or 'no'. ")
#         if (question=="yes"){
#           next
#         } else {
#           print(paste("Ending loop at args_short row ", i, ".", sep=""))
#           break
#     }
#   }
# 
# 
#Show 1 plot for immortal cohort.
# checks<-calib_plots(immort=0, dcanmort=1, ddemmort=1, SdemIRR=.3, SmortIRR=.3, p_S=0.3, cancertype="all")
# print(checks$cumsurv.calib.plot)
# print(checks$dem.calib.plot)
# print(checks$cancer.calib.plot)
