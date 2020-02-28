#------------------------------------------------------------------
# Title:Compartment model function
# Author: Eleanor Hayes-Larson
#------------------------------------------------------------------

#Input parameters of simulation
# immortal #This is a switch to make cohort immortal (when immortal==1)
# diff.cancer.mort #This is a switch to make mortality higher in those with cancer
# diff.dem.mort #This is a switch to make mortality higher in those with dementia
# S1.dem.rateratio #Annual rate ratio for protective effect of S charactistic on dementia
# S1.mort.rateratio #Annual rate ratio for protective effect of S characteristic on mortality in cancer+
# times  #Time steps for model in years. Request long time interval to ensure steady state (100% dead)
# pars #This is empty because all rates defined separately as functions of time, 
  # so no time-constant rates to input. But model won't run without it. 

compartmodel<-function(immortal,diff.cancer.mort,diff.dem.mort,S1.dem.rateratio, 
                       S1.mort.rateratio, Sprev, times, type, pars){

#------------------------------------------------------------------
# Part 1: Defining sources and calibration parameters for each cancer type
  #These are calibrations that I did manually 
  #They can be checked using the "Calibration checks" Rmd file.
#------------------------------------------------------------------
    if (type=="all"){Mortality<-Lifetables$Mortality_obs_all  
                  Rel_surv<-SEER_rel.surv$Rel.surv_all 
                  Cancer_inc<-SEER_inc_obs$Rate_all
                  Dem_inc<-AD_inc_obs$Rate_all
                  if (S1.mort.rateratio>0.7){
                      Mort_cal<-c(.9,0.75,0.7,0.7,0.7,0.7,0.7,0.6) 
                      demrate.adj.S_IRR0.3<-list(prev0.5=c(rep(1.5,3),1.65, 1.85, 2.25),
                                                 prev0.4=c(rep(1.45,3),1.5, 1.6, 1.95),
                                                 prev0.3=c(rep(1.3,3),1.35, 1.42, 1.65),
                                                 prev0.2=c(rep(1.15,3),1.18, 1.25, 1.39),
                                                 prev0.1=c(c(rep(1.06,3),1.1, 1.12, 1.17)))
                      demrate.adj.S_IRR0.5<-list(prev0.5=c(rep(1.3,3),1.37, 1.42,1.53),
                                                 prev0.4=c(rep(1.3,3),1.32, 1.34,1.4),
                                                 prev0.3=c(rep(1.2,3),1.22, 1.24,1.29),
                                                 prev0.2=c(rep(1.1,3),1.15, 1.14, 1.19),
                                                 prev0.1=c(c(rep(1.05,4),1.07, 1.09)))
                      demrate.adj.S_IRR0.7<-list(prev0.5=c(rep(1.2,3),1.17,1.2,1.2),
                                                 prev0.4=c(rep(1.2,3),1.17,1.18,1.18),
                                                 prev0.3=c(rep(1.15,3),1.11,1.11,1.12),
                                                 prev0.2=c(rep(1.08,3),1.08, 1.08, 1.09),
                                                 prev0.1=c(c(rep(1.04,5),1.03)))
                      demrate.adj.S_IRR0.9<-list(prev0.5=c(rep(1.05,5),1.05),
                                                 prev0.4=c(rep(1.05,5),1.05),
                                                 prev0.3=c(rep(1.03,5),1.03),
                                                 prev0.2=c(rep(1,3),1.03, 1.03, 1.03),
                                                 prev0.1=rep(1,6))
                   } else if (S1.mort.rateratio>0.5 & S1.mort.rateratio<=0.7){
                      Mort_cal<-c(.9,0.75,0.7,0.7,0.7,0.7,0.7,0.6) 
                      demrate.adj.S_IRR0.3<-list(prev0.5=c(rep(1.5,2),1.55,1.67, 1.9, 2.35),
                                                 prev0.4=c(rep(1.5,2),1.5, 1.56, 1.67, 2.07),
                                                 prev0.3=c(rep(1.3,2),1.33, 1.37, 1.47, 1.75),
                                                 prev0.2=c(rep(1.15,3),1.22, 1.27, 1.45),
                                                 prev0.1=c(c(rep(1.06,3),1.08, 1.12, 1.2)))
                      demrate.adj.S_IRR0.5<-list(prev0.5=c(rep(1.3,3),1.38, 1.45,1.57),
                                                 prev0.4=c(rep(1.3,3),1.32, 1.35, 1.45),
                                                 prev0.3=c(rep(1.22,3),1.22, 1.23, 1.33),
                                                 prev0.2=c(rep(1.12,3),1.15, 1.15, 1.2),
                                                 prev0.1=c(c(rep(1.05,4),1.07, 1.1)))
                      demrate.adj.S_IRR0.7<-list(prev0.5=c(rep(1.2,3),1.17,1.2,1.23),
                                                 prev0.4=c(rep(1.2,3),1.17,1.18,1.19),
                                                 prev0.3=c(rep(1.12,3),1.12,1.12,1.14),
                                                 prev0.2=c(rep(1.07,3),1.08, 1.08, 1.1),
                                                 prev0.1=c(c(rep(1.04,5),1.05)))
                      demrate.adj.S_IRR0.9<-list(prev0.5=c(rep(1.05,5),1.05),
                                                 prev0.4=c(rep(1.05,5),1.05),
                                                 prev0.3=c(rep(1.03,3),1.03, 1.03, 1.03),
                                                 prev0.2=c(rep(1,3),1.02, 1.02, 1.02),
                                                 prev0.1=rep(1,6))
                    } else if (S1.mort.rateratio>0.3 & S1.mort.rateratio<=0.5){ 
                      Mort_cal<-c(.9,0.8,0.75,0.75,0.75,0.75,0.75,0.65) 
                      demrate.adj.S_IRR0.3<-list(prev0.5=c(rep(1.52,3), 1.73, 1.98, 2.5),
                                                 prev0.4=c(rep(1.5,3), 1.6, 1.75, 2.25),
                                                 prev0.3=c(rep(1.32,3), 1.38, 1.5, 1.9),
                                                 prev0.2=c(rep(1.2,3),1.25, 1.3, 1.56),
                                                 prev0.1=c(c(rep(1.07,3),1.09, 1.13, 1.25)))
                      demrate.adj.S_IRR0.5<-list(prev0.5=c(rep(1.35,3),1.43, 1.48,1.62),
                                                 prev0.4=c(rep(1.33,3),1.34, 1.38,1.51),
                                                 prev0.3=c(rep(1.23,3),1.24, 1.28,1.38),
                                                 prev0.2=c(rep(1.15,3),1.15, 1.18, 1.26),
                                                 prev0.1=c(c(rep(1.05,3),1.07, 1.09, 1.13)))
                      demrate.adj.S_IRR0.7<-list(prev0.5=c(rep(1.2,3),1.19,1.22,1.25),
                                                 prev0.4=c(rep(1.2,3),1.18,1.20,1.21),
                                                 prev0.3=c(rep(1.12,3),1.13,1.14,1.16),
                                                 prev0.2=c(rep(1.1,3),1.1, 1.1, 1.11),
                                                 prev0.1=c(c(rep(1.05,5),1.06)))
                      demrate.adj.S_IRR0.9<-list(prev0.5=c(rep(1.05,5),1.07),
                                                 prev0.4=c(rep(1.05,5),1.07),
                                                 prev0.3=c(rep(1.03,5),1.05),
                                                 prev0.2=c(rep(1.01,5),1.03),
                                                 prev0.1=c(rep(1,4), 1.02, 1.02))
                    } else{ 
                      Mort_cal<-c(.9,0.8,0.75,0.75,0.8,0.75,0.75,0.65) 
                      demrate.adj.S_IRR0.3<-list(prev0.5=c(rep(1.52,2), 1.55, 1.78, 2.15, 2.75),
                                                   prev0.4=c(rep(1.48,2), 1.45, 1.6, 1.87, 2.53),
                                                   prev0.3=c(rep(1.38,2), 1.34, 1.42, 1.61, 2.19),
                                                   prev0.2=c(rep(1.2,3),1.22, 1.35, 1.79),
                                                   prev0.1=c(c(rep(1.07,3),1.09, 1.15, 1.39)))
                        demrate.adj.S_IRR0.5<-list(prev0.5=c(rep(1.35,3),1.43, 1.55,1.7),
                                                   prev0.4=c(rep(1.33,3),1.35, 1.42,1.63),
                                                   prev0.3=c(rep(1.23,3),1.22, 1.32,1.49),
                                                   prev0.2=c(rep(1.13,3),1.15, 1.21,1.36),
                                                   prev0.1=c(c(rep(1.05,3),1.08, 1.1, 1.19)))
                        demrate.adj.S_IRR0.7<-list(prev0.5=c(rep(1.2,3),1.19,1.22,1.29),
                                                   prev0.4=c(rep(1.2,3),1.18,1.2,1.25),
                                                   prev0.3=c(rep(1.15,3),1.15,1.15,1.21),
                                                   prev0.2=c(rep(1.1,3),1.1,1.1,1.15),
                                                   prev0.1=c(c(rep(1.05,5),1.09)))
                        demrate.adj.S_IRR0.9<-list(prev0.5=c(rep(1.05,5),1.07),
                                                   prev0.4=c(rep(1.05,5),1.07),
                                                   prev0.3=c(rep(1.03,5),1.05),
                                                   prev0.2=c(rep(1.02,5),1.03),
                                                   prev0.1=c(rep(1,4), 1.02, 1.02))
                    }
    } else if (type=="lung"){Mortality<-Lifetables$Mortality_obs_all 
                        Rel_surv<-SEER_rel.surv$Rel.surv_lung 
                        Cancer_inc<-SEER_inc_obs$Rate_lung
                        Dem_inc<-AD_inc_obs$Rate_all
                        if (S1.mort.rateratio>0.7){
                          Mort_cal<-c(.95,0.8,0.8,0.8,0.8,0.75,0.75,0.65) 
                          demrate.adj.S_IRR0.3<-list(prev0.5=c(rep(1.5,3),1.65, 1.85, 2.25),
                                                     prev0.4=c(rep(1.45,3),1.5, 1.6, 1.95),
                                                     prev0.3=c(rep(1.3,3),1.35, 1.42, 1.65),
                                                     prev0.2=c(rep(1.15,3),1.18, 1.25, 1.39),
                                                     prev0.1=c(c(rep(1.06,3),1.1, 1.12, 1.17)))
                          demrate.adj.S_IRR0.5<-list(prev0.5=c(rep(1.3,3),1.35, 1.4,1.53),
                                                     prev0.4=c(rep(1.3,3),1.32, 1.34,1.4),
                                                     prev0.3=c(rep(1.2,3),1.22, 1.24,1.29),
                                                     prev0.2=c(rep(1.1,3),1.15, 1.14, 1.19),
                                                     prev0.1=c(c(rep(1.05,4),1.07, 1.09)))
                          demrate.adj.S_IRR0.7<-list(prev0.5=c(rep(1.2,3),1.17,1.2,1.23),
                                                     prev0.4=c(rep(1.2,3),1.16,1.17,1.18),
                                                     prev0.3=c(rep(1.15,3),1.11,1.11,1.12),
                                                     prev0.2=c(rep(1.08,3),1.08, 1.08, 1.09),
                                                     prev0.1=c(c(rep(1.04,5),1.03)))
                          demrate.adj.S_IRR0.9<-list(prev0.5=c(rep(1.05,5),1.05),
                                                     prev0.4=c(rep(1.05,5),1.05),
                                                     prev0.3=c(rep(1.03,5),1.03),
                                                     prev0.2=c(rep(1,3),1.03, 1.03, 1.03),
                                                     prev0.1=rep(1,6))
                        } else if (S1.mort.rateratio>0.5 & S1.mort.rateratio<=0.7){
                            Mort_cal<-c(.95,0.8,0.8,0.8,0.8,0.75,0.75,0.65) 
                            demrate.adj.S_IRR0.3<-list(prev0.5=c(rep(1.5,3),1.65, 1.85, 2.25),
                                                       prev0.4=c(rep(1.4,2),1.45, 1.5, 1.63, 1.97),
                                                       prev0.3=c(rep(1.3,2),1.33, 1.35, 1.42, 1.67),
                                                       prev0.2=c(rep(1.15,3),1.21, 1.25, 1.4),
                                                       prev0.1=c(c(rep(1.06,3),1.08, 1.12, 1.18)))
                            demrate.adj.S_IRR0.5<-list(prev0.5=c(rep(1.3,3),1.38, 1.42,1.53),
                                                       prev0.4=c(rep(1.3,3),1.31, 1.33, 1.42),
                                                       prev0.3=c(rep(1.22,3),1.22, 1.23, 1.3),
                                                       prev0.2=c(rep(1.12,3),1.15, 1.15, 1.2),
                                                       prev0.1=c(c(rep(1.05,4),1.07, 1.1)))
                            demrate.adj.S_IRR0.7<-list(prev0.5=c(rep(1.2,3),1.17,1.2,1.23),
                                                       prev0.4=c(rep(1.16,3),1.16,1.16,1.18),
                                                       prev0.3=c(rep(1.12,3),1.12,1.12,1.14),
                                                       prev0.2=c(rep(1.07,3),1.08, 1.08, 1.1),
                                                       prev0.1=c(c(rep(1.04,5),1.03)))
                            demrate.adj.S_IRR0.9<-list(prev0.5=c(rep(1.05,5),1.05),
                                                       prev0.4=c(rep(1.05,5),1.05),
                                                       prev0.3=c(rep(1.03,3),1.03, 1.03, 1.03),
                                                       prev0.2=c(rep(1,3),1.02, 1.02, 1.02),
                                                       prev0.1=rep(1,6))
                        }else if (S1.mort.rateratio>0.3 & S1.mort.rateratio<=0.5){ 
                          Mort_cal<-c(.9,0.8,0.75,0.75,0.75,0.75,0.75,0.65) 
                            demrate.adj.S_IRR0.3<-list(prev0.5=c(rep(1.52,3), 1.71, 1.85, 2.3),
                                                       prev0.4=c(rep(1.4,3), 1.52, 1.63, 1.98),
                                                       prev0.3=c(rep(1.32,3), 1.36, 1.44, 1.66),
                                                       prev0.2=c(rep(1.18,3),1.23, 1.25, 1.4),
                                                       prev0.1=c(c(rep(1.07,3),1.09, 1.13, 1.19)))
                            demrate.adj.S_IRR0.5<-list(prev0.5=c(rep(1.35,3),1.41, 1.43,1.55),
                                                       prev0.4=c(rep(1.32,3),1.31, 1.34,1.41),
                                                       prev0.3=c(rep(1.23,3),1.22, 1.25,1.3),
                                                       prev0.2=c(rep(1.13,3),1.13, 1.14, 1.2),
                                                       prev0.1=c(c(rep(1.05,3),1.07, 1.06, 1.08)))
                            demrate.adj.S_IRR0.7<-list(prev0.5=c(rep(1.2,3),1.19,1.2,1.21),
                                                       prev0.4=c(rep(1.16,3),1.16,1.17,1.18),
                                                       prev0.3=c(rep(1.12,3),1.13,1.11,1.12),
                                                       prev0.2=c(rep(1.08,3),1.08, 1.08, 1.08),
                                                       prev0.1=c(c(rep(1.05,5),1.06)))
                            demrate.adj.S_IRR0.9<-list(prev0.5=c(rep(1.05,5),1.07),
                                                       prev0.4=c(rep(1.05,5),1.05),
                                                       prev0.3=c(rep(1.03,5),1.04),
                                                       prev0.2=c(rep(1.01,5),1.03),
                                                       prev0.1=c(rep(1,4), 1.01, 1.01))
                          } else{ 
                            Mort_cal<-c(.9,0.85,0.8,0.8,0.75,0.75,0.75,0.65) 
                            demrate.adj.S_IRR0.3<-list(prev0.5=c(rep(1.52,2), 1.55, 1.72, 1.9, 2.35),
                                                       prev0.4=c(rep(1.42,2), 1.42, 1.52, 1.66, 2.01),
                                                       prev0.3=c(rep(1.25,2), 1.28, 1.35, 1.44, 1.68),
                                                       prev0.2=c(rep(1.2,3),1.22, 1.25, 1.43),
                                                       prev0.1=c(c(rep(1.07,3),1.09, 1.1, 1.19)))
                            demrate.adj.S_IRR0.5<-list(prev0.5=c(rep(1.35,3),1.4, 1.45,1.53),
                                                       prev0.4=c(rep(1.3,3),1.31, 1.34,1.43),
                                                       prev0.3=c(rep(1.23,3),1.22, 1.25,1.31),
                                                       prev0.2=c(rep(1.13,3),1.15, 1.16,1.21),
                                                       prev0.1=c(c(rep(1.05,3),1.08, 1.07, 1.1)))
                            demrate.adj.S_IRR0.7<-list(prev0.5=c(rep(1.2,3),1.19,1.22,1.23),
                                                       prev0.4=c(rep(1.17,3),1.17,1.18,1.19),
                                                       prev0.3=c(rep(1.12,3),1.12,1.12,1.14),
                                                       prev0.2=c(rep(1.08,3),1.08,1.08,1.09),
                                                       prev0.1=c(c(rep(1.05,5),1.06)))
                            demrate.adj.S_IRR0.9<-list(prev0.5=c(rep(1.05,5),1.07),
                                                       prev0.4=c(rep(1.05,5),1.05),
                                                       prev0.3=c(rep(1.03,5),1.03),
                                                       prev0.2=c(rep(1.02,5),1.02),
                                                       prev0.1=c(rep(1,4), 1.02, 1.02))
                          }
    
                        
    } else if (type=="breast"){Mortality<-Lifetables$Mortality_obs_female 
                          Rel_surv<-SEER_rel.surv$Rel.surv_breast 
                          Cancer_inc<-SEER_inc_obs$Rate_breast
                          Dem_inc<-AD_inc_obs$Rate_F
                          if (S1.mort.rateratio>0.7){
                            Mort_cal<-c(.95,0.9,0.85,0.85,0.8,0.8,0.75,0.65) 
                            demrate.adj.S_IRR0.3<-list(prev0.5=c(rep(1.5,3),1.65, 1.85, 2.25),
                                                       prev0.4=c(rep(1.45,3),1.5, 1.6, 1.95),
                                                       prev0.3=c(rep(1.3,3),1.35, 1.42, 1.65),
                                                       prev0.2=c(rep(1.15,3),1.18, 1.25, 1.39),
                                                       prev0.1=c(c(rep(1.06,3),1.1, 1.12, 1.17)))
                            demrate.adj.S_IRR0.5<-list(prev0.5=c(rep(1.3,3),1.35, 1.4,1.53),
                                                       prev0.4=c(rep(1.3,3),1.32, 1.34,1.4),
                                                       prev0.3=c(rep(1.2,3),1.22, 1.24,1.29),
                                                       prev0.2=c(rep(1.1,3),1.15, 1.14, 1.19),
                                                       prev0.1=c(c(rep(1.05,4),1.07, 1.09)))
                            demrate.adj.S_IRR0.7<-list(prev0.5=c(rep(1.2,3),1.17,1.2,1.23),
                                                       prev0.4=c(rep(1.2,3),1.16,1.17,1.18),
                                                       prev0.3=c(rep(1.15,3),1.11,1.11,1.12),
                                                       prev0.2=c(rep(1.08,3),1.08, 1.08, 1.09),
                                                       prev0.1=c(c(rep(1.04,5),1.03)))
                            demrate.adj.S_IRR0.9<-list(prev0.5=c(rep(1.05,5),1.05),
                                                       prev0.4=c(rep(1.05,5),1.05),
                                                       prev0.3=c(rep(1.03,5),1.03),
                                                       prev0.2=c(rep(1,3),1.03, 1.03, 1.03),
                                                       prev0.1=rep(1,6))
                          } else if  (S1.mort.rateratio>0.5 & S1.mort.rateratio<=0.7){ 
                            Mort_cal<-c(.95,0.9,0.85,0.8,0.8,0.8,0.75,0.65) 
                            demrate.adj.S_IRR0.3<-list(prev0.5=c(rep(1.5,3),1.65, 1.85, 2.3),
                                                       prev0.4=c(rep(1.4,2),1.45, 1.5, 1.63, 1.99),
                                                       prev0.3=c(rep(1.3,2),1.33, 1.35, 1.42, 1.69),
                                                       prev0.2=c(rep(1.15,3),1.21, 1.25, 1.4),
                                                       prev0.1=c(c(rep(1.06,3),1.08, 1.12, 1.18)))
                            demrate.adj.S_IRR0.5<-list(prev0.5=c(rep(1.3,3),1.38, 1.42,1.53),
                                                       prev0.4=c(rep(1.3,3),1.31, 1.33, 1.42),
                                                       prev0.3=c(rep(1.22,3),1.22, 1.23, 1.3),
                                                       prev0.2=c(rep(1.12,3),1.15, 1.15, 1.2),
                                                       prev0.1=c(c(rep(1.05,4),1.07, 1.1)))
                            demrate.adj.S_IRR0.7<-list(prev0.5=c(rep(1.2,3),1.17,1.2,1.23),
                                                       prev0.4=c(rep(1.16,3),1.16,1.16,1.18),
                                                       prev0.3=c(rep(1.12,3),1.12,1.12,1.14),
                                                       prev0.2=c(rep(1.07,3),1.08, 1.08, 1.1),
                                                       prev0.1=c(c(rep(1.04,5),1.03)))
                            demrate.adj.S_IRR0.9<-list(prev0.5=c(rep(1.05,5),1.05),
                                                       prev0.4=c(rep(1.05,5),1.05),
                                                       prev0.3=c(rep(1.03,3),1.03, 1.03, 1.03),
                                                       prev0.2=c(rep(1,3),1.02, 1.02, 1.02),
                                                       prev0.1=rep(1,6))
                          } else if  (S1.mort.rateratio>0.3 & S1.mort.rateratio<=0.5){ 
                            Mort_cal<-c(.95,0.9,0.85,0.8,0.8,0.8,0.75,0.65) 
                            demrate.adj.S_IRR0.3<-list(prev0.5=c(rep(1.52,3), 1.71, 1.88, 2.38),
                                                       prev0.4=c(rep(1.4,3), 1.52, 1.63, 2.05),
                                                       prev0.3=c(rep(1.32,3), 1.36, 1.44, 1.74),
                                                       prev0.2=c(rep(1.18,3),1.23, 1.27, 1.45),
                                                       prev0.1=c(c(rep(1.07,3),1.09, 1.13, 1.19)))
                            demrate.adj.S_IRR0.5<-list(prev0.5=c(rep(1.35,3),1.41, 1.43,1.55),
                                                       prev0.4=c(rep(1.32,3),1.31, 1.34,1.45),
                                                       prev0.3=c(rep(1.23,3),1.22, 1.25,1.34),
                                                       prev0.2=c(rep(1.13,3),1.13, 1.15, 1.23),
                                                       prev0.1=c(c(rep(1.05,3),1.07, 1.08, 1.12)))
                            demrate.adj.S_IRR0.7<-list(prev0.5=c(rep(1.2,3),1.19,1.2,1.24),
                                                       prev0.4=c(rep(1.16,3),1.16,1.17,1.18),
                                                       prev0.3=c(rep(1.12,3),1.13,1.11,1.15),
                                                       prev0.2=c(rep(1.08,3),1.08, 1.08, 1.1),
                                                       prev0.1=c(c(rep(1.05,5),1.06)))
                            demrate.adj.S_IRR0.9<-list(prev0.5=c(rep(1.05,5),1.07),
                                                       prev0.4=c(rep(1.05,5),1.05),
                                                       prev0.3=c(rep(1.03,5),1.04),
                                                       prev0.2=c(rep(1.01,5),1.03),
                                                       prev0.1=c(rep(1,4), 1.01, 1.01))
                          } else{ 
                            Mort_cal<-c(.95,0.9,0.85,0.8,0.8,0.8,0.75,0.65) 
                            demrate.adj.S_IRR0.3<-list(prev0.5=c(rep(1.52,2), 1.55, 1.72, 1.9, 2.48),
                                                       prev0.4=c(rep(1.42,2), 1.42, 1.52, 1.66, 2.16),
                                                       prev0.3=c(rep(1.25,2), 1.28, 1.35, 1.46, 1.87),
                                                       prev0.2=c(rep(1.2,3),1.22, 1.27, 1.53),
                                                       prev0.1=c(c(rep(1.07,3),1.09, 1.13, 1.25)))
                            demrate.adj.S_IRR0.5<-list(prev0.5=c(rep(1.35,3),1.4, 1.45, 1.62),
                                                       prev0.4=c(rep(1.3,3),1.31, 1.34, 1.49),
                                                       prev0.3=c(rep(1.23,3),1.22, 1.25,1.38),
                                                       prev0.2=c(rep(1.13,3),1.15, 1.16,1.26),
                                                       prev0.1=c(c(rep(1.05,3),1.08, 1.07, 1.13)))
                            demrate.adj.S_IRR0.7<-list(prev0.5=c(rep(1.2,3),1.19,1.22,1.27),
                                                       prev0.4=c(rep(1.17,3),1.17,1.18,1.21),
                                                       prev0.3=c(rep(1.12,3),1.12,1.12,1.16),
                                                       prev0.2=c(rep(1.08,3),1.08,1.08,1.12),
                                                       prev0.1=c(c(rep(1.05,5),1.06)))
                            demrate.adj.S_IRR0.9<-list(prev0.5=c(rep(1.05,5),1.07),
                                                       prev0.4=c(rep(1.05,5),1.05),
                                                       prev0.3=c(rep(1.03,5),1.03),
                                                       prev0.2=c(rep(1.02,5),1.02),
                                                       prev0.1=c(rep(1,4), 1.02, 1.02))
                          }
  
    } else if (type=="prostate"){Mortality<-Lifetables$Mortality_obs_male 
                            Rel_surv<-SEER_rel.surv$Rel.surv_prostate 
                            Cancer_inc<-SEER_inc_obs$Rate_prostate
                            Dem_inc<-AD_inc_obs$Rate_M
                            if (S1.mort.rateratio>0.7){
                              Mort_cal<-c(.95,0.9,0.85,0.85,0.8,0.8,0.75,0.65) 
                              demrate.adj.S_IRR0.3<-list(prev0.5=c(rep(1.5,3),1.65, 1.85, 2.25),
                                                         prev0.4=c(rep(1.45,3),1.5, 1.6, 1.95),
                                                         prev0.3=c(rep(1.3,3),1.35, 1.42, 1.65),
                                                         prev0.2=c(rep(1.15,3),1.18, 1.25, 1.39),
                                                         prev0.1=c(c(rep(1.06,3),1.1, 1.12, 1.17)))
                              demrate.adj.S_IRR0.5<-list(prev0.5=c(rep(1.3,3),1.35, 1.4,1.53),
                                                         prev0.4=c(rep(1.3,3),1.32, 1.34,1.4),
                                                         prev0.3=c(rep(1.2,3),1.22, 1.24,1.29),
                                                         prev0.2=c(rep(1.1,3),1.15, 1.14, 1.19),
                                                         prev0.1=c(c(rep(1.05,4),1.07, 1.09)))
                              demrate.adj.S_IRR0.7<-list(prev0.5=c(rep(1.2,3),1.17,1.2,1.23),
                                                         prev0.4=c(rep(1.2,3),1.16,1.17,1.18),
                                                         prev0.3=c(rep(1.15,3),1.11,1.11,1.12),
                                                         prev0.2=c(rep(1.08,3),1.08, 1.08, 1.09),
                                                         prev0.1=c(c(rep(1.04,5),1.03)))
                              demrate.adj.S_IRR0.9<-list(prev0.5=c(rep(1.05,5),1.05),
                                                         prev0.4=c(rep(1.05,5),1.05),
                                                         prev0.3=c(rep(1.03,5),1.03),
                                                         prev0.2=c(rep(1,3),1.03, 1.03, 1.03),
                                                         prev0.1=rep(1,6))
                            } else if (S1.mort.rateratio>0.5 & S1.mort.rateratio<=0.7){
                              Mort_cal<-c(.95,0.9,0.85,0.85,0.8,0.8,0.75,0.65) 
                              demrate.adj.S_IRR0.3<-list(prev0.5=c(rep(1.5,3),1.69, 1.88, 2.3),
                                                         prev0.4=c(rep(1.4,2),1.45, 1.5, 1.63, 1.99),
                                                         prev0.3=c(rep(1.3,2),1.33, 1.35, 1.42, 1.69),
                                                         prev0.2=c(rep(1.15,3),1.21, 1.25, 1.4),
                                                         prev0.1=c(c(rep(1.06,3),1.08, 1.12, 1.18)))
                              demrate.adj.S_IRR0.5<-list(prev0.5=c(rep(1.3,3),1.38, 1.42,1.53),
                                                         prev0.4=c(rep(1.3,3),1.31, 1.33, 1.42),
                                                         prev0.3=c(rep(1.22,3),1.22, 1.23, 1.32),
                                                         prev0.2=c(rep(1.12,3),1.15, 1.15, 1.2),
                                                         prev0.1=c(c(rep(1.05,4),1.07, 1.1)))
                              demrate.adj.S_IRR0.7<-list(prev0.5=c(rep(1.2,3),1.17,1.2,1.23),
                                                         prev0.4=c(rep(1.16,3),1.16,1.16,1.18),
                                                         prev0.3=c(rep(1.12,3),1.12,1.12,1.14),
                                                         prev0.2=c(rep(1.07,3),1.08, 1.08, 1.1),
                                                         prev0.1=c(c(rep(1.04,5),1.03)))
                              demrate.adj.S_IRR0.9<-list(prev0.5=c(rep(1.05,5),1.05),
                                                         prev0.4=c(rep(1.05,5),1.05),
                                                         prev0.3=c(rep(1.03,3),1.03, 1.03, 1.03),
                                                         prev0.2=c(rep(1,3),1.02, 1.02, 1.02),
                                                         prev0.1=c(rep(1,5),1.02))
                            } else if (S1.mort.rateratio>0.3 & S1.mort.rateratio<=0.5) { 
                              Mort_cal<-c(.9,0.85,0.8,0.8,0.75,0.75,0.75,0.65) 
                              demrate.adj.S_IRR0.3<-list(prev0.5=c(rep(1.52,3), 1.71, 1.88, 2.38),
                                                         prev0.4=c(rep(1.4,3), 1.52, 1.65, 2.07),
                                                         prev0.3=c(rep(1.32,3), 1.36, 1.46, 1.77),
                                                         prev0.2=c(rep(1.18,3),1.23, 1.27, 1.47),
                                                         prev0.1=c(c(rep(1.07,3),1.09, 1.13, 1.22)))
                              demrate.adj.S_IRR0.5<-list(prev0.5=c(rep(1.35,3),1.41, 1.43,1.58),
                                                         prev0.4=c(rep(1.32,3),1.31, 1.34,1.45),
                                                         prev0.3=c(rep(1.23,3),1.22, 1.25,1.34),
                                                         prev0.2=c(rep(1.13,3),1.13, 1.15, 1.23),
                                                         prev0.1=c(c(rep(1.05,3),1.07, 1.08, 1.12)))
                              demrate.adj.S_IRR0.7<-list(prev0.5=c(rep(1.2,3),1.19,1.2,1.24),
                                                         prev0.4=c(rep(1.16,3),1.16,1.17,1.18),
                                                         prev0.3=c(rep(1.12,3),1.13,1.11,1.15),
                                                         prev0.2=c(rep(1.08,3),1.08, 1.08, 1.1),
                                                         prev0.1=c(c(rep(1.05,5),1.06)))
                              demrate.adj.S_IRR0.9<-list(prev0.5=c(rep(1.05,5),1.07),
                                                         prev0.4=c(rep(1.05,5),1.05),
                                                         prev0.3=c(rep(1.03,5),1.04),
                                                         prev0.2=c(rep(1.01,5),1.03),
                                                         prev0.1=c(rep(1,4), 1.01, 1.01))
                            } else{ 
                              Mort_cal<-c(.9,0.9,0.85,0.85,0.75,0.75,0.75,0.65) 
                              demrate.adj.S_IRR0.3<-list(prev0.5=c(rep(1.52,2), 1.55, 1.72, 1.97, 2.55),
                                                         prev0.4=c(rep(1.42,2), 1.42, 1.52, 1.71, 2.24),
                                                         prev0.3=c(rep(1.25,2), 1.29, 1.37, 1.49, 1.91),
                                                         prev0.2=c(rep(1.2,3),1.22, 1.3, 1.6),
                                                         prev0.1=c(c(rep(1.07,3),1.09, 1.13, 1.27)))
                              demrate.adj.S_IRR0.5<-list(prev0.5=c(rep(1.35,3),1.4, 1.45,1.62),
                                                         prev0.4=c(rep(1.3,3),1.31, 1.36, 1.53),
                                                         prev0.3=c(rep(1.23,3),1.23, 1.28,1.42),
                                                         prev0.2=c(rep(1.13,3),1.15, 1.16,1.28),
                                                         prev0.1=c(c(rep(1.05,3),1.08, 1.07, 1.13)))
                              demrate.adj.S_IRR0.7<-list(prev0.5=c(rep(1.2,3),1.19,1.22,1.27),
                                                         prev0.4=c(rep(1.17,3),1.17,1.18,1.21),
                                                         prev0.3=c(rep(1.12,3),1.12,1.13,1.19),
                                                         prev0.2=c(rep(1.08,3),1.08,1.08,1.12),
                                                         prev0.1=c(c(rep(1.05,5),1.06)))
                              demrate.adj.S_IRR0.9<-list(prev0.5=c(rep(1.05,5),1.07),
                                                         prev0.4=c(rep(1.05,5),1.05),
                                                         prev0.3=c(rep(1.03,5),1.05),
                                                         prev0.2=c(rep(1.02,5),1.02),
                                                         prev0.1=c(rep(1,4), 1.02, 1.02))
                            }
    }

  
#------------------------------------------------------------------
# Part 2: Defining input rates based on time (i.e. age)
  # In this section, I define several functions for rates corresponding to 
  # different states (e.g. C0D1S0). I start with a base rate and apply 
  # rate ratios to get rates for other groups.   
#------------------------------------------------------------------
 
    #------------------------------------------------------------------
    # Defining mortality rates #
      # There are 8 mortality rates for the 8 states 
      # States are combos of Cancer, Dementia, and Survival characteristic
    #------------------------------------------------------------------
  
    # Define reference mortality rate for non-cancer, non-dementia as a 
      # function of time using lifetables data for whites
        mortality.rate.C0D0 <- function(time){
          if (immortal==1){0}
          else if (time<5){Mortality[2]*Mort_cal[1]}
          else if (5<=time & time<10){Mortality[3]*Mort_cal[2]}
          else if (10<=time & time<15){Mortality[4]*Mort_cal[3]}
          else if (15<=time & time<20){Mortality[5]*Mort_cal[4]}
          else if (20<=time & time<25){Mortality[6]*Mort_cal[5]}
          else if (25<=time & time<30){Mortality[7]*Mort_cal[6]}
          else if (30<=time & time<35){Mortality[8]*Mort_cal[7]}
          else if (35<=time){Mortality[9]*Mort_cal[8]}
        }
        
        
    # Use SEER data to determine relative mortality of cancer patients
    # The relative survival is probability of survival relative to non-cancer
      # probability of survival (taken from SEER)
        
        #This function defines mortality IRR for cancer using relative survival
        mort.cancer.irr<-function(time,rel.surv.5yrs){
          mort.cancer.rate = mortality.rate.C0D0(time) - 1/5*log(rel.surv.5yrs)
          mort.cancer.rate/mortality.rate.C0D0(time)
        }

        #This function defines mortality rates, applying IRR to base rate
        mortality.rate.C1D0 <- function(time){
          if (immortal==1){0
          } else if (diff.cancer.mort==0){mortality.rate.C0D0(time)
          } else if (diff.cancer.mort==1) {
              if (time<10) {mortality.rate.C0D0(time)*mort.cancer.irr(time,Rel_surv[1])
              } else if (time>=10) {mortality.rate.C0D0(time)*mort.cancer.irr(time,Rel_surv[2])} 
          }
        }
    
    
    # Use Mayeda et al. 2017 paper to assign relative mortality of dementia patients
      # The mortality rate is HR for mortality among dementa vs. no dementia (whites)
          # times mortality rate calculated from lifetables 
        mortality.rate.C0D1 <- function(time){
          if (immortal==1){0}
          else if (diff.dem.mort==0){mortality.rate.C0D0(time)}
          else if (diff.dem.mort==1) {
            if (time<5) {mortality.rate.C0D0(time)*AD_HR_obs$Dementia.HR[1]} 
            else if (time>=5 & time<10) {mortality.rate.C0D0(time)*AD_HR_obs$Dementia.HR[2]}
            else if (time>=10 & time<15) {mortality.rate.C0D0(time)*AD_HR_obs$Dementia.HR[3]}
            else if (time>=15 & time<20) {mortality.rate.C0D0(time)*AD_HR_obs$Dementia.HR[4]}
            else if (time>=20 & time<25) {mortality.rate.C0D0(time)*AD_HR_obs$Dementia.HR[5]}
            else if (time>=25) {mortality.rate.C0D0(time)*AD_HR_obs$Dementia.HR[6]}
          }
          }
  
    #Combine SEER and Mayeda data to assign mortality rate for 
      # cancer AND dementia patients, assuming multiplicative relationship.
        mortality.rate.C1D1 <- function(time){
          if (immortal==1){0}
          else if (diff.cancer.mort==0){mortality.rate.C0D1(time)}
          else if (diff.dem.mort==0){mortality.rate.C1D0(time)}
          else if (diff.cancer.mort==1 & diff.dem.mort==1) {
            if (time<10) {mortality.rate.C0D1(time)*mort.cancer.irr(time,Rel_surv[1])} 
            else if (time>=10) {mortality.rate.C0D1(time)*mort.cancer.irr(time,Rel_surv[2])} 
            }
        }
    
    
    #Define final mortality rates for all 8 states          
    
        # For cancer-free individuals, mortality rate is rate for their dementia 
        #   status, regardless of S
            mortality.rate.C0D0S0<-function(time){mortality.rate.C0D0(time)}
            mortality.rate.C0D0S1<-function(time){mortality.rate.C0D0(time)}
            mortality.rate.C0D1S0<-function(time){mortality.rate.C0D1(time)}
            mortality.rate.C0D1S1<-function(time){mortality.rate.C0D1(time)}
            
        # For those with cancer and no S, mortality rate is rate for their 
        #   cancer/dementia status 
            mortality.rate.C1D1S0<-function(time){mortality.rate.C1D1(time)}
            mortality.rate.C1D0S0<-function(time){mortality.rate.C1D0(time)}
            
        # For those with cancer and S=1, mortality rate is rate for their 
        #   cancer/dementia status  * protective effect of S
            mortality.rate.C1D1S1<-function(time){mortality.rate.C1D1(time)*S1.mort.rateratio}
            mortality.rate.C1D0S1<-function(time){mortality.rate.C1D0(time)*S1.mort.rateratio}
            
        
    #------------------------------------------------------------------
    # Defining cancer incidence rates based on time (i.e. age)
    #------------------------------------------------------------------
            
            #Define cancer incidence rate as a function of time using SEER data
            get.cancer <- function(time){
              if (time<5){Cancer_inc[1]}
              else if (5<=time & time<10){Cancer_inc[2]}
              else if (10<=time & time<15){Cancer_inc[3]}
              else if (15<=time & time<20){Cancer_inc[4]}
              else if (20<=time){Cancer_inc[5]}
            }
            
    #------------------------------------------------------------------
    # Defining dementia incidence rates based on time (i.e. age)
    #------------------------------------------------------------------

        Prevlist<-c(0.5, 0.4, 0.3, 0.2, 0.1) #List in the order specified above for the list of demrate_adj
        demratecorrection<- 
          if (S1.dem.rateratio==1 | Sprev==0) {rep(1,6)
          } else if (S1.dem.rateratio==0.3) {demrate.adj.S_IRR0.3[[which(Prevlist==Sprev)]]
          } else if (S1.dem.rateratio==0.5) {demrate.adj.S_IRR0.5[[which(Prevlist==Sprev)]]
          } else if (S1.dem.rateratio==0.7) {demrate.adj.S_IRR0.7[[which(Prevlist==Sprev)]]
          } else if (S1.dem.rateratio==0.9) {demrate.adj.S_IRR0.9[[which(Prevlist==Sprev)]]}

        
        #Define reference dementia incidence rate as a function of time using ACT data
          get.dem <- function(time){
              if (time<5){Dem_inc[2]*demratecorrection[1]}
              else if (5<=time & time<10){Dem_inc[3]*demratecorrection[2]}
              else if (10<=time & time<15){Dem_inc[4]*demratecorrection[3]}
              else if (15<=time & time<20){Dem_inc[5]*demratecorrection[4]}
              else if (20<=time & time<25){Dem_inc[6]*demratecorrection[5]}
              else if (25<=time){Dem_inc[7]*demratecorrection[6]}
            }  
        
        #For those with survival characteristic, reference rate * protective IRR
          get.dem.S1<-function(time){get.dem(time)*S1.dem.rateratio}
        
        #For those without survival characteristic, dementia incidence rate is reference rate 
          get.dem.S0<-function(time){get.dem(time)}  
          

#------------------------------------------------------------------
# Part 3: Specify model
#      #Define the differential equations for the model#
#      #Define initial conditions for the model#
#------------------------------------------------------------------
# diff.eq.model() is a function which defines the set of differential equations
  # inputs: time, states, parameters
      ## time: a vector of times at which you want the derivates
      ## states: a vector of number of people in each state with the entries labeled
      ## parameters: a vector of the time-constant parameters (rates) of the model. 
            #This vector is empty for this model, but the model won't run without it. 
  # outputs a list where the first entry is a vector derivatives and the second entry is NULL
      ## this output format is required by the differential equation solver
  
  diff.eq.model<- function(time,states, parameters){
    with(as.list(c(states, parameters)),
         {
           #Rates for 9 states in the model
           dC0D0S0 <- -get.cancer(time)*C0D0S0 - get.dem.S0(time)*C0D0S0 - 
                          mortality.rate.C0D0S0(time)*C0D0S0
           dC0D1S0 <- -get.cancer(time)*C0D1S0 + get.dem.S0(time)*C0D0S0 - 
                          mortality.rate.C0D1S0(time)*C0D1S0
           dC1D0S0 <-  get.cancer(time)*C0D0S0 - get.dem.S0(time)*C1D0S0 - 
                          mortality.rate.C1D0S0(time)*C1D0S0
           dC1D1S0 <-  get.cancer(time)*C0D1S0 + get.dem.S0(time)*C1D0S0 - 
                          mortality.rate.C1D1S0(time)*C1D1S0
           dC0D0S1 <- -get.cancer(time)*C0D0S1 - get.dem.S1(time)*C0D0S1 - 
                          mortality.rate.C0D0S1(time)*C0D0S1
           dC0D1S1 <- -get.cancer(time)*C0D1S1 + get.dem.S1(time)*C0D0S1 - 
                          mortality.rate.C0D1S1(time)*C0D1S1
           dC1D0S1 <-  get.cancer(time)*C0D0S1 - get.dem.S1(time)*C1D0S1 - 
                          mortality.rate.C1D0S1(time)*C1D0S1
           dC1D1S1 <-  get.cancer(time)*C0D1S1 + get.dem.S1(time)*C1D0S1 - 
                          mortality.rate.C1D1S1(time)*C1D1S1
           dDEAD <-  mortality.rate.C0D0S0(time)*C0D0S0 + 
                      mortality.rate.C1D0S0(time)*C1D0S0 + 
                      mortality.rate.C0D1S0(time)*C0D1S0 + 
                      mortality.rate.C1D1S0(time)*C1D1S0 + 
                      mortality.rate.C0D0S1(time)*C0D0S1 + 
                      mortality.rate.C1D0S1(time)*C1D0S1 + 
                      mortality.rate.C0D1S1(time)*C0D1S1 + 
                      mortality.rate.C1D1S1(time)*C1D1S1
           
           
           #The rest of these are for tallying results (e.g. cumul. incidence)
           dDEMNOC <- get.dem.S0(time)*C0D0S0 + get.dem.S1(time)*C0D0S1 
           dDEMWC <- get.dem.S0(time)*C1D0S0 + get.dem.S1(time)*C1D0S1
           dDEMENTIA <- get.dem.S0(time)*(C0D0S0 + C1D0S0) +  get.dem.S1(time)*(C0D0S1 + C1D0S1)
           dCANCER <- get.cancer(time)*(C0D0S0 + C0D1S0 + C0D0S1 + C0D1S1)
           dDEADC0D0S0 <- mortality.rate.C0D0S0(time)*C0D0S0
           dDEADC1D0S0 <- mortality.rate.C1D0S0(time)*C1D0S0
           dDEADC0D1S0 <- mortality.rate.C0D1S0(time)*C0D1S0
           dDEADC1D1S0 <- mortality.rate.C1D1S0(time)*C1D1S0
           dDEADC0D0S1 <- mortality.rate.C0D0S1(time)*C0D0S1
           dDEADC1D0S1 <- mortality.rate.C1D0S1(time)*C1D0S1
           dDEADC0D1S1 <- mortality.rate.C0D1S1(time)*C0D1S1
           dDEADC1D1S1 <- mortality.rate.C1D1S1(time)*C1D1S1
           dCBEFOREDEM <- get.cancer(time)*(C0D0S0+C0D0S1)
           dDEMNOSNOC <- get.dem.S0(time)*(C0D0S0) 
           dDEMWITHSNOC <-  get.dem.S1(time)*(C0D0S1)
           
           list(c(dC0D0S0,dC0D1S0,dC1D0S0,dC1D1S0,dC0D0S1,dC0D1S1,dC1D0S1,dC1D1S1,dDEAD,
                  dDEMNOC,dDEMWC,dDEMENTIA, dCANCER, dDEADC0D0S0, dDEADC1D0S0, 
                  dDEADC0D1S0, dDEADC1D1S0, dDEADC0D0S1, dDEADC1D0S1, dDEADC0D1S1, dDEADC1D1S1, dCBEFOREDEM, dDEMNOSNOC, dDEMWITHSNOC),NULL)
         })}


# gen.init() is a function which generates initial conditions of the model. 
gen.init <- function(parameters){
  with(as.list(parameters),
       {
         initial.conditions <- c(C0D0S0=1*(1-Sprev),C0D1S0=0,C1D0S0=0,C1D1S0=0,
                                 C0D0S1=Sprev,C0D1S1=0,C1D0S1=0,C1D1S1=0, DEAD=0, 
                                 DEMNOC=0,DEMWC=0, DEMENTIA=0, CANCER=0, 
                                 DEADC0D0S0=0, DEADC1D0S0=0, DEADC0D1S0=0, DEADC1D1S0=0,  
                                 DEADC0D0S1=0, DEADC1D0S1=0, DEADC0D1S1=0, DEADC1D1S1=0, 
                                 CBEFOREDEM=0, DEMNOSNOC=0, DEMWITHSNOC=0)	
         initial.conditions 
       })}



#------------------------------------------------------------------
# Part 4: Run simulations
#------------------------------------------------------------------
  # run.sim() is a function which numerically solves the differential equation
  # inputs: parameters of the model and times where you want to know the solution
  # lsoda automatically determines the solving method and the time step to get 
      # the solution to within some tolerance at the specified times
    run.sim <- function(parameters,times){
        lsoda(gen.init(parameters),times,diff.eq.model,parameters)
      }
  
  #RUN SIMULATION
      #THis runs the simulation and generates teh state distribution dataframe
      out.data <- run.sim(pars,times)
      out.data %<>% as.data.frame
      
      #THis saves the rate results
      rate.results<-matrix(NA, nrow=length(times), ncol=24)
  
      for (i in 1:length(times)){
              rate.results[i,1:24]<- 
                t(diff.eq.model(time=round(times[i],2),
                                out.data[round(out.data$time,2)==round(times[i],2),],
                                pars))[[1]]
            }
      
      rate.results<-cbind(times,rate.results)
      colnames(rate.results)<-c("time", 
                                "dC0D0S0","dC0D1S0","dC1D0S0","dC1D1S0",
                                "dC0D0S1","dC0D1S1","dC1D0S1","dC1D1S1", "dDEAD",
                                "dDEMNOC","dDEMWC","dDEMENTIA", "dCANCER", 
                                "dDEADC0D0S0", "dDEADC1D0S0", "dDEADC0D1S0", 
                                "dDEADC1D1S0", "dDEADC0D0S1", "dDEADC1D0S1", 
                                "dDEADC0D1S1", "dDEADC1D1S1", 
                               "dCBEFOREDEM", "dDEMNOSNOC", "dDEMWITHSNOC") #This order has to match as above. 

  #Outputs for full compartmodel function 
      # out.data is a matrix with each row for the # of people in each state at input times
      # rate.results is a matrix of rates for each states at each input time.
      
    return(list(out.data=out.data, rate.results=rate.results))
}


#------------------------------------------------------------------
# Part 5. Check to make sure model runs! 
# Comment out this section when using this script as a source file
#------------------------------------------------------------------
# library(here)
# source(here("Code","Multistate model", "Rate inputs.R"))
# testmodel<-compartmodel(0,1,1,0.5,0.5,0.5,c(seq(from=0, to=40, by=1/12),10000),"all",c())
# head(testmodel$out.data)
# head(testmodel$rate.results)
