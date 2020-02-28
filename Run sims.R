#------------------------------------------------------------------
# Title: Run simulations
# Author: Eleanor Hayes-Larson
#------------------------------------------------------------------

library(ggplot2)
library(here)

#------------------------------------------------------------------
# Load model 
#------------------------------------------------------------------
source(here("Code", "Multistate model","Rate inputs.R"))
source(here("Code", "Multistate model","Model.R"))

#------------------------------------------------------------------
# Create function that runs model and produces results. 
#------------------------------------------------------------------
multirun<-function(arg_S1.dem.rateratio, arg_S1.mort.rateratio, arg_Sprev, cancertype){
  modresults<-list()

#Run model
  modresults<-compartmodel(immortal=0, 
                   diff.cancer.mort = 1, 
                   diff.dem.mort = 1, 
                   S1.dem.rateratio = arg_S1.dem.rateratio, 
                   S1.mort.rateratio = arg_S1.mort.rateratio, 
                   Sprev = arg_Sprev,
                   times=c(seq(from=0, to=40, by=1/12),10000),
                   type= cancertype,
                   pars=c())

#Get model results

    #Format results
    out.data<-modresults$out.data
    rate.results<-data.frame(modresults$rate.results)
    
    #Get lifetime relative risk
    relativerisk<-data.frame(time=out.data$time, RR=((out.data$DEMWC/out.data$CBEFOREDEM)/(out.data$DEMNOC/(1-out.data$CBEFOREDEM))))
    lifetimeRR<-relativerisk$RR[relativerisk$time==10000] #request RR for very large time to ensure steady state--i.e. all have died.
    lifetimeRR
    
    #Get time-specific rate ratio from diffEQ model output
    IRR.months<-rate.results$dDEMWC/rate.results$dDEMNOC*((out.data$C0D0S0+out.data$C0D0S1)/(out.data$C1D0S0+out.data$C1D0S1))
    weights<-out.data$C0D0S0+out.data$C1D0S0+out.data$C0D0S1+out.data$C1D0S1
    
    overallIRR<-exp(weighted.mean(log(IRR.months[2:481]),weights[2:481]))
    overallIRR
    
    decade.IRRs<-rep(NA,4)
    for (i in 1:4){
      if (i==1){decade.IRRs[i]<-exp(weighted.mean(log(IRR.months[2:(i*120)]),weights[2:(i*120)]))}
      else {decade.IRRs[i]<-exp(weighted.mean(log(IRR.months[((i-1)*120+1):(i*120)]),weights[((i-1)*120+1):(i*120)]))}
    }
    decade.IRRs
    
#Return model parameter values and results
    return(list(cancertype=paste(cancertype), S1.dem.rateratio=arg_S1.dem.rateratio, S1.mort.rateratio=arg_S1.mort.rateratio, Sprev=arg_Sprev, 
            lifetimeRR=lifetimeRR,
            IRR.overall=overallIRR,
            IRR.6575=decade.IRRs[1],
            IRR.7585=decade.IRRs[2], 
            IRR.8595=decade.IRRs[3], 
            IRR.95105=decade.IRRs[4],
            outdata=out.data,
            rate.results=rate.results))
}


#------------------------------------------------------------------
# Define parameter values of interest 
#------------------------------------------------------------------
S1demvals<-c(1, 0.9, 0.7, 0.5, 0.3)
S1mortvals<-c(1, 0.9, 0.7, 0.5, 0.3)
Sprevvals<-c(0.5, 0.4, 0.3, 0.2, 0.1, 0)
typevals<-c("all", "lung", "breast", "prostate")

args<-expand.grid(S1demvals=S1demvals, S1mortvals=S1mortvals, Sprevvals=Sprevvals, typevals=typevals)


#------------------------------------------------------------------
# Run model across parameter values and store results 
#------------------------------------------------------------------
Allsims_CRSB<-mapply(FUN = multirun, arg_S1.dem.rateratio=args$S1demvals, arg_S1.mort.rateratio=args$S1mortvals, arg_Sprev=args$Sprevvals, cancertype=args$typevals)
                
save(Allsims_CRSB, file=here("Output", "Allsims_CRSB.Rdata"))





