##############################################################################
##### Simulation of functional additive and non-additive genetic effects #####
####### using statistical estimates from quantitative genetic models #########
############### by Thinh Tuan Chu, QGG, Aarhus University ####################
## Note that not all steps and inputs are used in a scenario. However, they ##
#are all needed due to laziness & minimum changes needed for diff. scenarios.#
##############################################################################

##############################################################################
## Example 3: Plot phenotypes from 3-way hybrid crop breeding using SS_NOIA ##
##############################################################################

rm(list = ls()); cat("\014")   #Clear lists and console 
# starting log file & remove old files from previous runs
outlog <- "out.log"
temp1 <- c(outlog,'tempSol1.txt','tempSol2.txt','temp.txt','varSchemTrue.res','varDMUparout1.res','varDMUparout2.res','varDMUparout3.res','varDMUparout4.res','varDMUparout5.res',
           'varFounder.res','varSchemTrueAvCPlotOne.res')
for (itemp in temp1) { 
  if (file.exists(itemp)) file.remove(itemp) } 
write('***************************SimPops************************',file=outlog,append=file.exists(outlog),sep='\t')
write('*This program simulates QTL effects for multi-populations*',file=outlog,append=file.exists(outlog),sep='\t')
write('*allowing genetic differences among pops by changing allel',file=outlog,append=file.exists(outlog),sep='\t')
write('*frequency of pop2-popn. Program written by T.T. Chu, QGG*',file=outlog,append=file.exists(outlog),sep='\t')
write('',file=outlog,append=file.exists(outlog),sep='\t')
write('',file=outlog,append=file.exists(outlog),sep='\t')

# read inputs from s0input1.R
filename <- "S0input1.R"
if (!file.exists(filename)) {
  write(paste('Error: The file',filename,'does not exists. Please make sure all inputs are given to',filename),file=outlog,append=file.exists(outlog))}
write(paste('Reading',filename,'. Please make sure all inputs given in',filename,'are correct.'),file=outlog,append=file.exists(outlog))
source(filename)

####### Load required modules for generating founder pops #########
# Read founder haplotype and genome structures
filename <- paste0(programdir,"S0readFounderGenome1.R")
if (!file.exists(filename)) {  write(paste('Error: The file',filename,'does not exists. This module is needed for',filename),file=outlog,append=file.exists(outlog))}
source(filename)

# Prior effects. Not necessary leading to correct results.
filename <- paste0(programdir,"S1CalculateQtlEffects10.R")
if (!file.exists(filename)) {  write(paste('Error: The file',filename,'does not exists. This module is needed for',filename),file=outlog,append=file.exists(outlog))}
source(filename)

# Generate prior phenotype for base populations
filename <- paste0(programdir,"S2BasePop1.R")
if (!file.exists(filename)) {  write(paste('Error: The file',filename,'does not exists. This module is needed for',filename),file=outlog,append=file.exists(outlog))}
source(filename)

# Load required functions for sampling
filename <- paste0(programdir,"S2SampleFuncs7.R")
if (!file.exists(filename)) {  write(paste('Error: The file',filename,'does not exists. This module is needed for',filename),file=outlog,append=file.exists(outlog))}
source(filename)

#save(list=ls(),file='fsave0.rda')
#load('fsave0.rda')

###########   Finding QTL effects to meet inputs of   #############
#### aRes_covar, aNRMS_covar,aaNRMSRes_covar & dNRMSRes_covar #####
###################################################################
#
# Generate genotypes of base population(s) used for finding QTL 
idLst <- pop[pop$fixedFactor1==0 & pop$currentpop %in% c(1),'id']
nAnimalTemp=length(idLst)
# Create hybrid pop for ref pop of 3-way hybrid
for (currenttime in 1:2) {
  for (currentpop in 1:2) {   # currentpop=1; currenttime=1
    idLst1 <- sample(which(pop$currentpop==currentpop & pop$fixedFactor1==currenttime-1),nAnimalTemp,replace=F)
    idLst2 <- sample(which(pop$currentpop==currentpop & pop$fixedFactor1==currenttime-1),nAnimalTemp,replace=F)
    
    idpoptime <- c((maxid+1):(maxid+length(idLst1)))
    # prestep to prepare for mapply 
    pop[(nrow(pop)+1):(nrow(pop)+length(idLst1)),]=0
    haplo_pat <- rbind(haplo_pat,matrix(ncol=ncol(haplo_pat),nrow=length(idLst1)))
    haplo_mat <- rbind(haplo_mat,matrix(ncol=ncol(haplo_mat),nrow=length(idLst1)))
    maxid=maxid+length(idLst1)
    # F1: mapply sample all offsprings at once
    invisible(mapply(sampleOffspring,idd=idpoptime,par1=idLst1,par2=idLst2,fixedFactors=list(c(currenttime,0)),popid=currentpop))
    
    #F2-Fn
    for (i in 1:nselfing) { #i=1
      invisible(mapply(sampleOffspring,idd=idpoptime,par1=idpoptime,par2=idpoptime,fixedFactors=list(c(currenttime,i+1)),popid=currentpop))			# mapply sample all offsprings at once
    }	
  }  
  if (currenttime>1) {
    # make two-way hybrids  
    currentpop=3
    
    idLst1 <- sample(which(pop$currentpop==2 & pop$fixedFactor1==currenttime),nAnimalTemp,replace=F)  		# select NR
    idLst2 <- sample(which(pop$currentpop==2 & pop$fixedFactor1==currenttime-1),nAnimalTemp,replace=F) 		# select MS
    
    idpoptime <- c((maxid+1):(maxid+length(idLst1)))
    # prestep to prepare for mapply 
    pop[(nrow(pop)+1):(nrow(pop)+length(idLst1)),]=0
    haplo_pat <- rbind(haplo_pat,matrix(ncol=ncol(haplo_pat),nrow=length(idLst1)))
    haplo_mat <- rbind(haplo_mat,matrix(ncol=ncol(haplo_mat),nrow=length(idLst1)))
    maxid=maxid+length(idLst1)
    # F1: mapply sample all offsprings at once
    invisible(mapply(sampleOffspring,idd=idpoptime,par1=idLst1,par2=idLst2,fixedFactors=list(c(currenttime,0)),popid=currentpop))
    
    # make two-way hybrids  
    currentpop=4
    idLst1 <- sample(which(pop$currentpop==1 & pop$fixedFactor1==currenttime),nAnimalTemp,replace=F)  		# select NR
    idLst2 <- sample(which(pop$currentpop==3 & pop$fixedFactor1==currenttime),nAnimalTemp,replace=F) 	  	# select MS
    
    idpoptime <- c((maxid+1):(maxid+length(idLst1)))
    # prestep to prepare for mapply 
    pop[(nrow(pop)+1):(nrow(pop)+length(idLst1)),]=0
    haplo_pat <- rbind(haplo_pat,matrix(ncol=ncol(haplo_pat),nrow=length(idLst1)))
    haplo_mat <- rbind(haplo_mat,matrix(ncol=ncol(haplo_mat),nrow=length(idLst1)))
    maxid=maxid+length(idLst1)
    # F1: mapply sample all offsprings at once
    invisible(mapply(sampleOffspring,idd=idpoptime,par1=idLst1,par2=idLst2,fixedFactors=list(c(currenttime,0)),popid=currentpop))
    
  }  
}

idLst <- pop[pop$fixedFactor1==0 & pop$currentpop %in% c(1),'id']
genoFounderRes=haplo_pat[idLst,]+haplo_mat[idLst,]-2

idLst <- pop[pop$fixedFactor1==0 & pop$currentpop %in% c(2),'id']
genoFounderNRMS=haplo_pat[idLst,]+haplo_mat[idLst,]-2

idLst <- pop[pop$fixedFactor1==2 & pop$currentpop %in% c(1),'id']
genoinbredRes=haplo_pat[idLst,]+haplo_mat[idLst,]-2

idLst <- pop[pop$fixedFactor1==2 & pop$currentpop %in% c(2),'id']
genoinbredNRMS=haplo_pat[idLst,]+haplo_mat[idLst,]-2

idLst <- pop[pop$fixedFactor1>0 & pop$currentpop %in% c(3),'id']
geno2wayNRMS=haplo_pat[idLst,]+haplo_mat[idLst,]-2

idLst <- pop[pop$fixedFactor1>0 & pop$currentpop %in% c(4),'id']
geno3way=haplo_pat[idLst,]+haplo_mat[idLst,]-2

# inputs
nSamplingRep=1
ndmuRep=1
noffspringpercross <- 5
ntime=4
nAnimalEachPop <- c(300,250,100,1000)  # pop1,2,3,4 
# # GA run
nstartpop=150			# actual number will be nstartpop*2
maxiterOpt=10			# testing 10, otherwise 2000 or more

iRepScheme=1

for (iSampling in 1:nSamplingRep)  {
  
  ## recalculating aQtleffects, domValues & episPairLociEffect, based on statistical inputs given  
  # variables needed to run optimization algorithm: all global variables + genoFounderRes, genoFounderNRMS, genoinbredRes, genoinbredNRMS, geno2wayNRMS, geno3way, nstartpop, maxiterOpt
  # Method of variance calculation by indivdiual or by locus (Refer to the paper for the method)
  
  filename <- paste0(programdir,"S3GAopt3_byid.R")   # S3GAopt4_bylocus.R
  if (!file.exists(filename)) {  write(paste('Error: The file',filename,'does not exists. This module is needed for',filename),file=outlog,append=file.exists(outlog))}
  source(filename)
  
  ################################################################
  ## Breeding schemeGenerate 3-way hybrid indivdiuals ############
  ## breeding programs for rye plant in inbred populations########
  ################################################################
  # The breeding scheme requires program DMU to estimate Variance component & EBV.
  # However, DMU is not provided together with this code. Thus, the code stops here.
if (DMUavailable) {   
  for (dmuRep in 1:ndmuRep) {
    
    ########## Plot phenotypes
    keepVariables <- ls()   
    
    # reset master dataset
    maxid=nFounderAnimal*nFounderPop
    pop=pop[1:maxid,]
    haplo_pat=haplo_pat[1:maxid,]
    haplo_mat=haplo_mat[1:maxid,]
    
    for (currenttime in 1:ntime) {
      for (currentpop in 1:2) {   # currentpop=1; currenttime=1
        idLst1 <- rep(sample(which(pop$currentpop==currentpop & pop$fixedFactor1==currenttime-1),floor(nAnimalEachPop[currentpop]/noffspringpercross),replace=T),noffspringpercross)		# select animals for crossing
        idLst2 <- rep(sample(which(pop$currentpop==currentpop & pop$fixedFactor1==currenttime-1),floor(nAnimalEachPop[currentpop]/noffspringpercross),replace=T),noffspringpercross)		# select animals for crossing
        
        idpoptime <- c((maxid+1):(maxid+length(idLst1)))
        # prestep to prepare for mapply 
        pop[(nrow(pop)+1):(nrow(pop)+length(idLst1)),]=0
        haplo_pat <- rbind(haplo_pat,matrix(ncol=ncol(haplo_pat),nrow=length(idLst1)))
        haplo_mat <- rbind(haplo_mat,matrix(ncol=ncol(haplo_mat),nrow=length(idLst1)))
        maxid=maxid+length(idLst1)
        # F1: mapply sample all offsprings at once
        invisible(mapply(sampleOffspring,idd=idpoptime,par1=idLst1,par2=idLst2,fixedFactors=list(c(currenttime,0)),popid=currentpop))
        
        #F2-Fn
        for (i in 1:nselfing) { #i=1
          invisible(mapply(sampleOffspring,idd=idpoptime,par1=idpoptime,par2=idpoptime,fixedFactors=list(c(currenttime,i+1)),popid=currentpop))			# mapply sample all offsprings at once
        }	
      }  
      if (currenttime>1) {
        # make two-way hybrids  
        currentpop=3
        
        idLst1 <- sample(which(pop$currentpop==2 & pop$fixedFactor1==currenttime),nAnimalEachPop[currentpop],replace=T)  		# select NR
        idLst2 <- sample(which(pop$currentpop==2 & pop$fixedFactor1==currenttime-1),nAnimalEachPop[currentpop],replace=T) 		# select MS
        
        idpoptime <- c((maxid+1):(maxid+length(idLst1)))
        # prestep to prepare for mapply 
        pop[(nrow(pop)+1):(nrow(pop)+length(idLst1)),]=0
        haplo_pat <- rbind(haplo_pat,matrix(ncol=ncol(haplo_pat),nrow=length(idLst1)))
        haplo_mat <- rbind(haplo_mat,matrix(ncol=ncol(haplo_mat),nrow=length(idLst1)))
        maxid=maxid+length(idLst1)
        # F1: mapply sample all offsprings at once
        invisible(mapply(sampleOffspring,idd=idpoptime,par1=idLst1,par2=idLst2,fixedFactors=list(c(currenttime,0)),popid=currentpop))
        
        # make two-way hybrids  
        currentpop=4
        #idLst1 <- rep(sample(which(pop$currentpop==1 & pop$fixedFactor1==currenttime),floor(nAnimalEachPop[currentpop]/noffspringpercross),replace=T),noffspringpercross)		# select NR
        #idLst2 <- rep(sample(which(pop$currentpop==3 & pop$fixedFactor1==currenttime),floor(nAnimalEachPop[currentpop]/noffspringpercross),replace=T),noffspringpercross)		# select MS
        idLst1 <- sample(which(pop$currentpop==1 & pop$fixedFactor1==currenttime),nAnimalEachPop[currentpop],replace=T)  		# select NR
        idLst2 <- sample(which(pop$currentpop==3 & pop$fixedFactor1==currenttime),nAnimalEachPop[currentpop],replace=T) 	  	# select MS
        
        idpoptime <- c((maxid+1):(maxid+length(idLst1)))
        # prestep to prepare for mapply 
        pop[(nrow(pop)+1):(nrow(pop)+length(idLst1)),]=0
        haplo_pat <- rbind(haplo_pat,matrix(ncol=ncol(haplo_pat),nrow=length(idLst1)))
        haplo_mat <- rbind(haplo_mat,matrix(ncol=ncol(haplo_mat),nrow=length(idLst1)))
        maxid=maxid+length(idLst1)
        # F1: mapply sample all offsprings at once
        invisible(mapply(sampleOffspring,idd=idpoptime,par1=idLst1,par2=idLst2,fixedFactors=list(c(currenttime,0)),popid=currentpop))
        
      }  
    }
    
    #
    idplotobs <- pop[pop$fixedFactor1>0 & pop$currentpop %in% c(4),'id']  # use for phenotyping obs as plot 
    filename <- paste0(programdir,"S5TestCheckDMU8.R")
    if (!file.exists(filename)) {  write(paste('Error: The file',filename,'does not exists. This module is needed for',filename),file=outlog,append=file.exists(outlog))}
    source(filename)
    
    # reset new replicates
    rm(list=setdiff(ls(), keepVariables))
    
    
    ###### Given the same founder population and functional effects, check variances as if Individual phenotypes (instead of plot phenotypes).
    keepVariables <- ls()
    
    # reset master dataset
    maxid=nFounderAnimal*nFounderPop
    pop=pop[1:maxid,]
    haplo_pat=haplo_pat[1:maxid,]
    haplo_mat=haplo_mat[1:maxid,]
    
    for (currenttime in 1:ntime) {
      for (currentpop in 1:2) {   # currentpop=1; currenttime=1
        idLst1 <- rep(sample(which(pop$currentpop==currentpop & pop$fixedFactor1==currenttime-1),floor(nAnimalEachPop[currentpop]/noffspringpercross),replace=T),noffspringpercross)		# select animals for crossing
        idLst2 <- rep(sample(which(pop$currentpop==currentpop & pop$fixedFactor1==currenttime-1),floor(nAnimalEachPop[currentpop]/noffspringpercross),replace=T),noffspringpercross)		# select animals for crossing
        
        idpoptime <- c((maxid+1):(maxid+length(idLst1)))
        # prestep to prepare for mapply 
        pop[(nrow(pop)+1):(nrow(pop)+length(idLst1)),]=0
        haplo_pat <- rbind(haplo_pat,matrix(ncol=ncol(haplo_pat),nrow=length(idLst1)))
        haplo_mat <- rbind(haplo_mat,matrix(ncol=ncol(haplo_mat),nrow=length(idLst1)))
        maxid=maxid+length(idLst1)
        # F1: mapply sample all offsprings at once
        invisible(mapply(sampleOffspring,idd=idpoptime,par1=idLst1,par2=idLst2,fixedFactors=list(c(currenttime,0)),popid=currentpop))
        
        #F2-Fn
        for (i in 1:nselfing) { #i=1
          invisible(mapply(sampleOffspring,idd=idpoptime,par1=idpoptime,par2=idpoptime,fixedFactors=list(c(currenttime,i+1)),popid=currentpop))			# mapply sample all offsprings at once
        }	
      }  
      if (currenttime>1) {
        # make two-way hybrids  
        currentpop=3
        
        idLst1 <- sample(which(pop$currentpop==2 & pop$fixedFactor1==currenttime),nAnimalEachPop[currentpop],replace=T)  		# select NR
        idLst2 <- sample(which(pop$currentpop==2 & pop$fixedFactor1==currenttime-1),nAnimalEachPop[currentpop],replace=T) 		# select MS
        
        idpoptime <- c((maxid+1):(maxid+length(idLst1)))
        # prestep to prepare for mapply 
        pop[(nrow(pop)+1):(nrow(pop)+length(idLst1)),]=0
        haplo_pat <- rbind(haplo_pat,matrix(ncol=ncol(haplo_pat),nrow=length(idLst1)))
        haplo_mat <- rbind(haplo_mat,matrix(ncol=ncol(haplo_mat),nrow=length(idLst1)))
        maxid=maxid+length(idLst1)
        # F1: mapply sample all offsprings at once
        invisible(mapply(sampleOffspring,idd=idpoptime,par1=idLst1,par2=idLst2,fixedFactors=list(c(currenttime,0)),popid=currentpop))
        
        # make two-way hybrids  
        currentpop=4
        idLst1 <- rep(sample(which(pop$currentpop==1 & pop$fixedFactor1==currenttime),floor(nAnimalEachPop[currentpop]/noffspringpercross),replace=T),noffspringpercross)		# select NR
        idLst2 <- rep(sample(which(pop$currentpop==3 & pop$fixedFactor1==currenttime),floor(nAnimalEachPop[currentpop]/noffspringpercross),replace=T),noffspringpercross)		# select MS
        #idLst1 <- sample(which(pop$currentpop==1 & pop$fixedFactor1==currenttime),nAnimalEachPop[currentpop],replace=T)  		# select NR
        #idLst2 <- sample(which(pop$currentpop==3 & pop$fixedFactor1==currenttime),nAnimalEachPop[currentpop],replace=T) 	  	# select MS
        
        idpoptime <- c((maxid+1):(maxid+length(idLst1)))
        # prestep to prepare for mapply 
        pop[(nrow(pop)+1):(nrow(pop)+length(idLst1)),]=0
        haplo_pat <- rbind(haplo_pat,matrix(ncol=ncol(haplo_pat),nrow=length(idLst1)))
        haplo_mat <- rbind(haplo_mat,matrix(ncol=ncol(haplo_mat),nrow=length(idLst1)))
        maxid=maxid+length(idLst1)
        # F1: mapply sample all offsprings at once
        invisible(mapply(sampleOffspring,idd=idpoptime,par1=idLst1,par2=idLst2,fixedFactors=list(c(currenttime,0)),popid=currentpop))
        
      }  
    }
    
    # commented out, so individuals' observations 
    idplotobs <- NULL
    #idplotobs <- pop[pop$fixedFactor1>0 & pop$currentpop %in% c(4),'id']  # use for phenotyping obs as plot 
    filename <- paste0(programdir,"S5TestCheckDMU8Individuals.R")
    if (!file.exists(filename)) {  write(paste('Error: The file',filename,'does not exists. This module is needed for',filename),file=outlog,append=file.exists(outlog))}
    source(filename)
    
    # reset new replicates
	rm(list=setdiff(ls(), keepVariables))
    
    iRepScheme=iRepScheme+1
  }
}

}
