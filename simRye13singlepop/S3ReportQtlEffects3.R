if (!exists('iRepScheme')) iRepScheme=1

# Report distribution of solutions for biological effects
for (iTrait in 1:(nTrait)) {
  hist(aQtleffects[,iTrait],main=paste('Histogram of biological aQtleffects for Trait',iTrait))
  hist(episPairLociEffect[,iTrait],main=paste('Histogram of biological episPairLociEffect for Trait',iTrait))
  hist(domValues[,iTrait],main=paste('Histogram of biological domValues for Trait',iTrait))
  hist(domdegreeValues[,iTrait],main=paste('Histogram of biological domdegreeValues for Trait',iTrait))
}	
domdegreeValues <- c() # not use anymore, remove for computation efficiency.

# list of files to be reported.
# check if these genotype files are available. Remember these files must have the same dimensions
lstgenotypefiles <- c('genoFounderRes', 'genoFounderNRMS', 'genoinbredRes', 'genoinbredNRMS', 'geno2wayNRMS', 'geno3way')
temp1 <- c()
for (i in lstgenotypefiles) {
  if (!exists(i)) {      temp1 <- c(temp1,i)	  	}
}
lstgenotypefiles <- lstgenotypefiles[!lstgenotypefiles %in% temp1]
if (length(lstgenotypefiles)<1) stop (paste('No genotype file found for reporting.'))
if (length(dim(get(lstgenotypefiles[1])))!=2) stop (paste('wrong dimensions for the genotype file.'))
temp1 <- dim(get(lstgenotypefiles[1]))[1]
for (i in lstgenotypefiles) {
  if (dim(get(i))[1] !=temp1) stop (paste('wrong dimensions for the genotype file',i))
}

# report statistical variances of pops
filename='varFounder.res'
if (!file.exists(filename)) write(paste('iRep','pop','freqPop','effect','covname','method',paste0(paste0('val',c(1:(nTrait*nTrait))),collapse ='\t')),file=filename,append=F)

print('*** Statistical variances of pop based on current pop.***')
counter=0
for (ifile in lstgenotypefiles) {
  if (!exists(ifile)) next
  genoTemp = get(ifile)        # frequencyBasedPop
  temp0 <- calculateStatEffectOne(genodf=cbind(1:nrow(genoTemp),genoTemp),returnvalues='variances')
  print(paste('Pop',ifile,':'))
  print(temp0)
  temp0 <- data.frame(iRepScheme,pop=ifile,freqPop=ifile,temp0)
  write.table(temp0,file=filename,sep = "\t",row.names=F,col.names=!file.exists(filename),append=file.exists(filename))
  counter=counter+1
}
if (counter==0) print('*** No genotype files found for reporting variances. ***')

# report biological variances of pops
print('*** Biological variances of pop based on current pop.***')
counter=0
for (ifile in lstgenotypefiles) {
  if (!exists(ifile)) next
  genoTemp = get(ifile)        # frequencyBasedPop
  temp0 <- calculateBiolEffectOne(genodf=cbind(1:nrow(genoTemp),genoTemp),returnvalues='variances')
  print(paste('Pop',ifile,':'))
  print(temp0)
  temp0 <- data.frame(iRepScheme,pop=ifile,freqPop=ifile,temp0)
  write.table(temp0,file=filename,sep = "\t",row.names=F,col.names=!file.exists(filename),append=file.exists(filename))
  counter=counter+1
}
if (counter==0) print('*** No genotype files found for reporting variances. ***')

# report statistical variances of pops
print('*** Statistical variances of pop based on pop 3way.***')
counter=0
for (ifile in lstgenotypefiles[lstgenotypefiles!='geno3way']) {
  if (!exists(ifile) | !exists('geno3way')) next
  
  genoTemp = rbind(get(ifile),geno3way)        # frequencyBasedPop
  temp1 =cbind(1:nrow(genoTemp),genoTemp)
  temp2=(1:(dim(temp1)[1]/2)); temp3=(((dim(temp1)[1]/2)+1):dim(temp1)[1])
  temp0 <- calculateStatEffectOne(genodf=temp1,returnvalues='variances',idtobeused=temp2,idfrequencyCalculation=temp3)
  print(paste('Pop',ifile,':'))
  print(temp0)
  temp0 <- data.frame(iRepScheme,pop=ifile,freqPop='geno3way',temp0)
  write.table(temp0,file=filename,sep = "\t",row.names=F,col.names=!file.exists(filename),append=file.exists(filename))
  counter=counter+1
}
if (counter==0) print('*** No genotype files found for reporting variances. ***')

#### plot
if (exists('genoinbredRes') & exists('geno2wayNRMS')) { 
  trueVarCalPlot <- calculateStatEffectPlotOne(genoPar1=genoinbredRes,genoPar2=geno2wayNRMS,returnvalues='variances')
  trueVarCalPlot <- data.frame(iRep=iRepScheme,trueVarCalPlot)
  filename='varSchemTrueAvCPlotOne.res'
  write.table(trueVarCalPlot,file=filename,sep = "\t",row.names=F,col.names=!file.exists(filename),append=file.exists(filename))
}

#dev.off()
print('Report completed.')

######################################################################################################################
########## Update phenotype of base animals with new aQtleffects, domValues & episPairLociEffect generated ###########
######################################################################################################################

for (idd in 1:maxid) { # idd=1
  
  # true genetic values of aniaml id
  temp1=haplo_pat[idd,SNPinfo[,4]=='T']+haplo_mat[idd,SNPinfo[,4]=='T']-3
  for (iTrait in 1:(nTrait)) { #iTrait=1
    pop[idd,paste0('trait_av',iTrait)]=sum(aQtleffects[,iTrait] * temp1)
    pop[idd,paste0('trait_dv',iTrait)]=sum(domValues[temp1==0,iTrait]) 
  }
  
  # true epistasis genetic values of aniaml id
  if (!is.null(episPairLociEffect)) {
    
    temp1 <- haplo_pat[idd,episPairLoci[,1]]
    temp2 <- haplo_mat[idd,episPairLoci[,1]]
    xAk <- temp1 + temp2 - 3                  # scaled genotype dosages at locus k
    temp1 <- haplo_pat[idd,episPairLoci[,2]]
    temp2 <- haplo_mat[idd,episPairLoci[,2]]
    xAl <- temp1 + temp2 - 3                  # scaled genotype dosages at locus l
    for (iTrait in 1:(nTrait)) {
      pop[idd,paste0('trait_episgv',iTrait)]=sum(episPairLociEffect[1:nEpisPairLoci,iTrait]*xAk*xAl)
    }	  
  }  
  
  # sample residuals effects of aniaml id
  #pop[idd,paste0('trait_residual',c(1:nRes))]=mvrnorm(1,rep(0,nRes),varResiduals)  # not needed as it has been done before.
  
  # trait_tgv=trait_av+trait_dv+trait_episgv
  for (iTrait in 1:(nTrait)) {  
    pop[idd,paste0('trait_tgv',iTrait)]=pop[idd,paste0('trait_av',iTrait)]+
      pop[idd,paste0('trait_dv',iTrait)]+pop[idd,paste0('trait_episgv',iTrait)]
  }
  # sample obs. for now, obs is simply=trait_tgv+trait_residual
  for (iTrait in 1:(nObs)) {  
    pop[idd,paste0('trait_obs',iTrait)]=pop[idd,paste0('trait_tgv',iTrait)]+
      pop[idd,paste0('trait_residual',iTrait)]
  }
}
