########           generate Base Pops                   ##########
##################################################################
haplo_pat=haplo_mat=matrix(NA,nrow=nFounderAnimal*nFounderPop,ncol=nSnp)
haplo_pat[,]=allPop[seq(1,(nFounderAnimal*nFounderPop*2),by=2),]   # haplo_pat hold paternal genotypes of all animals
haplo_mat[,]=allPop[seq(2,(nFounderAnimal*nFounderPop*2),by=2),]   # haplo_mat hold maternal genotypes of all animals

# set up base population
maxid=nFounderAnimal*nFounderPop

pop <- data.frame(id=1:maxid,par1=0,par2=0,birthpop=0,currentpop=0)
for (iTrait in 1:(nTrait)) {
  pop[,paste0('trait_av',iTrait)]=0
  pop[,paste0('trait_dv',iTrait)]=0
  pop[,paste0('trait_tgv',iTrait)]=0	
  pop[,paste0('trait_episgv',iTrait)]=0	
}
for (iTrait in 1:(nObs)) {
  pop[,paste0('trait_obs',iTrait)]=0
}
for (iTrait in 1:(nEbv)) {
  pop[,paste0('trait_polyEbv',iTrait)]=0
  pop[,paste0('trait_genoEbv',iTrait)]=0
}
for (iTrait in 1:(nRes)) {
  pop[,paste0('trait_residual',iTrait)]=0
}
#add other factors
for (iFactor in 1:nfixedFactors) {
  pop[,paste0('fixedFactor',iFactor)]=0
}
for (iFactor in 1:nranFactors) {
  pop[,paste0('randomFactor',iFactor)]=-999.0
}

for (idd in 1:maxid) { # idd=1
  pop$birthpop[idd]=pop$currentpop[idd]=ceiling(idd/nFounderAnimal)
  pop$par1[idd]=pop$currentpop[idd]*(-1)
  pop$par2[idd]=pop$currentpop[idd]*(-1)
  
  # true genetic values of aniaml id
  iQtl=0
  temp1=haplo_pat[idd,SNPinfo[,4]=='T']+haplo_mat[idd,SNPinfo[,4]=='T']-3
  for (iTrait in 1:(nTrait)) { #iTrait=1
      pop[idd,paste0('trait_av',iTrait)]=pop[idd,paste0('trait_av',iTrait)]+
	    sum(aQtleffects[,iTrait] * temp1)
	  pop[idd,paste0('trait_dv',iTrait)]=pop[idd,paste0('trait_dv',iTrait)]+
	   sum(domValues[temp1==0,iTrait]) 
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
          pop[idd,paste0('trait_episgv',iTrait)]=pop[idd,paste0('trait_episgv',iTrait)]+
            sum(episPairLociEffect[1:nEpisPairLoci,iTrait]*xAk*xAl)
      }
	  
  }  
  
  # sample residuals effects of aniaml id
  pop[idd,paste0('trait_residual',c(1:nRes))]=mvrnorm(1,rep(0,nRes),varResiduals)

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

# remove unneeded variables
keepVariables <- c(keepVariables,'maxid','pop','haplo_pat','haplo_mat')
rm(list=setdiff(ls(), keepVariables))
