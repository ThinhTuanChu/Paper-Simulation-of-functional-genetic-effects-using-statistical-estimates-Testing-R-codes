######### Starting population, calculate genetic effects ############
#####################################################################
# Simulate functional effects, Using GA method

#list of library packages required
temp1 <- c('MASS','GA') 			# c("readxl","plyr","ggplot2",'GA')
temp2 <- temp1 %in% rownames(installed.packages()) # Install packages not yet installed
if (any(temp2 == FALSE)) {
  install.packages(temp1[!temp2])}
#load packages
invisible(lapply(temp1, library, character.only = TRUE))

################### Pre-requisite ##########################
## need many functions from S1CalculateQtlEffects*.R #######

#####################################################################
############# Function(s) to simulate genetic effects ###############
# required global variables: nQtl, nTrait, covarTraits, SNPinfo, episcovarTraits,domcovarTraits, domdegreemean, domdegreevar, nEpisPairLoci, episPairLoci
# check if these genotype files are available. Remember these files must have the same dimensions
lstgenotypefiles <- c('genoFounderRes', 'genoFounderNRMS', 'genoinbredRes', 'genoinbredNRMS', 'geno2wayNRMS', 'geno3way')
if (length(dim(get(lstgenotypefiles[1])))!=2) stop (paste('wrong dimensions for the genotype file.'))
temp1 <- dim(get(lstgenotypefiles[1]))[1]
for (i in lstgenotypefiles) {
  if (!exists(i)) stop (paste('genodf file',i,'does not exists. This genotype file is required.'))
  if (dim(get(i))[1] !=temp1) stop (paste('wrong dimensions for the genotype file',i))
}
# genodf=Pops[[2]][seq(1,200,by=2),]+Pops[[2]][seq(2,200,by=2),]-2
# covarTraits=matrix(c(400, 0.00,0.0, 400),nrow = nTrait,byrow = T)     # testing
# episcovarTraits=matrix(c(100, 0.00,0.0, 100),nrow = nTrait,byrow = T) # testing
# domcovarTraits=matrix(c(100, 0.00,0.0, 100),nrow = nTrait,byrow = T)  # testing
# idtobeused=NULL;idfrequencyCalculation=NULL;infuncParameters='biological';reportresult=TRUE;domdegreemean=NULL;domdegreevar=NULL

generateStartPop <- function(genodf=NULL,idtobeused=NULL,idfrequencyCalculation=NULL,infuncParameters='biological',nIndividuals=100,
                             covarTraits=NULL,episcovarTraits=NULL,domcovarTraits=NULL,domdegreemean=NULL,domdegreevar=NULL,
                             methodCalculateVar='byid') {
  #genodf=NULL;idtobeused=NULL;idfrequencyCalculation=NULL;infuncParameters='biological';reportresult=TRUE;
  #covarTraits=NULL;episcovarTraits=NULL;domcovarTraits=NULL;domdegreemean=NULL;domdegreevar=NULL
  #genodf=geno3way;infuncParameters='statisticalTwo'
  if (is.null(genodf)) {  stop ('input: genodf is needed to run this function.') }
  if (length(dim(genodf))!=2|(length(genodf[1,]))<1 | (length(genodf[,1]))<1) {  stop ('input: genodf could be empty or invalid matrix.') }
  genodf <- as.matrix(genodf)   
  if ((length(genodf[1,]))==nSnp) genodf = cbind(matrix(c(1:length(genodf[,1])),nrow=length(genodf[,1])),genodf)  # add column ID
  if (is.null(idtobeused))             {  idtobeused = genodf[,1]           }
  if (is.null(idfrequencyCalculation)) {  idfrequencyCalculation=idtobeused }
  if ((length(genodf[,1]))<length(unique(c(idtobeused,idfrequencyCalculation)))) stop(paste('Dimension of genodf is not right. Nrow must be at least:',length(unique(idtobeused,idfrequencyCalculation))))
  if ((length(genodf[1,])-1)!=nSnp) stop(paste('Dimension of genodf is not right. Ncolumn must be:',nSnp,'or',nSnp+1))
  if (any(!genodf[,2:length(genodf[1,])] %in% c(0,1,2))) stop('Genotype codes are not in correct format. It must be 0,1,2.')
  if (length(idtobeused)!=length(unique(idtobeused))) {
    idtobeused=unique(idtobeused); 	  print('Repeated ID are removed.')
  }
  # variances input
  if (is.null(covarTraits)) stop ('input: covarTraits is needed to run this function.')
  domsimulated <- TRUE; epissimulated <- TRUE
  if (is.null(domdegreemean) & is.null(domcovarTraits) & is.null(domdegreevar)) {   domsimulated <- FALSE  }
  if (is.null(episcovarTraits)) {   epissimulated <- FALSE  }
  if (is.null(domdegreevar))  {   domdegreevar <- diag(rep(sqrt(0.0973),nTrait))  }
  if (is.null(domdegreemean)) {   domdegreemean <- rep(0.193,nTrait)  }
  if (is.null(domdegreemean) & is.null(domcovarTraits) & is.null(domdegreevar)) {   domdegreemean <- rep(0.193,nTrait)  }
  
  # keep needed genotypes only
  allID <- unique(c(idtobeused,idfrequencyCalculation))
  allID <- allID[order(allID)]   
  genodf <- genodf[genodf[,1] %in% allID,]
  
  # set rownames of geno as id 
  geno <- genodf[,-1]
  rownames(geno) <- genodf[,1] ; genodf <- NULL
  idtobeused=as.character(idtobeused)
  nAnimal=length(idtobeused)
  idfrequencyCalculation=as.character(idfrequencyCalculation)
  
  ### calculate frequency of alleles and genotypes
  genoTemp <- geno[idfrequencyCalculation,]
  # alleles
  freqAlleles <- array(data = 0,dim=c(nSnp,2))
  MAF <- function(x){sum(x,na.rm=T)/(sum(!is.na(x))*2)}
  freqAlleles[,2]=apply(genoTemp, 2, MAF)   # q
  freqAlleles[,1]=1-freqAlleles[,2]         # p
  # genotypes
  freqGenotypes <- array(data = 0,dim=c(nSnp,3))
  GenoFreq <- function(x,y){sum(x==y)/length(x)}
  counter=0
  for (iAllele in 1:2) {
    for (jAllele in iAllele:2) {
      counter=counter+1
      freqGenotypes[,counter] <- apply(genoTemp, 2, GenoFreq,y=(counter-1))
    }
  }
  genoTemp <- NULL
  geno <- geno[idtobeused,]
  rownames(geno) <- NULL   # No need rownames
  
  # start sampling
  nLociforAlgorithm= nTrait*nQtl +  nTrait*nQtl + nTrait*nEpisPairLoci
  startpop <- matrix(NA, nrow = nIndividuals, ncol = nLociforAlgorithm)
  
  for (idstartpop in 1:nIndividuals) {  
    ##### simulate prior for epistasis effects
    episPairLociEffect <- array(data = 0,dim=c(nEpisPairLoci,nTrait))	
    if (epissimulated) {
      if (sum(diag(episcovarTraits))>0) {
        
        episPairLociEffect[,] <- mvrnorm(nEpisPairLoci,rep(0,nTrait),episcovarTraits/nEpisPairLoci)
        # if epis dominance not simulated, episPairLociEffect = statistical_epis, 
        # but the matrix to link to true biological phenotype and estimate phenotype are different between two.
        # The following is for estimate phenotype for variance components...
        
        # calculate variances by animals
        if (methodCalculateVar=='byid') {
          idgv <- array(data = 0,dim=c(nAnimal,nTrait))
          if (infuncParameters %in% c('statisticalOne','statisticalTwo','statisticalThree')) {
            temp1 <- freqGenotypes[episPairLoci[,1],]
            pxAk <- lapply(seq_len(nrow(temp1)), function(i) temp1[i,])
            temp1 <- freqGenotypes[episPairLoci[,2],]
            pxAl <- lapply(seq_len(nrow(temp1)), function(i) temp1[i,])
            for (idd in 1:nAnimal){
              xAk <- geno[idd,episPairLoci[,1]]
              xAl <- geno[idd,episPairLoci[,2]]
              haihaiVect <- mapply(findhaihai,pk=pxAk,ack=xAk,pl=pxAl,acl=xAl)		
              for (iTrait in 1:(nTrait)) {
                idgv[idd,iTrait]=sum(haihaiVect*episPairLociEffect[,iTrait])    # episPairLociEffect = statistical_epis
              }
            }
          }
          if (infuncParameters %in% c('biological')) {
            for (idd in 1:nAnimal){
              xAk <- geno[idd,episPairLoci[,1]]-1
              xAl <- geno[idd,episPairLoci[,2]]-1
              for (iTrait in 1:(nTrait)) {
                idgv[idd,iTrait]=sum(episPairLociEffect[1:nEpisPairLoci,iTrait]*xAk*xAl)
              }
            }
          }
          v_temp=var(idgv)   # calculate covariances for traits
        }
        
        if (methodCalculateVar=='bylocus') {
          if (infuncParameters %in% c('biological')) {
            varbylocus <- array(data = 0,dim=c(nTrait,nTrait,nEpisPairLoci))
            tempfunc <- function(itemp) {
              varbylocus[,,itemp] <<- calculateVarbyLocusEpis(EffInput=episPairLociEffect[itemp,],
                                                              pk=freqGenotypes[episPairLoci[itemp,1],],pl=freqGenotypes[episPairLoci[itemp,2],],
                                                              hk=c(-1,0,1), hl=c(-1,0,1))
            }
            invisible(mapply(tempfunc,itemp=1:nEpisPairLoci))    # using apply
            v_temp=array(data = 0,dim=c(nTrait,nTrait))
            for (iTrait in 1:(nTrait)) {
              for (jTrait in 1:(nTrait)) {
                v_temp[iTrait,jTrait]=sum(varbylocus[iTrait,jTrait,])
              }
            } 
          }
          if (infuncParameters %in% c('statisticalOne','statisticalTwo','statisticalThree')) {
            
            varbylocus <- array(data = 0,dim=c(nTrait,nTrait,nEpisPairLoci))
            tempfunc <- function(itemp) {
              varbylocus[,,itemp] <<- calculateVarbyLocusEpis(EffInput=episPairLociEffect[itemp,],
                                                              pk=freqGenotypes[episPairLoci[itemp,1],],pl=freqGenotypes[episPairLoci[itemp,2],],
                                                              hk=c(findhai(px=freqGenotypes[episPairLoci[itemp,1],],allelecount=0),
                                                                   findhai(px=freqGenotypes[episPairLoci[itemp,1],],allelecount=1),
                                                                   findhai(px=freqGenotypes[episPairLoci[itemp,1],],allelecount=2)),
                                                              hl=c(findhai(px=freqGenotypes[episPairLoci[itemp,2],],allelecount=0),
                                                                   findhai(px=freqGenotypes[episPairLoci[itemp,2],],allelecount=1),
                                                                   findhai(px=freqGenotypes[episPairLoci[itemp,2],],allelecount=2)))
            }
            invisible(mapply(tempfunc,itemp=1:nEpisPairLoci))    # using apply
            v_temp=array(data = 0,dim=c(nTrait,nTrait))
            for (iTrait in 1:(nTrait)) {
              for (jTrait in 1:(nTrait)) {
                v_temp[iTrait,jTrait]=sum(varbylocus[iTrait,jTrait,])
              }
            } 	
          }
        }
        
        #cov2cor(v_temp); cov2cor(episcovarTraits); episcovarTraits
        Vmat=solve(v_temp)
        Smat <- episcovarTraits %*% Vmat
        # square root matrix 
        ei <- eigen(Smat)
        Vect <- ei$vectors
        res <- Vect %*% diag(sqrt(ei$values)) %*% solve(Vect) # 
        
        # rescale epistasis effects 
        tempRescale = episPairLociEffect
        rescaleFun <- function(iEpisPairLoci) {
          episPairLociEffect[iEpisPairLoci,] <<- res %*% tempRescale[iEpisPairLoci,]}
        invisible(mapply(rescaleFun,iEpisPairLoci=1:nEpisPairLoci)) # for (iEpisPairLoci in 1:nEpisPairLoci) {   episPairLociEffect[iEpisPairLoci,] <- res %*% episPairLociEffect[iEpisPairLoci,]  }   # using apply instead
      }  
    }
    
    ##### simulate prior for a effects
    aQtleffects <- array(data = 0,dim=c(nQtl,nTrait))
    aQtleffects[1:nQtl,] <- mvrnorm(nQtl,rep(0,nTrait),covarTraits/nQtl)  # nr sampling, mean vector, covariance matrix
    
    # calculate tbv by animals
    if (methodCalculateVar=='byid') {
      idgv <- array(data = 0,dim=c(nAnimal,nTrait)) 
      if (infuncParameters %in% c('statisticalOne','statisticalTwo','statisticalThree')) {
        temp1 <- freqGenotypes[SNPinfo[,4]=='T',]
        px <- lapply(seq_len(nrow(temp1)), function(i) temp1[i,])
        for (idd in 1:nAnimal){
          genotypecode=geno[idd,SNPinfo[,4]=='T']
          haiVect <- mapply(findhai,px=px,allelecount=genotypecode)
          for (iTrait in 1:(nTrait)) {
            idgv[idd,iTrait] <- sum(haiVect*aQtleffects[,iTrait])
          }	
        }
      }
      if (infuncParameters %in% c('biological')) {
        
        for (idd in 1:nAnimal){
          genotypecode=geno[idd,SNPinfo[,4]=='T']
          for (iTrait in 1:(nTrait)) {
            idgv[idd,iTrait]=sum(aQtleffects[,iTrait] * (genotypecode-1))
          }
        }
      }
      v_temp=var(idgv)   # calculate covariances for traits
    }
    
    if (methodCalculateVar=='bylocus') {
      
      if (infuncParameters %in% c('biological')) {
        varbylocus <- array(data = 0,dim=c(nTrait,nTrait,nQtl))
        px <- freqGenotypes[SNPinfo[,4]=='T',]
        tempfunc <- function(itemp) {
          varbylocus[,,itemp] <<- calculateVarbyLocus(EffInput=aQtleffects[itemp,],px=px[itemp,],hx=c(-1,0,1))
        }
        invisible(mapply(tempfunc,itemp=1:nQtl))    # using apply
        v_temp=array(data = 0,dim=c(nTrait,nTrait))
        for (iTrait in 1:(nTrait)) {
          for (jTrait in 1:(nTrait)) {
            v_temp[iTrait,jTrait]=sum(varbylocus[iTrait,jTrait,])
          }
        } 
      }
      
      if (infuncParameters %in% c('statisticalOne','statisticalTwo','statisticalThree')) {
        varbylocus <- array(data = 0,dim=c(nTrait,nTrait,nQtl))
        px <- freqGenotypes[SNPinfo[,4]=='T',]
        tempfunc <- function(itemp) {
          varbylocus[,,itemp] <<- calculateVarbyLocus(EffInput=aQtleffects[itemp,],px=px[itemp,],
                                                      hx=c(findhai(px=px[itemp,],allelecount=0),findhai(px=px[itemp,],allelecount=1),findhai(px=px[itemp,],allelecount=2)))
        }
        invisible(mapply(tempfunc,itemp=1:nQtl))    # using apply
        v_temp=array(data = 0,dim=c(nTrait,nTrait))
        for (iTrait in 1:(nTrait)) {
          for (jTrait in 1:(nTrait)) {
            v_temp[iTrait,jTrait]=sum(varbylocus[iTrait,jTrait,])
          }
        } 
      }	
    }
    
    #cov2cor(v_temp); cov2cor(covarTraits); covarTraits
    Vmat=solve(v_temp)
    Smat <- covarTraits %*% Vmat
    # square root matrix 
    ei <- eigen(Smat)
    Vect <- ei$vectors
    res <- Vect %*% diag(sqrt(ei$values)) %*% solve(Vect) # 
    
    # rescale QTl effects 
    tempRescale = aQtleffects
    rescaleFun <- function(iQtl) {
      aQtleffects[iQtl,] <<- res %*% tempRescale[iQtl,]}
    invisible(mapply(rescaleFun,iQtl=1:nQtl)) # for (iQtl in 1:nQtl) {      aQtleffects[iQtl,] <- res %*% aQtleffects[iQtl,]     }   # using apply instead
    
    ##### simulate prior for d effects.
    domValues <- array(data = 0,dim=c(nQtl,nTrait))
    domdegreeValues <- array(data = 0,dim=c(nQtl,nTrait))
    if (domsimulated) {
      
      domdegreeValues[,] <- mvrnorm(nQtl,domdegreemean,domdegreevar)
      
      # caculate biological dominance from domdegree
      for (iTrait in 1:nTrait) {
        domValues[,iTrait] <- abs(aQtleffects[1:nQtl,iTrait])*domdegreeValues[1:nQtl,iTrait]
      }
      
      if (!is.null(domcovarTraits)) {
        
        # calculate dv by animals
        if (methodCalculateVar=='byid') {
          idgv <- array(data = 0,dim=c(nAnimal,nTrait)) 
          if (infuncParameters %in% c('statisticalOne','statisticalTwo','statisticalThree')) {
            temp1 <- freqGenotypes[SNPinfo[,4]=='T',]
            px <- lapply(seq_len(nrow(temp1)), function(i) temp1[i,])
            for (idd in 1:nAnimal){
              genotypecode=geno[idd,SNPinfo[,4]=='T']
              hdiVect <- mapply(findhdi,px=px,allelecount=genotypecode)
              
              for (iTrait in 1:(nTrait)) {
                idgv[idd,iTrait] <- sum(hdiVect*domValues[,iTrait])
              }
            }
          }
          if (infuncParameters %in% c('biological')) {
            for (idd in 1:nAnimal){
              genotypecode=geno[idd,SNPinfo[,4]=='T']
              for (iTrait in 1:(nTrait)) {
                idgv[idd,iTrait]=sum(domValues[genotypecode==1,iTrait])
              }	  
            }
          }
          v_temp=var(idgv)   # calculate covariances for traits
        }
        if (methodCalculateVar=='bylocus') {
          if (infuncParameters %in% c('biological')) {
            varbylocus <- array(data = 0,dim=c(nTrait,nTrait,nQtl))
            px <- freqGenotypes[SNPinfo[,4]=='T',]
            tempfunc <- function(itemp) {
              varbylocus[,,itemp] <<- calculateVarbyLocus(EffInput=domValues[itemp,],px=px[itemp,],hx=c(0,1,0))
            }
            invisible(mapply(tempfunc,itemp=1:nQtl))    # using apply
            v_temp=array(data = 0,dim=c(nTrait,nTrait))
            for (iTrait in 1:(nTrait)) {
              for (jTrait in 1:(nTrait)) {
                v_temp[iTrait,jTrait]=sum(varbylocus[iTrait,jTrait,])
              }
            } 
          }
          if (infuncParameters %in% c('statisticalOne','statisticalTwo','statisticalThree')) {
            varbylocus <- array(data = 0,dim=c(nTrait,nTrait,nQtl))
            px <- freqGenotypes[SNPinfo[,4]=='T',]
            tempfunc <- function(itemp) {
              varbylocus[,,itemp] <<- calculateVarbyLocus(EffInput=domValues[itemp,],px=px[itemp,],
                                                          hx=c(findhdi(px=px[itemp,],allelecount=0),findhdi(px=px[itemp,],allelecount=1),findhdi(px=px[itemp,],allelecount=2)))
            }
            invisible(mapply(tempfunc,itemp=1:nQtl))    # using apply
            v_temp=array(data = 0,dim=c(nTrait,nTrait))
            for (iTrait in 1:(nTrait)) {
              for (jTrait in 1:(nTrait)) {
                v_temp[iTrait,jTrait]=sum(varbylocus[iTrait,jTrait,])
              }
            } 
          }	  
        }
        # rescaling
        Vmat=solve(v_temp)
        Smat <- domcovarTraits %*% Vmat
        # square root matrix 
        ei <- eigen(Smat)
        Vect <- ei$vectors
        res <- Vect %*% diag(sqrt(ei$values)) %*% solve(Vect) # 
        # rescale QTl effects 
        tempRescale = domValues
        rescaleFun <- function(iQtl) {
          domValues[iQtl,] <<- res %*% tempRescale[iQtl,]}
        invisible(mapply(rescaleFun,iQtl=1:nQtl)) # for (iQtl in 1:nQtl) {      domValues[iQtl,] <- res %*% domValues[iQtl,]     }   # using apply instead
        
        # recalculate domdegreeValues given new rescaled domValues 
        #    for (iTrait in 1:nTrait) {
        #      temp1=abs(aQtleffects[1:nQtl,iTrait])
        #      temp1[temp1<0.00001]=0.00001   # avoid zero numbers lead to infinite values
        #      domdegreeValues[,iTrait] <- domValues[1:nQtl,iTrait]/temp1
        #    } 
      }
    }
    
    #### Calculate functional effects from statistical effects
    if (infuncParameters %in% c('statisticalOne')) {
      scaleSucceeded <- TRUE  # define variable
      
      # calculate statisticalepis_a. We need statisticalepis_a to rescale aQtleffects because effects from episPairLociEffect change aQtleffects no matter 'biological' or 'statisticalOne'
      statisticalepis_a <- array(data = 0,dim=c(nQtl,nTrait))
      
      # function to set up bkl 9x1 matrix (see Duenk 2020) to bkl of all pairs. This function requires b_kl array to be defined. 
      tepiscontrast <- kronecker(matrix(c(-1,0,1),nrow=3,ncol = 1),matrix(c(-1,0,1),nrow=3,ncol = 1))
      b_kl <- array(data = NA,dim=c(nrow(tepiscontrast),nEpisPairLoci,nTrait))	
      
      setupbklmat <- function(traitno,pairno,pk,pl) {
        ckl=tepiscontrast*episPairLociEffect[pairno,traitno]
        #Dkl=setupDklmat(pk,pl)       #Dkl is canceled out in the equation. So, it is not important if it is present or not. Will remove them later.
        
        Wk <- Wl <- matrix(nrow=3,ncol=3)
        Wk[,1]= Wl[,1]= 1.0
        Wk[,2]=c((0.0-pk[2]-2.0*pk[3]),
                 (1.0-pk[2]-2.0*pk[3]),
                 (2.0-pk[2]-2.0*pk[3]))
        Wk[,3]= Wl[,3]=0  # no dominance simulated. So, using 0 instead of below column.
        Wl[,2]=c((0.0-pl[2]-2.0*pl[3]),
                 (1.0-pl[2]-2.0*pl[3]),
                 (2.0-pl[2]-2.0*pl[3]))
        
        Wkl <- kronecker(Wl,Wk) # Wl%x%Wk #
        
        temp1 <- t(Wkl)%*%Wkl
        if (det(temp1)<0.0000000001) diag(temp1)=diag(temp1)+0.00000001   # this might be needed for non-invertable matrix
        bkl=solve(temp1)%*%t(Wkl)%*%ckl
        b_kl[,pairno,traitno] <<- bkl
      }
      
      # calculate alpha from epistasis effects
      # preparation for apply functions
      temp1 <- freqGenotypes[episPairLoci[,1],]
      pk <- lapply(seq_len(nrow(temp1)), function(i) temp1[i,])
      temp1 <- freqGenotypes[episPairLoci[,2],]
      pl <- lapply(seq_len(nrow(temp1)), function(i) temp1[i,])
      # apply funtions for caculate bkl
      for (iTrait in 1:(nTrait)) {
        invisible(mapply(setupbklmat,traitno=rep(iTrait,nEpisPairLoci),pairno=1:nEpisPairLoci,pk=pk,pl=pl))
      }
      
      for (iTrait in 1:(nTrait)) {
        temp1=rep(0,nSnp)
        temp1[episPairLoci[,1]]=temp1[episPairLoci[,1]]+b_kl[2,,iTrait] #+b_kl[1,,iTrait]
        temp1[episPairLoci[,2]]=temp1[episPairLoci[,2]]+b_kl[4,,iTrait] #+b_kl[1,,iTrait]
        statisticalepis_a[,iTrait]=statisticalepis_a[,iTrait] + temp1[SNPinfo[,4]=='T']
      }	    
      
      ## With alpha=a+d(1-2p)+epis or (rewrite: alpha=a+d1+e1), find Va by solving following equation for variance:
      #  V(alpha)=V(a) + 2*sqrt(V(a))(cor(a,d1)*sqrt(V(d1)) + cor(a,e)*sqrt(V(e1))) + V(d1) + 2*cor(d1,e1)*sqrt(V(d1)*V(e1)) + V(e1).   # Note: this formula derived because cor is constant if a rescaled.
      V_alpha=covarTraits
      # calculate V(d(1-2p)), V(a) and V(statisticalepis_a): additive genetic variance due to domValues[,iTrait]*(1.0-2.0*freqAlleles[SNPinfo[,4]=='T',1])
      V_d12p = V_a = V_statisticalepis_a = array(data = 0,dim=c(nTrait,nTrait))
      idgv_d12p <- idgv_a <- idgv_statisticalepis_a <- array(data = 0,dim=c(nAnimal,nTrait)) 
      
      temp1 <- freqGenotypes[SNPinfo[,4]=='T',]
      px <- lapply(seq_len(nrow(temp1)), function(i) temp1[i,])
      for (idd in 1:nAnimal){
        genotypecode=geno[idd,SNPinfo[,4]=='T']
        haiVect <- mapply(findhai,px=px,allelecount=genotypecode)
        for (iTrait in 1:(nTrait)) {
          idgv_d12p[idd,iTrait] <- sum(haiVect*domValues[,iTrait]*(1.0-2.0*freqAlleles[SNPinfo[,4]=='T',1]))
          idgv_a[idd,iTrait] <- sum(haiVect*aQtleffects[,iTrait])
          idgv_statisticalepis_a[idd,iTrait] <- sum(haiVect*statisticalepis_a[,iTrait])
        }	
      }
      V_d12p=var(idgv_d12p)   # calculate covariances for traits
      V_a=var(idgv_a)
      V_statisticalepis_a=var(idgv_statisticalepis_a)
      
      # cor(a,d(1-2p))
      cor_ad12p = cor(idgv_d12p,idgv_a)
      cor_astatisticalepis_a = cor(idgv_statisticalepis_a,idgv_a)
      cor_d12pstatisticalepis_a = cor(idgv_d12p,idgv_statisticalepis_a)
      
      # solve the equation. For now, only Variance rescaled, but not covariance, meaning cor are kept. May do in matrix way, if needed.
      tempa=1
      for (iTrait in 1:(nTrait)) {
        tempb=2*cor_ad12p[iTrait,iTrait]*sqrt(V_d12p[iTrait,iTrait])+ 2*cor_astatisticalepis_a[iTrait,iTrait]*sqrt(V_statisticalepis_a[iTrait,iTrait])
        tempc=V_d12p[iTrait,iTrait] - V_alpha[iTrait,iTrait] + V_statisticalepis_a[iTrait,iTrait] + 
          2*cor_d12pstatisticalepis_a[iTrait,iTrait]*sqrt(V_d12p[iTrait,iTrait]*V_statisticalepis_a[iTrait,iTrait])
        tempd=tempb^2-4*tempa*tempc  # this needed to be >0
        if (tempd <0) {
          scaleSucceeded <- FALSE
          break
        }
        x0=c((-tempb-sqrt(tempd))/(2*tempa),(-tempb+sqrt(tempd))/(2*tempa))
        x0=x0[x0>0]            #; print(paste('X0=',x0))
        if (length(x0)>0) {    # take 1 valid result only
          x0=x0[1]
          aQtleffects[,iTrait]=aQtleffects[,iTrait]*x0/sqrt(V_a[iTrait,iTrait])
        } else {
          scaleSucceeded <- FALSE
          break
        }
      }
      
      if (!scaleSucceeded) {print('Scaling to meet statistical variance of additive genetics was not successful.')}
      if (methodCalculateVar=='bylocus') { print('methodCalculateVar of bylocus does not work with infuncParameters of statisticalOne. methodCalculateVar=byid is automatically used.') }
    }
    
    #### Calculate functional effects from statistical effects
    if (infuncParameters %in% c('statisticalTwo')) {
      
      # setup bkl using ivares-Castro and Calborg2007 method, to calculate functional genetic effects from statistical effects of a, d and e
      # only work if (nQtlEachQtlInteractwith==1)
      if (nQtlEachQtlInteractwith!=1) stop ('this method works with nQtlEachQtlInteractwith=1 only.') 
      
      bkl <- array(data = NA,dim=c(6,nEpisPairLoci,nTrait))
      Wfkl <- setupWfklmat1()  # put here so code more efficient
      
      setupbkl <- function(traitno,pairno,pk,pl) {
        #traitno=1;pairno=1;pk=freqGenotypes[episPairLoci[pairno,1],];pl=freqGenotypes[episPairLoci[pairno,2],]
        
        qtlno_k=sum(SNPinfo[1:episPairLoci[pairno,1],4]=='T')
        qtlno_l=sum(SNPinfo[1:episPairLoci[pairno,2],4]=='T')
        
        # setup Es manually like this or do it similarly as in setup Wsklmat or Wfkl for equivalent positions
        Es=matrix(0,nrow = 6,ncol=1)
        Es[,1]=c(1,aQtleffects[qtlno_k,traitno],aQtleffects[qtlno_l,traitno],domValues[qtlno_k,traitno],domValues[qtlno_l,traitno],episPairLociEffect[pairno,traitno])
        # Warning note: need to check if mean affect Ef
        
        Wskl <- setupWsklmat(pk,pl)
        #Wfkl <- setupWfklmat1()
        
        temp1 <- t(Wfkl)%*%Wfkl
        if (det(temp1)<0.0000000001) diag(temp1)=diag(temp1)+0.00000001   # this might be needed for non-invertable matrix
        Ef=solve(temp1)%*%t(Wfkl)%*%Wskl%*%Es
        bkl[,pairno,traitno] <<- Ef
      }
      # preparation for apply functions
      temp1 <- freqGenotypes[episPairLoci[,1],]
      pk <- lapply(seq_len(nrow(temp1)), function(i) temp1[i,])
      temp1 <- freqGenotypes[episPairLoci[,2],]
      pl <- lapply(seq_len(nrow(temp1)), function(i) temp1[i,])
      # apply funtions for caculate bkl
      for (iTrait in 1:(nTrait)) {
        invisible(mapply(setupbkl,traitno=rep(iTrait,nEpisPairLoci),pairno=1:nEpisPairLoci,pk=pk,pl=pl))
      }
      
      # transfer effects to the right array	
      for (iTrait in 1:(nTrait)) {
        temp1=rep(0,nSnp)
        temp1[episPairLoci[,1]]=temp1[episPairLoci[,1]]+bkl[2,,iTrait]
        temp1[episPairLoci[,2]]=temp1[episPairLoci[,2]]+bkl[3,,iTrait]
        aQtleffects[,iTrait]=temp1[SNPinfo[,4]=='T']
        
        temp1=rep(0,nSnp)
        temp1[episPairLoci[,1]]=temp1[episPairLoci[,1]]+bkl[4,,iTrait]
        temp1[episPairLoci[,2]]=temp1[episPairLoci[,2]]+bkl[5,,iTrait]
        domValues[,iTrait]=temp1[SNPinfo[,4]=='T']
        
        episPairLociEffect[,iTrait]=bkl[6,,iTrait]
      }	    
    }
    
    #### Calculate functional effects from statistical effects
    if (infuncParameters %in% c('statisticalThree')) {
      
      # setup bkl using ivares-Castro and Calborg2007 method, to calculate functional genetic effects from statistical effects of a, d and e
      # only work if (nQtlEachQtlInteractwith==1)
      if (nQtlEachQtlInteractwith!=1) stop ('this method works with nQtlEachQtlInteractwith=1 only.') 
      
      bkl <- array(data = NA,dim=c(6,nEpisPairLoci,nTrait))
      
      setupbkl <- function(traitno,pairno,pk,pl) {
        #traitno=1;pairno=1;pk=freqGenotypes[episPairLoci[pairno,1],];pl=freqGenotypes[episPairLoci[pairno,2],]
        
        qtlno_k=sum(SNPinfo[1:episPairLoci[pairno,1],4]=='T')
        qtlno_l=sum(SNPinfo[1:episPairLoci[pairno,2],4]=='T')
        
        # setup Es manually like this or do it similarly as in setup Wsklmat or Wfkl for equivalent positions
        Es=matrix(0,nrow = 6,ncol=1)
        Es[,1]=c(1,aQtleffects[qtlno_k,traitno],aQtleffects[qtlno_l,traitno],domValues[qtlno_k,traitno],domValues[qtlno_l,traitno],episPairLociEffect[pairno,traitno])
        # Warning note: need to check if mean affect Ef
        
        Wskl <- setupWsklmat(pk,pl)
        Wfkl <- setupWfklmat2(tx1=c(0,1,2),p1_k=pk,p1_l=pl)
        
        temp1 <- t(Wfkl)%*%Wfkl
        if (det(temp1)<0.0000000001) diag(temp1)=diag(temp1)+0.00000001   # this might be needed for non-invertable matrix
        Ef=solve(temp1)%*%t(Wfkl)%*%Wskl%*%Es
        bkl[,pairno,traitno] <<- Ef
      }
      # preparation for apply functions
      temp1 <- freqGenotypes[episPairLoci[,1],]
      pk <- lapply(seq_len(nrow(temp1)), function(i) temp1[i,])
      temp1 <- freqGenotypes[episPairLoci[,2],]
      pl <- lapply(seq_len(nrow(temp1)), function(i) temp1[i,])
      # apply funtions for caculate bkl
      for (iTrait in 1:(nTrait)) {
        invisible(mapply(setupbkl,traitno=rep(iTrait,nEpisPairLoci),pairno=1:nEpisPairLoci,pk=pk,pl=pl))
      }
      
      # transfer effects to the right array	
      for (iTrait in 1:(nTrait)) {
        temp1=rep(0,nSnp)
        temp1[episPairLoci[,1]]=temp1[episPairLoci[,1]]+bkl[2,,iTrait]
        temp1[episPairLoci[,2]]=temp1[episPairLoci[,2]]+bkl[3,,iTrait]
        aQtleffects[,iTrait]=temp1[SNPinfo[,4]=='T']
        
        temp1=rep(0,nSnp)
        temp1[episPairLoci[,1]]=temp1[episPairLoci[,1]]+bkl[4,,iTrait]
        temp1[episPairLoci[,2]]=temp1[episPairLoci[,2]]+bkl[5,,iTrait]
        domValues[,iTrait]=temp1[SNPinfo[,4]=='T']
        
        episPairLociEffect[,iTrait]=bkl[6,,iTrait]
      }	    
    }
    
    # recalculate domdegreeValues given new rescaled domValues 
    for (iTrait in 1:nTrait) {
      temp1=abs(aQtleffects[1:nQtl,iTrait])
      temp1[temp1<0.00001]=0.00001   # avoid zero numbers lead to infinite values
      domdegreeValues[,iTrait] <- domValues[1:nQtl,iTrait]/temp1
    } 
    
    ######### add priors to startpop ### remember following loops regards orders of effects in startpop
    counter=0
    # additive
    for (iTrait in 1:(nTrait)) {
      startpop[idstartpop,(counter+1):(counter+nQtl)] <- aQtleffects[1:nQtl,iTrait]
      counter=counter+nQtl
    }
    # dominance
    for (iTrait in 1:(nTrait)) {
      startpop[idstartpop,(counter+1):(counter+nQtl)] <- domdegreeValues[1:nQtl,iTrait]
      counter=counter+nQtl
    }
    # epistasis
    for (iTrait in 1:(nTrait)) {
      startpop[idstartpop,(counter+1):(counter+nEpisPairLoci)] <- episPairLociEffect[1:nEpisPairLoci,iTrait]
      counter=counter+nEpisPairLoci
    }
    if (!counter==nLociforAlgorithm) {stop('Something wrong when assigning priors to startpop')}
    #note: The traits in two populations using the same biolocial allele effects. 
    ###### remember above loops regards orders of effects in startpop####################
  }
  # return
  return(startpop)
}

print('Prior of startpop loading...')
startpop <- generateStartPop(genodf=geno3way,nIndividuals=nstartpop,infuncParameters='biological',
                             covarTraits=aRes_covar,episcovarTraits=aaNRMSRes_covar,domcovarTraits=dNRMSRes_covar)   # aRes_covar assuming this is the variance of founder Pop1 (HWE).
addingpop <- generateStartPop(genodf=geno3way,nIndividuals=nstartpop,infuncParameters='statisticalTwo',
                              covarTraits=aRes_covar,episcovarTraits=aaNRMSRes_covar,domcovarTraits=dNRMSRes_covar)
startpop <- rbind(startpop,addingpop)

################################################
######## Set up functions to maximise ##########
########    accept solutions from GA  ##########
########        Improve speed         ##########
################################################

#### Required variables: startpop, aRes_covar, aNRMS_covar, aaNRMSRes_covar, dNRMSRes_covar, critWeight, critWeightcorAbio
#### Genotyped Pops needed: genoRes, genoNRMS and geno3way

#### uppervalues & lowervalues
nLociforAlgorithm=nTrait*nQtl +  nTrait*nQtl + nTrait*nEpisPairLoci
uppervalues=lowervalues=rep(NA,nTrait*nQtl +  nTrait*nQtl + nTrait*nEpisPairLoci)
counter=0
# additive
for (iTrait in 1:(nTrait)) {
  uppervalues[(counter+1):(counter+nQtl)] <- apply(startpop[,(counter+1):(counter+nQtl)], 2, max)
  lowervalues[(counter+1):(counter+nQtl)] <- apply(startpop[,(counter+1):(counter+nQtl)], 2, min)
  counter=counter+nQtl
}
# dominance
for (iTrait in 1:(nTrait)) {
  uppervalues[(counter+1):(counter+nQtl)] <- apply(startpop[,(counter+1):(counter+nQtl)], 2, quantile,probs=0.95)
  lowervalues[(counter+1):(counter+nQtl)] <- apply(startpop[,(counter+1):(counter+nQtl)], 2, quantile,probs=0.05)
  counter=counter+nQtl
}
# epistasis
for (iTrait in 1:(nTrait)) {
  uppervalues[(counter+1):(counter+nEpisPairLoci)] <- apply(startpop[,(counter+1):(counter+nEpisPairLoci)], 2, max)
  lowervalues[(counter+1):(counter+nEpisPairLoci)] <- apply(startpop[,(counter+1):(counter+nEpisPairLoci)], 2, min)
  counter=counter+nEpisPairLoci
}

# remove outliners for starting pop
for (i in 1:nLociforAlgorithm) {
  temp1 <- uppervalues[i] < max(startpop[,i])
  if (temp1) startpop[startpop[,i]>uppervalues[i],i]=uppervalues[i]
  temp1 <- lowervalues[i] > min(startpop[,i])
  if (temp1) startpop[startpop[,i]<lowervalues[i],i]=lowervalues[i]
}

############# Take these variables outside of function criteriaAccept ###########
### this would improve significant speed

###### calculate allele frequency
MAF <- function(x){sum(x,na.rm=T)/(sum(!is.na(x))*2)}
GenoFreq <- function(x,y){sum(x==y)/length(x)}
# Frequency based on genotypefiles
for (ifile in lstgenotypefiles) {
  genoTemp = get(ifile)        # frequencyBasedPop
  # alleles
  freqAlleles <- array(data = 0,dim=c(nSnp,2))
  #MAF <- function(x){sum(x,na.rm=T)/(sum(!is.na(x))*2)}
  freqAlleles[,2]=apply(genoTemp, 2, MAF)   # q
  freqAlleles[,1]=1-freqAlleles[,2]         # p
  # genotypes
  freqGenotypes <- array(data = 0,dim=c(nSnp,3))
  #GenoFreq <- function(x,y){sum(x==y)/length(x)}
  counter=0
  for (iAllele in 1:2) {
    for (jAllele in iAllele:2) {
      counter=counter+1
      freqGenotypes[,counter] <- apply(genoTemp, 2, GenoFreq,y=(counter-1))
    }
  }
  assign(paste0('freqAlleles_',ifile), freqAlleles)
  assign(paste0('freqGenotypes_',ifile), freqGenotypes)
  genoTemp <- freqAlleles <- freqGenotypes <- NULL
}

########## 

# temp variables needed
nAnimal=dim(genoFounderRes)[1]

# setup hai matrices for a1 par1 pops contributed to plot 
tempfreqGeno <- freqGenotypes_genoFounderRes[SNPinfo[,4]=='T',]
hai_a1 <- temp1 <- genoFounderRes[,SNPinfo[,4]=='T']
tempfunc <- function(iQtl) {
  for (i1 in c(0,1,2)) {
    hai_a1[temp1[,iQtl]==i1,iQtl] <<- findhai(px=tempfreqGeno[iQtl,],allelecount=i1)  }
}
invisible(mapply(tempfunc,iQtl=1:ncol(hai_a1)))    # using apply
#hai_a1=hai_a1/2								   # devided by 2 because half of genes transfered. Not devided by two because here aiming to get variance of HWE population

# setup hai matrices for a2 par2 pops contributed to plot 
tempfreqGeno <- freqGenotypes_genoFounderNRMS[SNPinfo[,4]=='T',]
hai_a2 <- temp1 <- genoFounderNRMS[,SNPinfo[,4]=='T']
tempfunc <- function(iQtl) {
  for (i1 in c(0,1,2)) {
    hai_a2[temp1[,iQtl]==i1,iQtl] <<- findhai(px=tempfreqGeno[iQtl,],allelecount=i1)  }
}
invisible(mapply(tempfunc,iQtl=1:ncol(hai_a2)))    # using apply
#hai_a2=hai_a2/2								   # devided by 2 because half of genes transfered. Not devided by two because here aiming to get variance of HWE population

# setup hdi matrices
tempallelfreqpar1 <- freqAlleles_genoFounderRes[SNPinfo[,4]=='T',1]
tempallelfreqpar2 <- freqAlleles_genoFounderNRMS[SNPinfo[,4]=='T',1]
temppar1 <- genoinbredRes[,SNPinfo[,4]=='T']
temppar2 <- geno2wayNRMS[,SNPinfo[,4]=='T']
hdi3way <- geno3way[,SNPinfo[,4]=='T']
tempfunc <- function(iQtl) {
  for (i1 in c(0,1,2)) {
    for (i2 in c(0,1,2)) {
      hdi3way[temppar1[,iQtl]==i1 & temppar2[,iQtl]==i2,iQtl] <<- 
        findhdiPlot(p1=tempallelfreqpar1[iQtl],p2=tempallelfreqpar2[iQtl],allelecount1=i1,allelecount2=i2)
    }
  }
}
invisible(mapply(tempfunc,iQtl=1:ncol(hdi3way)))    # using apply

# setup haihai matrices for geno3way
haihai3way <- array(data = 0,dim=c(nAnimal,nEpisPairLoci))
tempallelfreqpar1 <- freqAlleles_genoFounderRes[,1]
tempallelfreqpar2 <- freqAlleles_genoFounderNRMS[,1]
temppar1 <- genoinbredRes[,]
temppar2 <- geno2wayNRMS[,]
tempfunc <- function(idd) {
  haihaiVect <- mapply(findhaihaiPlot,
                       ap1_k=tempallelfreqpar1[episPairLoci[,1]],ap2_k=tempallelfreqpar2[episPairLoci[,1]],ac1_k=temppar1[idd,episPairLoci[,1]],ac2_k=temppar2[idd,episPairLoci[,1]],
                       ap1_l=tempallelfreqpar1[episPairLoci[,2]],ap2_l=tempallelfreqpar2[episPairLoci[,2]],ac1_l=temppar1[idd,episPairLoci[,2]],ac2_l=temppar2[idd,episPairLoci[,2]])
  haihai3way[idd,] <<- haihaiVect
}
invisible(mapply(tempfunc,idd=1:nAnimal))    # using apply instead of loop for

# setup bklOther
bklOther <- array(data = NA,dim=c(8,8,nEpisPairLoci))
tempallelfreqpar1 <- freqAlleles_genoFounderRes[,1]
tempallelfreqpar2 <- freqAlleles_genoFounderNRMS[,1]
Wfkl=setupWfklmat1(tx1=c(0,1,2),tx2=c(0,1,2))
tempfunc <- function(pairno) {
  #pairno=1
  Wskl <- setupWsklPlotmat(tx1=c(0,1,2),tx2=c(0,1,2),
                           ap1_k=tempallelfreqpar1[episPairLoci[pairno,1]],ap2_k=tempallelfreqpar2[episPairLoci[pairno,1]],
                           ap1_l=tempallelfreqpar1[episPairLoci[pairno,2]],ap2_l=tempallelfreqpar2[episPairLoci[pairno,2]])
  temp1 <- t(Wskl)%*%Wskl
  if (det(temp1)<0.0000000001) diag(temp1)=diag(temp1)+0.00000001   # this might be needed for non-invertable matrix
  bklOther[,,pairno] <<- solve(temp1)%*%t(Wskl)%*%Wfkl
}
# apply funtions for caculate bklOther
invisible(mapply(tempfunc,pairno=1:nEpisPairLoci))

# other variables needed. Put outside function to improve speed
# diag(aRes_covar)/2 because this is the contribution of Res to plot phenotype. assuming aRes_covar estimated using Gmatrix /2 due to inbred population.
# (diag(aNRMS_covar)/2)*4 because this is the contribution of NR to plot phenotype. assuming aNRMS_covar estimated using Gmatrix /2 due to inbred population.
varvectExpected <- c(diag(aRes_covar),diag(aNRMS_covar),diag(aaNRMSRes_covar),diag(dNRMSRes_covar))  
critWeightcorAbio <- abs(critWeightcorAbio)
critWeight <- abs(critWeight)

# now set up functions to maximise # using with global variables 
criteriaAccept <- function(tempSol) {  #tempSol=startpop[1,];tempSol=opt@solution
  # read tempSol to arrays 
  # biolocial qtl effects
  if (length(tempSol)>nLociforAlgorithm) {
    print('Something is wrong.')
    stop('Dimension is not correct.')
  }
  # biological additive
  aQtleffects <- array(data = 0,dim=c(nQtl,nTrait))
  counter=0
  for (iTrait in 1:nTrait) {
    aQtleffects[1:nQtl,iTrait] <- tempSol[(counter+1):(counter+nQtl)]
    counter=counter+nQtl
  }
  # biological dominance degree
  domdegreeValues <- array(data = 0,dim=c(nQtl,nTrait))
  for (iTrait in 1:nTrait) {
    domdegreeValues[1:nQtl,iTrait] <- tempSol[(counter+1):(counter+nQtl)]
    counter=counter+nQtl
  }
  # biological epistasis
  episPairLociEffect <- array(data = 0,dim=c(nEpisPairLoci,nTrait))	
  for (iTrait in 1:nTrait) {
    episPairLociEffect[1:nEpisPairLoci,iTrait] <- tempSol[(counter+1):(counter+nEpisPairLoci)]
    counter=counter+nEpisPairLoci
  }
  
  # caculate biological dominance from domdegree
  domValues <- array(data = 0,dim=c(nQtl,nTrait))
  for (iTrait in 1:nTrait) {
    domValues[1:nQtl,iTrait] <- abs(aQtleffects[1:nQtl,iTrait])*domdegreeValues[1:nQtl,iTrait]
  }
  domdegreeValues<-NULL
  
  #### Calculate statistical effects given functional effects simulated
  # setup bkl using ivares-Castro and Calborg2007 method, to calculate functional genetic effects from statistical effects of a, d and e
  # only work if (nQtlEachQtlInteractwith==1) # if (nQtlEachQtlInteractwith!=1) stop ('this method works with nQtlEachQtlInteractwith=1 only.') 
  bkl <- array(data = NA,dim=c(8,nEpisPairLoci,nTrait))
  Ef=matrix(0,nrow = 8,ncol=1)
  setupbkl <- function(traitno,pairno) {
    #traitno=1;pairno=1
    qtlno_k=sum(SNPinfo[1:episPairLoci[pairno,1],4]=='T')
    qtlno_l=sum(SNPinfo[1:episPairLoci[pairno,2],4]=='T')
    # setup Ef manually like this or do it similarly as in setup Wsklmat or Wfkl for equivalent positions
    Ef[,1]=c(1,aQtleffects[qtlno_k,traitno],aQtleffects[qtlno_l,traitno],aQtleffects[qtlno_k,traitno],aQtleffects[qtlno_l,traitno],
             domValues[qtlno_k,traitno],domValues[qtlno_l,traitno],episPairLociEffect[pairno,traitno])
    bkl[,pairno,traitno] <<- bklOther[,,pairno] %*% Ef
  }
  # apply funtions for caculate bkl
  for (iTrait in 1:(nTrait)) {
    invisible(mapply(setupbkl,traitno=rep(iTrait,nEpisPairLoci),pairno=1:nEpisPairLoci))
  }
  
  # transfer effects to the right arrays
  a1Qtleffects_Stat <- a2Qtleffects_Stat <- domValues_Stat <- array(data = 0,dim=c(nQtl,nTrait))
  episPairLociEffect_Stat <- array(data = 0,dim=c(nEpisPairLoci,nTrait))
  
  for (iTrait in 1:(nTrait)) {
    temp1=rep(0,nSnp)
    temp1[episPairLoci[,1]]=temp1[episPairLoci[,1]]+bkl[2,,iTrait]
    temp1[episPairLoci[,2]]=temp1[episPairLoci[,2]]+bkl[3,,iTrait]
    a1Qtleffects_Stat[,iTrait]=temp1[SNPinfo[,4]=='T']
    
    temp1=rep(0,nSnp)
    temp1[episPairLoci[,1]]=temp1[episPairLoci[,1]]+bkl[4,,iTrait]
    temp1[episPairLoci[,2]]=temp1[episPairLoci[,2]]+bkl[5,,iTrait]
    a2Qtleffects_Stat[,iTrait]=temp1[SNPinfo[,4]=='T']
    
    temp1=rep(0,nSnp)
    temp1[episPairLoci[,1]]=temp1[episPairLoci[,1]]+bkl[6,,iTrait]
    temp1[episPairLoci[,2]]=temp1[episPairLoci[,2]]+bkl[7,,iTrait]
    domValues_Stat[,iTrait]=temp1[SNPinfo[,4]=='T']
    
    episPairLociEffect_Stat[,iTrait]=bkl[8,,iTrait]
  }
  
  ### calculate statistical variances of par1 pops contributed to plot 
  idgv <- array(data = 0,dim=c(nAnimal,nTrait)) 
  for (iTrait in 1:(nTrait)) {
    idgv[,iTrait] <- hai_a1 %*% a1Qtleffects_Stat[,iTrait]
  }
  v_a1Stat=var(idgv)   # calculate covariances for traits
  
  ### calculate statistical variances of par1 pops contributed to plot 
  idgv <- array(data = 0,dim=c(nAnimal,nTrait)) 
  for (iTrait in 1:(nTrait)) {
    idgv[,iTrait] <- hai_a2 %*% a2Qtleffects_Stat[,iTrait]
  }
  v_a2Stat=var(idgv)   # calculate covariances for traits
  
  #### calculate statistical epistasis by animals for geno3way
  idgv <- array(data = 0,dim=c(nAnimal,nTrait))
  for (iTrait in 1:(nTrait)) {
    idgv[,iTrait] <- haihai3way %*% episPairLociEffect_Stat[,iTrait]            # episPairLociEffect = episPairLociEffect_Stat
  }
  v_aaStat=var(idgv)   # calculate covariances for traits
  
  #### calculate statistical additive genetic effects by animals for geno3way
  idgv <- array(data = 0,dim=c(nAnimal,nTrait)) 
  for (iTrait in 1:(nTrait)) {
    idgv[,iTrait] <- hdi3way %*% domValues_Stat[,iTrait]                       # domValues_Stat = domValues
  }
  v_dStat=var(idgv)   # calculate covariances for traits
  
  #criteria for acceptance of solutions (convergence)
  crit =0
  # correlation between traits  # make sure multiple traits are correlated.
  crit=crit+sum(abs(round(cov2cor(aRes_covar)-cov2cor(v_a1Stat),digits=1))*critWeightcorAbio)
  crit=crit+sum(abs(round(cov2cor(aNRMS_covar)-cov2cor(v_a2Stat),digits=1))*critWeightcorAbio)
  crit=crit+sum(abs(round(cov2cor(aaNRMSRes_covar)-cov2cor(v_aaStat),digits=1))*critWeightcorAbio)
  crit=crit+sum(abs(round(cov2cor(dNRMSRes_covar)-cov2cor(v_dStat),digits=1))*critWeightcorAbio)
  
  # variance
  varvectObserved <- c(diag(v_a1Stat),diag(v_a2Stat),diag(v_aaStat),diag(v_dStat))
  
  crit=crit+sum(abs(round(varvectExpected-varvectObserved,digits=1))*critWeight)
  
  # maximize crit
  crit = -1*round(crit,digits=0)
  #print(crit) #
  
  return(crit)
}

# run algorthm to optimize solutions
print('Optimization running...')
opt <- ga(type = "real-valued", 
          fitness = criteriaAccept,
          lower = lowervalues,
          upper = uppervalues,
          popSize = dim(startpop)[1],
          maxiter = maxiterOpt,
          suggestions = startpop,
            maxFitness=0 )
# more running
if (opt@fitnessValue < -2) {
  for (itemp in 1:5) {
    if (opt@fitnessValue > -2) next
    print(paste('Optimization rerun:',itemp,'th'))
    rerunpop <- opt@population
    opt <- ga(type = "real-valued", 
              fitness = criteriaAccept,
              lower = lowervalues,
              upper = uppervalues,
              popSize = dim(rerunpop)[1],
              maxiter = maxiterOpt,
              suggestions = rerunpop,
              maxFitness=0 )
  }
}
if (opt@fitnessValue < -2) print(paste('Optimization not successful for rep:',iRepScheme))

print(summary(opt)); plot(opt)
print('Optimization completed.')

#####################################################################
######## Check optimum solutions & keep as a,d and aa effects #######
#####################################################################
print('******** Solution reports:')
# check solutions
tempSol=opt@solution

if (!is.null(dim(tempSol))) {
  print(paste('Dim of Sol:',paste(dim(tempSol),collapse=' ')))
  #write(tempSol, file = 'tempSol1.txt')
  if (dim(tempSol)[1]>1) tempSol=tempSol[1,]
  write(tempSol, file = 'tempSol2.txt')
}

# biological additive
aQtleffects <- array(data = 0,dim=c(nQtl,nTrait))
counter=0
for (iTrait in 1:nTrait) {
  aQtleffects[1:nQtl,iTrait] <- tempSol[(counter+1):(counter+nQtl)]
  counter=counter+nQtl
}
# biological dominance degree
domdegreeValues <- array(data = 0,dim=c(nQtl,nTrait))
for (iTrait in 1:nTrait) {
  domdegreeValues[1:nQtl,iTrait] <- tempSol[(counter+1):(counter+nQtl)]
  counter=counter+nQtl
}
# biological epistasis
episPairLociEffect <- array(data = 0,dim=c(nEpisPairLoci,nTrait))	
for (iTrait in 1:nTrait) {
  episPairLociEffect[1:nEpisPairLoci,iTrait] <- tempSol[(counter+1):(counter+nEpisPairLoci)]
  counter=counter+nEpisPairLoci
}
if ((!counter==nLociforAlgorithm)|(!counter==length(tempSol))) {stop('Something wrong with solution from opt.')}

# caculate biological dominance from domdegree
domValues <- array(data = 0,dim=c(nQtl,nTrait))
for (iTrait in 1:nTrait) {
  domValues[1:nQtl,iTrait] <- abs(aQtleffects[1:nQtl,iTrait])*domdegreeValues[1:nQtl,iTrait]
}
# Report distribution of solutions for biological effects
for (iTrait in 1:(nTrait)) {
  hist(aQtleffects[,iTrait],main=paste('Histogram of biological aQtleffects for Trait',iTrait))
  hist(episPairLociEffect[,iTrait],main=paste('Histogram of biological episPairLociEffect for Trait',iTrait))
  hist(domValues[,iTrait],main=paste('Histogram of biological domValues for Trait',iTrait))
  hist(domdegreeValues[,iTrait],main=paste('Histogram of biological domdegreeValues for Trait',iTrait))
}	
domdegreeValues <- c() # not use anymore, remove for computation efficiency.
if (!exists('iRepScheme')) iRepScheme=1

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
  write.table(temp0,file=filename,sep = "\t",row.names=F,col.names=!file.exists(filename),quote = FALSE,append=file.exists(filename))
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
  write.table(temp0,file=filename,sep = "\t",row.names=F,col.names=!file.exists(filename),quote = FALSE,append=file.exists(filename))
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
  write.table(temp0,file=filename,sep = "\t",row.names=F,col.names=!file.exists(filename),quote = FALSE,append=file.exists(filename))
  counter=counter+1
}
if (counter==0) print('*** No genotype files found for reporting variances. ***')

#### plot
if (exists('genoinbredRes') & exists('geno2wayNRMS')) { 
  trueVarCalPlot <- calculateStatEffectPlotOne(genoPar1=genoinbredRes,genoPar2=geno2wayNRMS,returnvalues='variances')
  trueVarCalPlot <- data.frame(iRep=iRepScheme,trueVarCalPlot)
  filename='varSchemTrueAvCPlotOne.res'
  write.table(trueVarCalPlot,file=filename,sep = "\t",row.names=F,col.names=!file.exists(filename),quote = FALSE,append=file.exists(filename))
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
