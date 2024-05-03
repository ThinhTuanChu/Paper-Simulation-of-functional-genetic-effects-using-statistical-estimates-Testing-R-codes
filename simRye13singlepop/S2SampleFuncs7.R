# sampling offsprings given two parents 
sampleOffspring <- function(idd=NULL,par1,par2,popid=NULL,fixedFactors=NULL,randomFactors=NULL) {
  # idd=NULL;par1=10;par2=13;popid=NULL;fixedFactors=NULL;randomFactors=NULL
  
  # get genotype of parents 
  genoPar1=genoPar2=matrix(NA, nrow = 2, ncol= length(SNPinfo[,1])) 
  genoPar1[1,] <- haplo_pat[par1,]
  genoPar1[2,] <- haplo_mat[par1,]
  genoPar2[1,] <- haplo_pat[par2,]
  genoPar2[2,] <- haplo_mat[par2,]
  
  # id of offspring
  if (is.null(idd)) {
    maxid=maxid+1
    assign("maxid", maxid, envir = .GlobalEnv)
    idd=maxid
  }
  # check input
  if (par1>maxid | par2>maxid | idd>(maxid+1)) stop('wrong input for this function. par1,par2,idd must be < maxid.')
  
  # sampling genomic genotypes of offsprings given genotypes of two parents 
  lastLocus=0
  temp_pat=temp_mat=rep(NA,length(SNPinfo[,1]))
  for (ichrom in 1:nChr) { #ichrom=1
    
    firstLocus=lastLocus+1
    lastLocus=firstLocus+chrom[[ichrom]]$nSnpEachChr-1
    if (chrom[[ichrom]]$nSnpEachChr==0) next
    for (iparent in 1:2) { #iparent=1
      
      nRecomb=rpois(n=1, lambda=chrom[[ichrom]]$chromlength/100.0)
      if (nRecomb>20) nRecomb=20
      if (nRecomb<0) nRecomb=0
      
      ipchr0=rbinom(n=1, size=1, 0.5)   # check this line
      
      posRecomb=runif(n=nRecomb, min = 0.0, max = chrom[[ichrom]]$chromlength)
      posRecomb=posRecomb[order(posRecomb)]
      if (nRecomb==0) {
        
        ipchr=ipchr0+1
        
        if (iparent==1) temp_pat[firstLocus:lastLocus]=genoPar1[ipchr,firstLocus:lastLocus]
        
        if (iparent==2) temp_mat[firstLocus:lastLocus]=genoPar2[ipchr,firstLocus:lastLocus]
      } else 	{
        
        posinfo=SNPinfo[SNPinfo[,2]==ichrom,3]
        frompos=0; fromloc=firstLocus
        
        for (il in 1:(nRecomb+1)) { #il=1
          
          if (il<(nRecomb+1)) {
            topos=posRecomb[il]
          } else { 			   topos=chrom[[ichrom]]$chromlength+1  			}
          
          ntoloc=sum(posinfo>=frompos & posinfo<topos)
          
          if (ntoloc==0) { 
            fromloc=fromloc+ntoloc
            frompos=topos
            next
          }
          
          ipchr=(((il-1)+ipchr0)%%2)+1
          
          if (iparent==1) temp_pat[fromloc:(fromloc+ntoloc-1)]=genoPar1[ipchr,fromloc:(fromloc+ntoloc-1)]
          if (iparent==2) temp_mat[fromloc:(fromloc+ntoloc-1)]=genoPar2[ipchr,fromloc:(fromloc+ntoloc-1)]
          
          fromloc=fromloc+ntoloc
          frompos=topos
        }
        
        
        # for (iloc in firstLocus:lastLocus) {
        # 
        #   currentRecomb=sum(posRecomb[1:nRecomb]<SNPinfo[iloc,3])  # count true
        #   ipchr=((currentRecomb+ipchr0)%%2)+1
        # 
        #   if (iparent==1) temp_pat[iloc]=genoPar1[ipchr,iloc]
        # 
        #   if (iparent==2) temp_mat[iloc]=genoPar2[ipchr,iloc]
        # }  
        
      }
    }
  }
  #temp1=data.frame(temp_pat,temp_mat,genoPar1[1,],genoPar1[2,],genoPar2[1,],genoPar2[2,])
  
  if (anyNA(temp_pat) | anyNA(temp_mat)) stop('something is wrong. temp_pat or temp_mat is not sampled correctly.')
  
  # add genotypes of the idd to haplo_pat & haplo_mat
  temp_pat=matrix(temp_pat,nrow=1)
  temp_mat=matrix(temp_mat,nrow=1)
  
  haplo_pat1 <- haplo_pat
  haplo_mat1 <- haplo_mat
  temp1=dim(haplo_pat1)[1]; temp2=dim(haplo_mat1)[1]
  if (temp1!=temp2) stop('something is wrong. haplo_pat1 and haplo_mat1 has unequal sizes.')
  
  if (idd>temp1) {
    haplo_pat1=rbind(haplo_pat1,temp_pat)
    haplo_mat1=rbind(haplo_mat1,temp_mat)
  } else {
    haplo_pat1[idd,]=temp_pat
    haplo_mat1[idd,]=temp_mat
  }
  assign("haplo_pat", haplo_pat1, envir = .GlobalEnv)
  assign("haplo_mat", haplo_mat1, envir = .GlobalEnv)
  
  # sample infor genetic values, observations of the offspring
  if (par1==par2) {
    par1temp=pop$par1[par1]
    par2temp=pop$par2[par2]
  } else {
    par1temp=par1
    par2temp=par2
  }
  temppop <- data.frame(id=idd,par1=par1temp,par2=par2temp,birthpop=0,currentpop=0)
  for (iTrait in 1:(nTrait)) {
    temppop[,paste0('trait_av',iTrait)]=0   # dont need these
    temppop[,paste0('trait_dv',iTrait)]=0   # dont need these
    temppop[,paste0('trait_tgv',iTrait)]=0	
    temppop[,paste0('trait_episgv',iTrait)]=0	
  }
  for (iTrait in 1:(nObs)) {
    temppop[,paste0('trait_obs',iTrait)]=0
  }
  for (iTrait in 1:(nEbv)) {
    temppop[,paste0('trait_polyEbv',iTrait)]=0
    temppop[,paste0('trait_genoEbv',iTrait)]=0
  }
  for (iTrait in 1:(nRes)) {
    temppop[,paste0('trait_residual',iTrait)]=0
  }
  if (is.null(popid)) {
    popid=pop$birthpop[par2]
  } 
  temppop$birthpop=popid
  temppop$currentpop=popid
  
  #add other factors
  for (iFactor in 1:nfixedFactors) {
    temppop[,paste0('fixedFactor',iFactor)]=ifelse(is.null(fixedFactors)|length(fixedFactors)!=nfixedFactors,0,fixedFactors[iFactor])
  }
  for (iFactor in 1:nranFactors) {
    temppop[,paste0('randomFactor',iFactor)]=ifelse(is.null(randomFactors)|length(randomFactors)!=nranFactors,-999.0,randomFactors[iFactor])
  }
  
  # true genetic values of aniaml id
  temp1 <- temp_pat[1,SNPinfo[,4]=='T'] + temp_mat[1,SNPinfo[,4]=='T'] - 3
  for (iTrait in 1:(nTrait)) {
    temppop[1,paste0('trait_av',iTrait)]=temppop[1,paste0('trait_av',iTrait)]+
      sum(aQtleffects[,iTrait] * temp1)
    temppop[1,paste0('trait_dv',iTrait)]=temppop[1,paste0('trait_dv',iTrait)]+
      sum(domValues[temp1==0,iTrait])
  }
  
  # true epistasis genetic values of aniaml id
  if (!is.null(episPairLociEffect)) {
    
    temp1 <- temp_pat[1,episPairLoci[,1]]
    temp2 <- temp_mat[1,episPairLoci[,1]]
    xAk <- temp1 + temp2 - 3                  # scaled genotype dosages at locus k
    temp1 <- temp_pat[1,episPairLoci[,2]]
    temp2 <- temp_mat[1,episPairLoci[,2]]
    xAl <- temp1 + temp2 - 3                  # scaled genotype dosages at locus l
    for (iTrait in 1:(nTrait)) {
      temppop[1,paste0('trait_episgv',iTrait)]=temppop[1,paste0('trait_episgv',iTrait)]+
        sum(episPairLociEffect[1:nEpisPairLoci,iTrait]*xAk*xAl)
    }
  }
  
  # sample residuals effects of aniaml id
  temppop[1,paste0('trait_residual',c(1:nRes))]=mvrnorm(1,rep(0,nRes),varResiduals)
  
  # trait_tgv=trait_av+trait_dv+trait_episgv
  for (iTrait in 1:(nTrait)) {  
    temppop[1,paste0('trait_tgv',iTrait)]=temppop[1,paste0('trait_av',iTrait)]+
      temppop[1,paste0('trait_dv',iTrait)]+temppop[1,paste0('trait_episgv',iTrait)]
  }
  # sample obs. for now, obs is simply=trait_tgv+trait_residual
  for (iTrait in 1:(nObs)) {  
    temppop[1,paste0('trait_obs',iTrait)]=temppop[1,paste0('trait_tgv',iTrait)]+
      temppop[1,paste0('trait_residual',iTrait)]
  }
  
  # add results to global enviroment pop
  pop1 <- pop
  if (idd > dim(pop)[1]) {
    pop1 <- rbind(pop1,temppop)
  } else {
    pop1[idd,] <- temppop
  }
  assign("pop", pop1, envir = .GlobalEnv) 
  
}

# Functions to generate data for BLUP prediction
generatePhenoData <- function(idtobeused=NULL,writedata=FALSE,returndata=TRUE) {
  
  ####### phenotype files
  # animals used
  if (is.null(idtobeused)) {
    idtobeused = pop$id
  }
  idtobeused <- idtobeused[order(idtobeused)]
  
  temppop <- pop[pop$id %in% idtobeused,]
  
  if (length(idtobeused)<2 | length(temppop[,1])<2) {     stop(paste('Error: no id matches criteria.'))   }  
  
  colnamesused <- c('id','par1','par2','birthpop','currentpop')
  
  #add other factors
  for (i in 1:nfixedFactors) {
    colnamesused <- c(colnamesused,paste0('fixedFactor',i))
  }
  for (i in 1:nObs) {
    colnamesused <- c(colnamesused,paste0('trait_obs',i))
  }
  
  for (i in 1:nranFactors) {
    colnamesused <- c(colnamesused,paste0('randomFactor',i))
  }
  temppop <- temppop[,colnamesused]
  
  if (writedata) write.table(temppop,file ='dmudat',col.names = F,row.names = F,quote = FALSE,na="-9999",sep =" ")
  
  if (returndata) return(temppop)
}
######### Genotype files
generateGenoData <- function(idtobeused=NULL,genoDataFormat=NULL,locitype='all',writedata=FALSE,returndata=TRUE) {  
  #options for locitype: 'all', 'qtl','markers'
  # animals used
  if (is.null(idtobeused)) {
    idtobeused = pop$id
  }  
  if (length(idtobeused)<2) {     stop(paste('Error: no id matches criteria.'))   }
  idtobeused <- unique(idtobeused)
  if (any(idtobeused>maxid)) {
    idtobeused <- idtobeused[idtobeused<=maxid]
    warning('Warning: some ID dont exist, thus are removed.')
  }
  idtobeused <- idtobeused[order(idtobeused)]
  
  # extract markers and genotype
  if (! locitype %in% c('all', 'qtl','markers'))  {
    locitype='all'
    print("Invalid locitype given. Locitype=all is used")	}
  if (locitype=='qtl') {
    temp_pat <- haplo_pat[idtobeused,SNPinfo[,4]=="T"]
    temp_mat <- haplo_mat[idtobeused,SNPinfo[,4]=="T"] 
    nLociUsed <- nQtl
  }
  if (locitype=='markers') {
    temp_pat <- haplo_pat[idtobeused,SNPinfo[,4]=="F"]
    temp_mat <- haplo_mat[idtobeused,SNPinfo[,4]=="F"]  
    nLociUsed <- nMarkers
  }
  if (locitype=='all') {
    temp_pat <- haplo_pat[idtobeused,]
    temp_mat <- haplo_mat[idtobeused,]  
    nLociUsed <- nSnp
  }
  
  genodf=NULL
  
  if (!is.null(genoDataFormat)) {
    if (genoDataFormat==2) {
      if (nAlleles>2) { stop(paste('Error: genoDataFormat',genoDataFormat,'does not work with nAlleles of',nAlleles))}
      genodf=matrix(NA, nrow=length(idtobeused),ncol=nLociUsed+1)
      genodf[,1]=idtobeused
      for (id in 1:length(idtobeused)) {
        genodf[id,2:(nLociUsed+1)]=(temp_pat[id,]-1)+(temp_mat[id,]-1)
      }
    }
    if (genoDataFormat==3) {
      genodf=matrix(NA, nrow=length(idtobeused),ncol=nLociUsed+1)
      genodf[,1]=idtobeused
      for (ipop in 1:(nFounderPop-1)) {
        temp_pat[temp_pat==(2*ipop+1)]=1
        temp_pat[temp_pat==(2*ipop+2)]=2
        temp_mat[temp_mat==(2*ipop+1)]=1
        temp_mat[temp_mat==(2*ipop+2)]=2		
      }
      for (id in 1:length(idtobeused)) {
        genodf[id,2:(nLociUsed+1)]=(temp_pat[id,]-1)+(temp_mat[id,]-1)
      }
    }	
  }  
  
  if (is.null(genoDataFormat) | is.null(genodf)) {
    genodf=matrix(NA, nrow=length(idtobeused),ncol=nLociUsed*2+1)
    genodf[,1]=idtobeused
    for (id in 1:length(idtobeused)) {
      genodf[id,seq(2,nLociUsed*2+1,by=2)]=temp_pat[id,]
      genodf[id,seq(3,nLociUsed*2+1,by=2)]=temp_mat[id,]
    }
  }
  
  if (writedata) write.table(genodf,file ='marker.RawData',col.names = F,row.names = F,quote = FALSE,na="9",sep =" ")  
  
  if (returndata) return(genodf)
}

####### report genetic phenotypic trend
reportGeneticPhenotypicTrend <- function(rep_nr=NULL,timestep=NULL) {  
  
  if (is.null(rep_nr)) rep_nr='all'
  if (is.null(timestep)) timestep=unique(pop$fixedFactor1)
  
  popid <- unique(pop$birthpop)
  poptemp <- data.frame()
  
  for (itime in timestep) {
    for (ipop in popid) {
      
      temp1 <- pop[pop$birthpop==ipop & pop$fixedFactor1==itime,]
      
      poptemp1 <- data.frame(rep_nr,timestep=itime,pop=ipop)	  
      for (i in 1:nTrait) {
        poptemp1[,paste0('trait_tgv',i)]=mean(temp1[,paste0('trait_tgv',i)])	
      }
      for (i in 1:nObs) {
        poptemp1[,paste0('trait_obs',i)]=mean(temp1[,paste0('trait_obs',i)])	
      }
      poptemp=rbind(poptemp,poptemp1)
    }
  }
  
  filename <- paste0('geneticphenotypictrend.res')
  write.table(poptemp,file =filename,col.names =!file.exists(filename),row.names = F,quote = FALSE,append=file.exists(filename))
  
}

######################################################################################################
####### Functions to calculate statistical effects of additive genetics, dominance & epistasis #######
### given we have dataset haplo_mat, haplo_pat and arrays of biological effects (a,d,epis)############
################ domValues; aQtleffects; episPairLoci; episPairLociEffect ############################
calculateStatEffectOne <- function(genodf=NULL,idtobeused=NULL,idfrequencyCalculation=NULL,returnvalues='variances') {   
  # genodf=NULL;idtobeused=pop$id[pop$currentpop==1];idfrequencyCalculation=NULL;rounddigits=3; returnvalues='variances'
  # check inputs
  if (!returnvalues %in% c('variances','effects','both')) returnvalues='variances'
  
  if (!is.null(genodf)) {  # if genotype file given, check correct formats
    genodf <- as.matrix(genodf)
  } else {
    genodf = haplo_mat + haplo_pat - 2
    genodf = cbind(matrix(c(1:length(genodf[,1])),nrow=length(genodf[,1])),genodf)  # add column ID
  }
  if (is.null(idtobeused))             {  idtobeused = genodf[,1]           }
  if (is.null(idfrequencyCalculation)) {  idfrequencyCalculation=idtobeused }
  if ((length(genodf[,1]))<length(unique(c(idtobeused,idfrequencyCalculation)))) stop(paste('Dimension of genodf is not right. Nrow must be at least:',length(unique(c(idtobeused,idfrequencyCalculation)))))
  if ((length(genodf[1,])-1)!=nSnp) stop(paste('Dimension of genodf is not right. Ncolumn must be:',nSnp+1))
  if (any(!genodf[,2:length(genodf[1,])] %in% c(0,1,2))) stop('Genotype codes are not in correct format. It must be 0,1,2.')
  if (length(idtobeused)!=length(unique(idtobeused))) {
    idtobeused=unique(idtobeused); 	  print('Repeated ID are removed.')
  }   
  # keep needed genotypes only
  allID <- unique(c(idtobeused,idfrequencyCalculation))
  allID <- allID[order(allID)]   
  genodf <- genodf[genodf[,1] %in% allID,]
  
  # set rownames of geno as id 
  geno <- genodf[,-1]
  rownames(geno) <- genodf[,1] ; genodf <- NULL
  idtobeused=as.character(idtobeused)
  idfrequencyCalculation=as.character(idfrequencyCalculation)
  
  ### calculate frequency of alleles and genotypes
  genotemp <- geno[idfrequencyCalculation,]
  # alleles
  freqAlleles <- array(data = 0,dim=c(nSnp,2))
  MAF <- function(x){sum(x,na.rm=T)/(sum(!is.na(x))*2)}
  freqAlleles[,2]=apply(genotemp, 2, MAF)   # q
  freqAlleles[,1]=1-freqAlleles[,2]         # p
  # genotypes
  freqGenotypes <- array(data = 0,dim=c(nSnp,3))
  GenoFreq <- function(x,y){sum(x==y)/length(x)}
  counter=0
  for (iAllele in 1:2) {
    for (jAllele in iAllele:2) {
      counter=counter+1
      freqGenotypes[,counter] <- apply(genotemp, 2, GenoFreq,y=(counter-1))
    }
  }
  genotemp <- NULL
  
  ##### Functions and arrays needed to calculate stastical effects can be found in S1CalculateQtlEffects: 
  bkl <- array(data = NA,dim=c(6,nEpisPairLoci,nTrait))
  Wfkl <- setupWfklmat1()  # put here so code more efficient
  
  setupbkl <- function(traitno,pairno,pk,pl) {
    #traitno=1;pairno=1;pk=freqGenotypes[episPairLoci[pairno,1],];pl=freqGenotypes[episPairLoci[pairno,2],]
    
    qtlno_k=sum(SNPinfo[1:episPairLoci[pairno,1],4]=='T')
    qtlno_l=sum(SNPinfo[1:episPairLoci[pairno,2],4]=='T')
    
    # setup Ef manually like this or do it similarly as in setup Wsklmat or Wfkl for equivalent positions
    Ef=matrix(0,nrow = 6,ncol=1)
    Ef[,1]=c(1,aQtleffects[qtlno_k,traitno],aQtleffects[qtlno_l,traitno],domValues[qtlno_k,traitno],domValues[qtlno_l,traitno],episPairLociEffect[pairno,traitno])
    # Warning note: need to check if mean affect Es
    
    Wskl <- setupWsklmat(pk,pl)
    #Wfkl <- setupWfklmat1()
    
    temp1 <- t(Wskl)%*%Wskl
    #if (is.na(det(temp1))) {
    #	print(paste('One of loci is not segrating for qtl pair n.',pairno,'of qtl',qtlno_k, 'or', qtlno_l))
    #	stop('NOIA conversion is not possible. Because (px[3]+px[1]-((px[1]-px[3])^2))=0')
    #}
    if (det(temp1)<0.0000000001) diag(temp1)=diag(temp1)+0.00000001   # this might be needed for non-invertable matrix
    Es=solve(temp1)%*%t(Wskl)%*%Wfkl%*%Ef
    bkl[,pairno,traitno] <<- Es
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
  
  # transfer effects to the right arrays
  aQtleffects_Stat <- array(data = 0,dim=c(nQtl,nTrait))
  domValues_Stat <- array(data = 0,dim=c(nQtl,nTrait))
  episPairLociEffect_Stat <- array(data = 0,dim=c(nEpisPairLoci,nTrait))
  
  for (iTrait in 1:(nTrait)) {
    temp1=rep(0,nSnp)
    temp1[episPairLoci[,1]]=temp1[episPairLoci[,1]]+bkl[2,,iTrait]
    temp1[episPairLoci[,2]]=temp1[episPairLoci[,2]]+bkl[3,,iTrait]
    aQtleffects_Stat[,iTrait]=temp1[SNPinfo[,4]=='T']
    
    temp1=rep(0,nSnp)
    temp1[episPairLoci[,1]]=temp1[episPairLoci[,1]]+bkl[4,,iTrait]
    temp1[episPairLoci[,2]]=temp1[episPairLoci[,2]]+bkl[5,,iTrait]
    domValues_Stat[,iTrait]=temp1[SNPinfo[,4]=='T']
    
    episPairLociEffect_Stat[,iTrait]=bkl[6,,iTrait]
  }	
  
  ########### Report statistical genetic effects ################
  # stastical variances by id
  nAnimal= length(idtobeused)
  idalpha <- iddom <- idepis <- array(data = NA,dim=c(nAnimal,nTrait))
  v_alpha <- v_statdom <- v_statepis <- data.frame()
  # idalpha & iddom
  counter=0
  for (idd in idtobeused){
    counter=counter+1
    genotypecode=geno[idd,SNPinfo[,4]=='T']
    temp1 <- freqGenotypes[SNPinfo[,4]=='T',]
    px <- lapply(seq_len(nrow(temp1)), function(i) temp1[i,])
    haiVect <- mapply(findhai,px=px,allelecount=genotypecode)
    hdiVect <- mapply(findhdi,px=px,allelecount=genotypecode)
    for (iTrait in 1:(nTrait)) {
      idalpha[counter,iTrait]=sum(haiVect*aQtleffects_Stat[,iTrait])
      iddom[counter,iTrait]=sum(hdiVect*domValues_Stat[,iTrait])
    }		
    
    #idepis
    xAk <- geno[idd,episPairLoci[,1]]
    temp1 <- freqGenotypes[episPairLoci[,1],]
    pxAk <- lapply(seq_len(nrow(temp1)), function(i) temp1[i,])
    
    xAl <- geno[idd,episPairLoci[,2]]
    temp1 <- freqGenotypes[episPairLoci[,2],]
    pxAl <- lapply(seq_len(nrow(temp1)), function(i) temp1[i,])
    
    haihaiVect <- mapply(findhaihai,pk=pxAk,ack=xAk,pl=pxAl,acl=xAl)		
    for (iTrait in 1:(nTrait)) {
      idepis[counter,iTrait]=sum(haihaiVect*episPairLociEffect_Stat[,iTrait])
    }
  }
  # calculate covariances for traits
  v_alpha <- rbind(v_alpha,as.vector(var(idalpha)))
  v_statdom <- rbind(v_statdom,as.vector(var(iddom)))
  v_statepis <- rbind(v_statepis,as.vector(var(idepis)))
  
  # stastical variances by locus
  # additive
  varbylocus <- array(data = 0,dim=c(nTrait,nTrait,nQtl))
  px <- freqGenotypes[SNPinfo[,4]=='T',]
  tempfunc <- function(itemp) {
    varbylocus[,,itemp] <<- calculateVarbyLocus(EffInput=aQtleffects_Stat[itemp,],px=px[itemp,],
                                                hx=c(findhai(px=px[itemp,],allelecount=0),findhai(px=px[itemp,],allelecount=1),findhai(px=px[itemp,],allelecount=2)))
  }
  invisible(mapply(tempfunc,itemp=1:nQtl))    # using apply
  v_temp=array(data = 0,dim=c(nTrait,nTrait))
  for (iTrait in 1:(nTrait)) {
    for (jTrait in 1:(nTrait)) {
      v_temp[iTrait,jTrait]=sum(varbylocus[iTrait,jTrait,])
    }
  } 
  v_alpha <- rbind(v_alpha,as.vector(v_temp))
  # dominance
  varbylocus <- array(data = 0,dim=c(nTrait,nTrait,nQtl))
  px <- freqGenotypes[SNPinfo[,4]=='T',]
  tempfunc <- function(itemp) {
    varbylocus[,,itemp] <<- calculateVarbyLocus(EffInput=domValues_Stat[itemp,],px=px[itemp,],
                                                hx=c(findhdi(px=px[itemp,],allelecount=0),findhdi(px=px[itemp,],allelecount=1),findhdi(px=px[itemp,],allelecount=2)))
  }
  invisible(mapply(tempfunc,itemp=1:nQtl))    # using apply
  v_temp=array(data = 0,dim=c(nTrait,nTrait))
  for (iTrait in 1:(nTrait)) {
    for (jTrait in 1:(nTrait)) {
      v_temp[iTrait,jTrait]=sum(varbylocus[iTrait,jTrait,])
    }
  } 
  v_statdom <- rbind(v_statdom,as.vector(v_temp))
  # epis
  varbylocus <- array(data = 0,dim=c(nTrait,nTrait,nEpisPairLoci))
  tempfunc <- function(itemp) {
    varbylocus[,,itemp] <<- calculateVarbyLocusEpis(EffInput=episPairLociEffect_Stat[itemp,],
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
  v_statepis <- rbind(v_statepis,as.vector(v_temp))
  
  # naming
  names(v_alpha) <- paste0('val',1:dim(v_alpha)[2])
  v_alpha = data.frame(effect='stat',covname='v_alpha',method=1:dim(v_alpha)[1],v_alpha)
  names(v_statdom) <- paste0('val',1:dim(v_statdom)[2])
  v_statdom = data.frame(effect='stat',covname='v_statdom',method=1:dim(v_statdom)[1],v_statdom)
  names(v_statepis) <- paste0('val',1:dim(v_statepis)[2])
  v_statepis = data.frame(effect='stat',covname='v_statepis',method=1:dim(v_statepis)[1],v_statepis)
  
  # return
  if (returnvalues=='variances') 	{   return(do.call("rbind", list(v_alpha,v_statdom,v_statepis))) 	}  
  if (returnvalues=='effects') 	    {   return(do.call("cbind", list(as.numeric(idtobeused),idalpha,iddom,idepis)))	}
  if (returnvalues=='both')         {   return(list(do.call("rbind", list(v_alpha,v_statdom,v_statepis)),
                                                    do.call("cbind", list(as.numeric(idtobeused),idalpha,iddom,idepis))))  }
  
}

######################################################################################################
####### Functions to calculate biological effects of additive genetics, dominance & epistasis ########
### given we have dataset haplo_mat, haplo_pat and arrays of biological effects (a,d,epis)############
################ domValues; aQtleffects; episPairLoci; episPairLociEffect ############################
calculateBiolEffectOne <- function(genodf=NULL,idtobeused=NULL,returnvalues='variances') {   
  # genodf=NULL;idtobeused=pop$id[pop$currentpop==1];rounddigits=3; returnvalues='variances'
  # check inputs
  if (!returnvalues %in% c('variances','effects','both')) returnvalues='variances'
  
  if (!is.null(genodf)) {  # if genotype file given, check correct formats
    genodf <- as.matrix(genodf)
  } else {
    genodf = haplo_mat + haplo_pat - 2
    genodf = cbind(matrix(c(1:length(genodf[,1])),nrow=length(genodf[,1])),genodf)  # add column ID
  }
  if (is.null(idtobeused))             {  idtobeused = genodf[,1]           }
  if ((length(genodf[1,])-1)!=nSnp) stop(paste('Dimension of genodf is not right. Ncolumn must be:',nSnp+1))
  if (any(!genodf[,2:length(genodf[1,])] %in% c(0,1,2))) stop('Genotype codes are not in correct format. It must be 0,1,2.')
  if (length(idtobeused)!=length(unique(idtobeused))) {
    idtobeused=unique(idtobeused); 	  print('Repeated ID are removed.')
  }   
  # keep needed genotypes only
  allID <- unique(idtobeused)
  allID <- allID[order(allID)]   
  genodf <- genodf[genodf[,1] %in% allID,]
  
  # set rownames of geno as id 
  geno <- genodf[,-1]
  rownames(geno) <- genodf[,1] ; genodf <- NULL
  idtobeused=as.character(idtobeused)
  
  ### calculate frequency of alleles and genotypes
  genotemp <- geno
  # alleles
  freqAlleles <- array(data = 0,dim=c(nSnp,2))
  MAF <- function(x){sum(x,na.rm=T)/(sum(!is.na(x))*2)}
  freqAlleles[,2]=apply(genotemp, 2, MAF)   # q
  freqAlleles[,1]=1-freqAlleles[,2]         # p
  # genotypes
  freqGenotypes <- array(data = 0,dim=c(nSnp,3))
  GenoFreq <- function(x,y){sum(x==y)/length(x)}
  counter=0
  for (iAllele in 1:2) {
    for (jAllele in iAllele:2) {
      counter=counter+1
      freqGenotypes[,counter] <- apply(genotemp, 2, GenoFreq,y=(counter-1))
    }
  }
  genotemp <- NULL
  
  ########### Report biological genetic effects ################
  nAnimal= length(idtobeused)
  ida <- iddom <- idepis <- idtgv <- array(data = NA,dim=c(nAnimal,nTrait))
  v_a <- v_dom <- v_epis <- data.frame()
  
  counter=0
  for (idd in idtobeused){
    counter=counter+1
    genotypecode=geno[idd,SNPinfo[,4]=='T']
    # a
    for (iTrait in 1:(nTrait)) {
      ida[counter,iTrait]=sum(aQtleffects[,iTrait] * (genotypecode-1))
    }
    # dom
    for (iTrait in 1:(nTrait)) {
      iddom[counter,iTrait]=sum(domValues[genotypecode==1,iTrait])
    }
    
    #idepis
    xAk <- geno[idd,episPairLoci[,1]]-1
    xAl <- geno[idd,episPairLoci[,2]]-1
    for (iTrait in 1:(nTrait)) {
      idepis[counter,iTrait]=sum(episPairLociEffect[1:nEpisPairLoci,iTrait]*xAk*xAl)
    }
    
    #idtgv
    #for (iTrait in 1:(nTrait)) {
    #  idtgv[counter,iTrait]=ida[counter,iTrait]+iddom[counter,iTrait]+idepis[counter,iTrait]
    #}	  
  }
  # calculate covariances for traits
  v_a <- rbind(v_a,as.vector(var(ida)))
  v_dom <- rbind(v_dom,as.vector(var(iddom)))
  v_epis <- rbind(v_epis,as.vector(var(idepis)))
  
  # stastical variances by locus
  # additive
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
  v_a <- rbind(v_a,as.vector(v_temp))
  # dominance
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
  v_dom <- rbind(v_dom,as.vector(v_temp))
  # epis
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
  v_epis <- rbind(v_epis,as.vector(v_temp))
  
  # naming
  names(v_a) <- paste0('val',1:dim(v_a)[2])
  v_a = data.frame(effect='biol',covname='v_a',method=1:dim(v_a)[1],v_a)
  names(v_dom) <- paste0('val',1:dim(v_dom)[2])
  v_dom = data.frame(effect='biol',covname='v_dom',method=1:dim(v_dom)[1],v_dom)
  names(v_epis) <- paste0('val',1:dim(v_epis)[2])
  v_epis = data.frame(effect='biol',covname='v_epis',method=1:dim(v_epis)[1],v_epis)
  
  # return
  if (returnvalues=='variances') 	{   return(do.call("rbind", list(v_a,v_dom,v_epis))) 	}  
  if (returnvalues=='effects') 	    {   return(do.call("cbind", list(as.numeric(idtobeused),ida,iddom,idepis)))	}
  if (returnvalues=='both')         {   return(list(do.call("rbind", list(v_a,v_dom,v_epis)),
                                                    do.call("cbind", list(as.numeric(idtobeused),ida,iddom,idepis))))  }
}

######################################################################################################
## Plot obs: Functions to calculate statistical effects of additive genetics, dominance & epistasis ##
### given we have dataset haplo_mat, haplo_pat and arrays of biological effects (a,d,epis)############
################ domValues; aQtleffects; episPairLoci; episPairLociEffect ############################
calculateStatEffectPlotOne <- function(genoPar1=NULL,idPar1=NULL,genoPar2=NULL,idPar2=NULL,returnvalues='variances',
                                       genoPar1Freq=NULL,idPar1Freq=NULL,genoPar2Freq=NULL,idPar2Freq=NULL) {
  # genoPar1=NULL;genoPar2=NULL;returnvalues='variances';genoPar1Freq=NULL;genoPar2Freq=NULL;idPar1=NULL;idPar2=NULL;idPar1Freq=NULL;idPar2Freq=NULL
  # check inputs
  if (!returnvalues %in% c('variances','effects','both')) returnvalues='variances'
  
  # genotype of parents
  if ((is.null(genoPar1) & is.null(idPar1))| (is.null(genoPar2) & is.null(idPar2))) {
    stop(paste('Genotype files or id of parents must be provided.'))   }  
  if (is.null(genoPar1)) {		genoPar1 <- haplo_mat[idPar1,] + haplo_pat[idPar1,] - 2  }
  if (is.null(genoPar2)) {		genoPar2 <- haplo_mat[idPar2,] + haplo_pat[idPar2,] - 2  }
  genoPar1 <- as.matrix(genoPar1); genoPar2 <- as.matrix(genoPar2)
  
  # genotype for freq calculation
  if (is.null(genoPar1Freq) & is.null(idPar1Freq)) genoPar1Freq = genoPar1
  if (is.null(genoPar2Freq) & is.null(idPar2Freq)) genoPar2Freq = genoPar2
  if (is.null(genoPar1Freq)) genoPar1Freq = haplo_mat[idPar1Freq,] + haplo_pat[idPar1Freq,] - 2
  if (is.null(genoPar2Freq)) genoPar2Freq = haplo_mat[idPar2Freq,] + haplo_pat[idPar2Freq,] - 2
  genoPar1Freq <- as.matrix(genoPar1Freq); genoPar2Freq <- as.matrix(genoPar2Freq)
  
  if (any(!unique(genoPar1Freq) %in% c(0,1,2))) stop(paste('Genotype code of genoPar1Freq can be 0,1,2 only..'))
  if (any(!unique(genoPar2Freq) %in% c(0,1,2))) stop(paste('Genotype code of genoPar2Freq can be 0,1,2 only..'))
  
  if ((length(genoPar1[1,]))!=nSnp) stop(paste('Dimension of genoPar1 is not right. Ncolumn must be:',nSnp))
  if ((length(genoPar2[1,]))!=nSnp) stop(paste('Dimension of genoPar2 is not right. Ncolumn must be:',nSnp))
  
  if (dim(genoPar1)[1]!=dim(genoPar2)[1]) stop(paste('Dimension of genoPar1 must be the same as genoPar2.'))
  nAnimal=dim(genoPar1)[1]
  if (nAnimal<2) stop(paste('Genotype file is empty.'))
  
  Va_Par1 <- Va_Par2 <- Vd_plot <- Vaa_plot <- data.frame()
  
  ###### calculate allele frequency
  MAF <- function(x){sum(x,na.rm=T)/(sum(!is.na(x))*2)}
  GenoFreq <- function(x,y){sum(x==y)/length(x)}
  # Frequency based on genoPar1Freq
  genotemp = genoPar1Freq   # frequencyBasedPop
  # alleles
  freqAlleles <- array(data = 0,dim=c(nSnp,2))
  #MAF <- function(x){sum(x,na.rm=T)/(sum(!is.na(x))*2)}
  freqAlleles[,2]=apply(genotemp, 2, MAF)   # q
  freqAlleles[,1]=1-freqAlleles[,2]         # p
  # genotypes
  freqGenotypes <- array(data = 0,dim=c(nSnp,3))
  #GenoFreq <- function(x,y){sum(x==y)/length(x)}
  counter=0
  for (iAllele in 1:2) {
    for (jAllele in iAllele:2) {
      counter=counter+1
      freqGenotypes[,counter] <- apply(genotemp, 2, GenoFreq,y=(counter-1))
    }
  }
  freqAlleles_Par1 <- freqAlleles
  freqGenotypes_Par1 <- freqGenotypes
  genotemp <- freqAlleles <- freqGenotypes <- NULL
  
  # Frequency based on genoPar2Freq
  genotemp = genoPar2Freq   # frequencyBasedPop
  # alleles
  freqAlleles <- array(data = 0,dim=c(nSnp,2))
  #MAF <- function(x){sum(x,na.rm=T)/(sum(!is.na(x))*2)}
  freqAlleles[,2]=apply(genotemp, 2, MAF)   # q
  freqAlleles[,1]=1-freqAlleles[,2]         # p
  # genotypes
  freqGenotypes <- array(data = 0,dim=c(nSnp,3))
  #GenoFreq <- function(x,y){sum(x==y)/length(x)}
  counter=0
  for (iAllele in 1:2) {
    for (jAllele in iAllele:2) {
      counter=counter+1
      freqGenotypes[,counter] <- apply(genotemp, 2, GenoFreq,y=(counter-1))
    }
  }
  freqAlleles_Par2 <- freqAlleles
  freqGenotypes_Par2 <- freqGenotypes
  genotemp <- freqAlleles <- freqGenotypes <- NULL
  
  ########## 
  
  # setup bklplot1 using ivares-Castro and Calborg2007 method, to calculate statistical genetic effects of alpha from functional effects of a, d and e
  # only work if (nQtlEachQtlInteractwith==1)
  if (nQtlEachQtlInteractwith!=1) stop ('this method works with nQtlEachQtlInteractwith=1 only.') 
  Wfkl=setupWfklmat1(tx1=c(0,1,2),tx2=c(0,0.5,1,1.5,2))
  bklplot1 <- array(data = NA,dim=c(8,nEpisPairLoci,nTrait))
  setupbklmatPlot <- function(traitno,pairno,ap1_k, ap2_k,ap1_l, ap2_l) {
    #traitno=1;pairno=1;ap1_k=freqAlleles_Par1[episPairLoci[pairno,1],1]; ap2_k=freqAlleles_Par2[episPairLoci[pairno,1],1]
    #ap1_l=freqAlleles_Par1[episPairLoci[pairno,2],1]; ap2_l=freqAlleles_Par2[episPairLoci[pairno,2],1]
    
    qtlno_k=sum(SNPinfo[1:episPairLoci[pairno,1],4]=='T')
    qtlno_l=sum(SNPinfo[1:episPairLoci[pairno,2],4]=='T')
    
    Ef=matrix(0,nrow = 8,ncol=1) # Try ref point of 1. But ref point did not change anything, so leave 0 here.
    Wx_k=c(1,aQtleffects[qtlno_k,traitno],aQtleffects[qtlno_k,traitno],domValues[qtlno_k,traitno])
    Wx_l=c(1,aQtleffects[qtlno_l,traitno],aQtleffects[qtlno_l,traitno],domValues[qtlno_l,traitno])
    # a1
    temp1=kronecker(Wx_l[c(1,2)],Wx_k[c(1,2)]) 
    Ef[2:3,]=temp1[2:3]
    # a2
    temp1=kronecker(Wx_l[c(1,3)],Wx_k[c(1,3)]) 
    Ef[4:5,]=temp1[2:3]
    # d3way
    temp1=kronecker(Wx_l[c(1,4)],Wx_k[c(1,4)]) 
    Ef[6:7,]=temp1[2:3]
    # aa3way
    Ef[8,]= episPairLociEffect[pairno,traitno]
    
    #Wfkl=setupWfklmat1()
    Wskl=setupWsklPlotmat(tx1=c(0,1,2),ap1_k,ap2_k,tx2=c(0,0.5,1,1.5,2),ap1_l,ap2_l)
    temp1 <- t(Wskl)%*%Wskl
    if (det(temp1)<0.0000000001) diag(temp1)=diag(temp1)+0.00000001   # this might be needed for non-invertable matrix
    bklmet=solve(temp1)%*%t(Wskl)%*%Wfkl%*%Ef
    bklplot1[,pairno,traitno] <<- bklmet
    #temp2=Wfkl%*%Ef
    #temp3=Wskl%*%bklmet
    #temp4=data.frame(temp2,temp3)
  }
  # preparation for apply functions
  ap1_k <- freqAlleles_Par1[episPairLoci[,1],1]    # not sure 1 or 2, try 1 here
  ap2_k <- freqAlleles_Par2[episPairLoci[,1],1]    # not sure 1 or 2, try 1 here
  
  ap1_l <- freqAlleles_Par1[episPairLoci[,2],1]    # not sure 1 or 2, try 1 here
  ap2_l <- freqAlleles_Par2[episPairLoci[,2],1]    # not sure 1 or 2, try 1 here
  
  # apply funtions for caculate bklplot1
  for (iTrait in 1:(nTrait)) {
    invisible(mapply(setupbklmatPlot,traitno=rep(iTrait,nEpisPairLoci),pairno=1:nEpisPairLoci,
                     ap1_k=ap1_k, ap2_k=ap2_k, ap1_l=ap1_l, ap2_l=ap2_l))
  }
  
  # Transfer bklplot1 to easier handle vector of stastitical effects 
  statistical_a1 = statistical_a2 = statistical_d12 = array(data = 0,dim=c(nQtl,nTrait))
  statistical_aa12 = array(data = 0,dim=c(nEpisPairLoci,nTrait))  
  for (iTrait in 1:(nTrait)) {
    temp1=rep(0,nSnp)
    temp1[episPairLoci[,1]]=temp1[episPairLoci[,1]]+bklplot1[2,,iTrait]
    temp1[episPairLoci[,2]]=temp1[episPairLoci[,2]]+bklplot1[3,,iTrait]
    statistical_a1[,iTrait]=temp1[SNPinfo[,4]=='T']
    
    temp1=rep(0,nSnp)
    temp1[episPairLoci[,1]]=temp1[episPairLoci[,1]]+bklplot1[4,,iTrait]
    temp1[episPairLoci[,2]]=temp1[episPairLoci[,2]]+bklplot1[5,,iTrait]
    statistical_a2[,iTrait]=temp1[SNPinfo[,4]=='T']
    
    temp1=rep(0,nSnp)
    temp1[episPairLoci[,1]]=temp1[episPairLoci[,1]]+bklplot1[6,,iTrait]
    temp1[episPairLoci[,2]]=temp1[episPairLoci[,2]]+bklplot1[7,,iTrait]
    statistical_d12[,iTrait]=temp1[SNPinfo[,4]=='T']
    
    statistical_aa12[,iTrait]=bklplot1[8,,iTrait]
    
  }	
  # Note that Here, statistical_aa12= episPairLociEffect; and statistical_d12=domValues. But use statistical_aa12 & statistical_d12 for ease to understand
  
  ### calculate statistical variances of par1 pops contributed to plot 
  idgv <- array(data = NA,dim=c(nAnimal,nTrait))
  counter=0
  for (idd in 1:nAnimal) {
    counter=counter+1
    genotypecode=(genoPar1[idd,SNPinfo[,4]=='T'])
    temp1 <- freqGenotypes_Par1[SNPinfo[,4]=='T',]
    px <- lapply(seq_len(nrow(temp1)), function(i) temp1[i,])
    haiVect <- mapply(findhai,px=px,allelecount=genotypecode)/2       # devided by 2 because half of genes transfered. 
    for (iTrait in 1:(nTrait)) {
      idgv[counter,iTrait]=sum(haiVect*statistical_a1[,iTrait])
    }	  
  }
  Va_Par1 <- rbind(Va_Par1,as.vector(var(idgv)))
  
  ### calculate statistical variances of par2 pops contributed to plot 
  idgv <- array(data = NA,dim=c(nAnimal,nTrait))
  counter=0
  for (idd in 1:nAnimal) {
    counter=counter+1
    genotypecode=(genoPar2[idd,SNPinfo[,4]=='T'])
    temp1 <- freqGenotypes_Par2[SNPinfo[,4]=='T',]
    px <- lapply(seq_len(nrow(temp1)), function(i) temp1[i,])
    haiVect <- mapply(findhai,px=px,allelecount=genotypecode)/2       # devided by 2 because half of genes transfered. 
    for (iTrait in 1:(nTrait)) {
      idgv[counter,iTrait]=sum(haiVect*statistical_a2[,iTrait])
    }	  
  }
  Va_Par2 <- rbind(Va_Par2,as.vector(var(idgv)))
  
  #### Dominance
  idgv <- array(data = NA,dim=c(nAnimal,nTrait))
  counter=0
  for (idd in 1:nAnimal) {
    counter=counter+1
    allelecount1=genoPar1[idd,SNPinfo[,4]=='T']
    allelecount2=genoPar2[idd,SNPinfo[,4]=='T']
    p1=freqAlleles_Par1[SNPinfo[,4]=='T',1]    # p or q, not sure.
    p2=freqAlleles_Par2[SNPinfo[,4]=='T',1]    # p or q, not sure.
    hdiVect <- mapply(findhdiPlot,p1=p1,p2=p2,
                      allelecount1=allelecount1,allelecount2=allelecount2)
    for (iTrait in 1:(nTrait)) {
      idgv[counter,iTrait]=sum(hdiVect*statistical_d12[,iTrait])  # domValues=statistical_d12, so dont need to be complicated
    }	  
  }
  Vd_plot <- rbind(Vd_plot,as.vector(var(idgv)))
  
  #### epistasis (aa)
  idgv <- array(data = NA,dim=c(nAnimal,nTrait))
  # preparation for apply functions
  ap1_k <- freqAlleles_Par1[episPairLoci[,1],1]    # not sure 1 or 2, try 1 here
  ap2_k <- freqAlleles_Par2[episPairLoci[,1],1]    # not sure 1 or 2, try 1 here
  ap1_l <- freqAlleles_Par1[episPairLoci[,2],1]    # not sure 1 or 2, try 1 here
  ap2_l <- freqAlleles_Par2[episPairLoci[,2],1]    # not sure 1 or 2, try 1 here
  
  counter=0
  for (idd in 1:nAnimal) {
    counter=counter+1
    ac1_k=genoPar1[idd,episPairLoci[,1]]
    ac2_k=genoPar2[idd,episPairLoci[,1]]
    
    ac1_l=genoPar1[idd,episPairLoci[,2]]
    ac2_l=genoPar2[idd,episPairLoci[,2]]
    
    haihaiPlotVect <- mapply(findhaihaiPlot,ap1_k=ap1_k,ap2_k=ap2_k,ac1_k=ac1_k,ac2_k=ac2_k,
                             ap1_l=ap1_l,ap2_l=ap2_l,ac1_l=ac1_l,ac2_l=ac2_l)
    
    for (iTrait in 1:(nTrait)) {
      idgv[counter,iTrait]=sum(haihaiPlotVect*statistical_aa12[,iTrait])  # Note that Here, statistical_aa12= episPairLociEffect
    }	  
  }
  Vaa_plot <- rbind(Vaa_plot,as.vector(var(idgv)))
  
  # stastical variances by locus
  # additive par1
  varbylocus <- array(data = 0,dim=c(nTrait,nTrait,nQtl))
  px <- freqGenotypes_Par1[SNPinfo[,4]=='T',]
  tempfunc <- function(itemp) {
    varbylocus[,,itemp] <<- calculateVarbyLocus(EffInput=statistical_a1[itemp,],px=px[itemp,],
                                                hx=c(findhai(px=px[itemp,],allelecount=0),findhai(px=px[itemp,],allelecount=1),findhai(px=px[itemp,],allelecount=2))/2)
  }
  invisible(mapply(tempfunc,itemp=1:nQtl))    # using apply
  v_temp=array(data = 0,dim=c(nTrait,nTrait))
  for (iTrait in 1:(nTrait)) {
    for (jTrait in 1:(nTrait)) {
      v_temp[iTrait,jTrait]=sum(varbylocus[iTrait,jTrait,])
    }
  } 
  Va_Par1 <- rbind(Va_Par1,as.vector(v_temp))
  
  # additive par2
  varbylocus <- array(data = 0,dim=c(nTrait,nTrait,nQtl))
  px <- freqGenotypes_Par2[SNPinfo[,4]=='T',]
  tempfunc <- function(itemp) {
    varbylocus[,,itemp] <<- calculateVarbyLocus(EffInput=statistical_a1[itemp,],px=px[itemp,],
                                                hx=c(findhai(px=px[itemp,],allelecount=0),findhai(px=px[itemp,],allelecount=1),findhai(px=px[itemp,],allelecount=2))/2)
  }
  invisible(mapply(tempfunc,itemp=1:nQtl))    # using apply
  v_temp=array(data = 0,dim=c(nTrait,nTrait))
  for (iTrait in 1:(nTrait)) {
    for (jTrait in 1:(nTrait)) {
      v_temp[iTrait,jTrait]=sum(varbylocus[iTrait,jTrait,])
    }
  } 
  Va_Par2 <- rbind(Va_Par2,as.vector(v_temp))
  
  #### Dominance
  varbylocus <- array(data = 0,dim=c(nTrait,nTrait,nQtl))
  p1 <- freqGenotypes_Par1[SNPinfo[,4]=='T',]
  p2 <- freqGenotypes_Par2[SNPinfo[,4]=='T',]
  ap1 <- freqAlleles_Par1[SNPinfo[,4]=='T',1]
  ap2 <- freqAlleles_Par2[SNPinfo[,4]=='T',1]
  tempfunc <- function(itemp) {
    varbylocus[,,itemp] <<- calculateStatVarDPlotbyLocus(EffInput=statistical_d12[itemp,],ap1=ap1[itemp],ap2=ap2[itemp],
                                                         p1=p1[itemp,],p2=p2[itemp,])
  }
  invisible(mapply(tempfunc,itemp=1:nQtl))    # using apply
  v_temp=array(data = 0,dim=c(nTrait,nTrait))
  for (iTrait in 1:(nTrait)) {
    for (jTrait in 1:(nTrait)) {
      v_temp[iTrait,jTrait]=sum(varbylocus[iTrait,jTrait,])
    }
  } 
  Vd_plot <- rbind(Vd_plot,as.vector(v_temp))
  
  #### epistasis (aa)
  varbylocus <- array(data = 0,dim=c(nTrait,nTrait,nEpisPairLoci))
  tempfunc <- function(itemp) {
    varbylocus[,,itemp] <<- calculateVarPlotbyLocusEpis(EffInput=statistical_aa12[itemp,],
                                                        apk1=freqAlleles_Par1[episPairLoci[itemp,1],1],apk2=freqAlleles_Par2[episPairLoci[itemp,1],1],
                                                        pk1=freqGenotypes_Par1[episPairLoci[itemp,1],],pk2=freqGenotypes_Par2[episPairLoci[itemp,1],],
                                                        apl1=freqAlleles_Par1[episPairLoci[itemp,2],1],apl2=freqAlleles_Par2[episPairLoci[itemp,2],1],
                                                        pl1=freqGenotypes_Par1[episPairLoci[itemp,2],],pl2=freqGenotypes_Par2[episPairLoci[itemp,2],])
  }
  invisible(mapply(tempfunc,itemp=1:nEpisPairLoci))    # using apply
  v_temp=array(data = 0,dim=c(nTrait,nTrait))
  for (iTrait in 1:(nTrait)) {
    for (jTrait in 1:(nTrait)) {
      v_temp[iTrait,jTrait]=sum(varbylocus[iTrait,jTrait,])
    }
  } 
  Vaa_plot <- rbind(Vaa_plot,as.vector(v_temp))
  
  #
  names(Va_Par1) <- paste0('val',1:dim(Va_Par1)[2])
  Va_Par1 = data.frame(effect='stat',covname='Va_Par1',method=1:dim(Va_Par1)[1],Va_Par1)
  names(Va_Par2) <- paste0('val',1:dim(Va_Par2)[2])
  Va_Par2 = data.frame(effect='stat',covname='Va_Par2',method=1:dim(Va_Par2)[1],Va_Par2)
  names(Vd_plot) <- paste0('val',1:dim(Vd_plot)[2])
  Vd_plot = data.frame(effect='stat',covname='Vd_plot',method=1:dim(Vd_plot)[1],Vd_plot)
  names(Vaa_plot) <- paste0('val',1:dim(Vaa_plot)[2])
  Vaa_plot = data.frame(effect='stat',covname='Vaa_plot',method=1:dim(Vaa_plot)[1],Vaa_plot)
  
  # return
  if (returnvalues=='variances') 	{   return(do.call("rbind", list(Va_Par1,Va_Par2,Vd_plot,Vaa_plot))) 	}  
}
