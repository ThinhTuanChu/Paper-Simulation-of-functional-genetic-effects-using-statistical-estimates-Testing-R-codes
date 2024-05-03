######### Generate DMU inputs for estimates ########
##### make phenotype files for testing with DMU after optimization done and new aQtleffects, domValues & episPairLociEffect generated
# update phenotype with new aQtleffects, domValues & episPairLociEffect generated

if (!exists('iRepScheme')) iRepScheme=1

# Use plot observation rather than individual obs. Plot observation is the average from parents, not its own phenotype.
if (!is.null(idplotobs)) {
  ntemp=length(idplotobs)
  idplotobs <- idplotobs[order(idplotobs)]
  idplotobs <- idplotobs[idplotobs %in% c((nFounderAnimal*nFounderPop+1):maxid)]
  idplotobs <- idplotobs[order(idplotobs)]
  if (ntemp != length(idplotobs)) {       print(paste('Number of ids that plot observation could not be calculated:',ntemp-length(idplotobs)))   }    
  if (length(idplotobs)>0) {
    for (idd in idplotobs) {
      # additive
      temp1=((haplo_pat[pop$par1[idd],SNPinfo[,4]=='T']+haplo_mat[pop$par1[idd],SNPinfo[,4]=='T']-2 +
                haplo_pat[pop$par2[idd],SNPinfo[,4]=='T']+haplo_mat[pop$par2[idd],SNPinfo[,4]=='T']-2)/2) - 1
      for (iTrait in 1:(nTrait)) { #iTrait=1
        pop[idd,paste0('trait_av',iTrait)]=sum(aQtleffects[,iTrait] * temp1)             
      }
      # dominance
      temp1=haplo_pat[pop$par1[idd],SNPinfo[,4]=='T']+haplo_mat[pop$par1[idd],SNPinfo[,4]=='T']-2
      temp2=haplo_pat[pop$par2[idd],SNPinfo[,4]=='T']+haplo_mat[pop$par2[idd],SNPinfo[,4]=='T']-2
      temp0=rep(NA,nQtl)
      temp0[temp1== 0 & temp2== 0] <- 0.0
      temp0[temp1== 0 & temp2== 1] <- 0.5 
      temp0[temp1== 0 & temp2== 2] <- 1.0 
      temp0[temp1== 1 & temp2== 0] <- 0.5 
      temp0[temp1== 1 & temp2== 1] <- 0.5 
      temp0[temp1== 1 & temp2== 2] <- 0.5 
      temp0[temp1== 2 & temp2== 0] <- 1.0
      temp0[temp1== 2 & temp2== 1] <- 0.5 
      temp0[temp1== 2 & temp2== 2] <- 0.0 		 
      pop[idd,paste0('trait_dv',iTrait)]=sum(domValues[,iTrait] * temp0) 
      
      # epistasis: a x a
      if (!is.null(episPairLociEffect)) {
        
        xAk <- (((haplo_pat[pop$par1[idd],episPairLoci[,1]] + haplo_mat[pop$par1[idd],episPairLoci[,1]] - 2)+
                   (haplo_pat[pop$par2[idd],episPairLoci[,1]] + haplo_mat[pop$par2[idd],episPairLoci[,1]] - 2))/2)-1
        xAl <- (((haplo_pat[pop$par1[idd],episPairLoci[,2]] + haplo_mat[pop$par1[idd],episPairLoci[,2]] - 2)+
                   (haplo_pat[pop$par2[idd],episPairLoci[,2]] + haplo_mat[pop$par2[idd],episPairLoci[,2]] - 2))/2)-1
        for (iTrait in 1:(nTrait)) {
          pop[idd,paste0('trait_episgv',iTrait)]=sum(episPairLociEffect[1:nEpisPairLoci,iTrait]*xAk*xAl)
        }	  
      }
      
      # recalculate phenotype
      for (iTrait in 1:(nTrait)) {
        pop[idd,paste0('trait_tgv',iTrait)]=pop[idd,paste0('trait_av',iTrait)]+pop[idd,paste0('trait_dv',iTrait)]+
          pop[idd,paste0('trait_episgv',iTrait)]
        pop[idd,paste0('trait_obs',iTrait)]=pop[idd,paste0('trait_tgv',iTrait)]+pop[idd,paste0('trait_residual',iTrait)]
      }
    }
    print('Plot observations realized.')
  }
}

reportVarScheme=FALSE
if (reportVarScheme) {
  # check true variances by pop
  filename='varSchemTrue.res'
  if (!file.exists(filename)) write(paste('iRep','pop','freqPop','effect','covname','method',paste0(paste0('val',c(1:(nTrait*nTrait))),collapse ='\t')),file=filename,append=F)
  for (ipop in unique(pop$currentpop)) {
    idLst <- pop[pop$fixedFactor1>0 & pop$currentpop %in% c(ipop),'id']
	if (length(idLst)<10) next
    covname <- c('trait_av','trait_dv','trait_episgv','trait_residual')
    for (icovname in covname) {
      covnameTemp <- c()
      for (i in 1:nTrait) {
        covnameTemp <- c(covnameTemp,paste0(icovname,i))
      }
      temp1 <- var(pop[idLst,covnameTemp])
      #cat(paste("\n All generation: Pop",ipop,"for biological",icovname,'\n'))
      #print(temp1)#;print(round(cov2cor(temp1),digits=3))
      write(paste(iRepScheme,ipop,ipop,'biol',icovname,1,paste0(temp1,collapse ='\t')),file=filename,append=file.exists(filename))
    }
  }
  
  # True Statistical Variances
  for (ipop in unique(pop$currentpop)) {
    #print(paste('True statistical variances of pop:',ipop))
    idLst <- pop[pop$fixedFactor1>0 & pop$currentpop %in% c(ipop),'id']
	if (length(idLst)<10) next
    # True Statistical Variances
    temp1=calculateStatEffectOne(idtobeused=idLst,returnvalues='variances') #; print(temp1)
    temp1 <- data.frame(iRepScheme,pop=ipop,freqPop=ipop,temp1)
	write.table(temp1,file=filename,sep = "\t",row.names=F,col.names=!file.exists(filename),quote = FALSE,append=file.exists(filename))
  }
  
  # True Statistical Variances based on pop4
  idLst2 <- pop[pop$fixedFactor1>0 & pop$currentpop %in% c(4),'id']
  for (ipop in 1:3) {
    #print(paste('True statistical variances of pop:',ipop))
    idLst1 <- pop[pop$fixedFactor1>0 & pop$currentpop %in% c(ipop),'id']
	if (length(idLst)<10 | length(idLst2)<10) next
    # True Statistical Variances
    temp1=calculateStatEffectOne(idtobeused=idLst1,idfrequencyCalculation=idLst2,returnvalues='variances')
    temp1 <- data.frame(iRepScheme,pop=ipop,freqPop=4,temp1)
	write.table(temp1,file=filename,sep = "\t",row.names=F,col.names=!file.exists(filename),quote = FALSE,append=file.exists(filename))
  }
}

########## make phenotype files, different Gmatrix for pops ##############
## load Gmat calculation
filename <- paste0('C:/ThinhInDenmark/Rpreloadfuncs/',"BlupPreloadfunctions2.R")
if (!file.exists(filename)) filename <- paste0('/usr/home/qgg/chuthinh/adam/',"BlupPreloadfunctions2.R")
if (!file.exists(filename)) {  write(paste('Error: The file',filename,'does not exists. This module is needed for',filename),file=outlog,append=file.exists(outlog))}
source(filename)

#pheno
idLst <- pop[pop$fixedFactor1>0 & pop$currentpop %in% c(4),'id']
pheno <- generatePhenoData(idtobeused=idLst)
temp1 <- unique(pheno$par2)
temp1 <- pop[pop$id %in% temp1,c('id','par1','par2')]
names(temp1) <- c('par2','par21','par22')
temp2 <- merge(pheno,temp1,by='par2',all.x = T)
temp1 <- c('id','par1','par2','par21','par22')        # id=3-way hybrid, par1=Res, par2=msNR, par21=NR, par22=ms
temp11 <- names(pheno)[!names(pheno) %in% temp1]
temp1 <- c(temp1,temp11)
pheno <- temp2[,temp1]
write.table(pheno,file =paste0('dmudatInbredPop',4,'rep',iRepScheme,'.obs'),col.names = F,row.names = F,quote = FALSE,na="-9999",sep =" ")

# geno
idLst <- pop[pop$fixedFactor1>0 & pop$currentpop %in% c(1:4),'id']
geno <- generateGenoData(idtobeused=idLst,genoDataFormat=3)

################################### Different GRM ###################################
## gmatrices in statistical models
# ms NR additive
idLst <- pop[pop$fixedFactor1>0 & pop$currentpop %in% c(2),'id']
#idFreq <- pop[pop$fixedFactor1==1 & pop$currentpop %in% c(2),'id']
gmatTemp <- calculateGmatbyRmths1(genodf=geno,idtobeused=idLst,igmatrixfile=paste0('gmatinvInbredPop',2,'.add'),
                                  methodgmat='VanRaden1.add',scaleCov=1,returngmat='none',scaleDiag=TRUE)
# Restorer
idLst <- pop[pop$fixedFactor1>0 & pop$currentpop %in% c(1),'id']
#idFreq <- pop[pop$fixedFactor1==1 & pop$currentpop %in% c(1),'id']
gmatTemp <- calculateGmatbyRmths1(genodf=geno,idtobeused=idLst,igmatrixfile=paste0('gmatinvInbredPop',1,'.add'),
                                  methodgmat='VanRaden1.add',scaleCov=1,returngmat='none',scaleDiag=TRUE)

# 3way hybrid ms NR Res: Epis
pedTemp <- pheno[,c('id','par1','par2','par21','par22')]
genoMean3way=(haplo_pat[pedTemp$par1,]+haplo_mat[pedTemp$par1,]-2)/2 +
  ((haplo_pat[pedTemp$par21,]+haplo_mat[pedTemp$par21,]-2)/2 + 
     (haplo_pat[pedTemp$par22,]+haplo_mat[pedTemp$par22,]-2)/2)/2
genoMean3way=cbind(pedTemp[,'id'],genoMean3way)
gmatTemp <- calculateGmatbyRmths1(genodf=genoMean3way,methodgmat='VanRaden1.add',scaleCov=1,returngmat='gmat',scaleDiag=TRUE)
calculateEpisGmatbyRmths1(gmat1=gmatTemp,igmatrixfile=paste0('gmatinvInbredPop4threewaymean.aa'),scaleDiag=TRUE)
gmatTemp <- NULL

# 3way hybrid ms NR Res: dominance by PSK
calculatedominancegmat3wayPSK(genodf=geno,pedigree=pedTemp,igmatrixfile=paste0('gmatinvInbredPop',4,'threeway.PSKdom'))

#################################### Running DMU #####################################
#Create DMU dir file
filename=paste0('inbredpop4Mod101','rep',iRepScheme,'.DIR')
dmufile=substring(filename, 1, last = (nchar(filename)-4))
write(paste('$COMMENT ADAM-DMU interface'),file=filename,append=F)
write(paste('y = mean + hys + ms_a + NR_a + Res_a + msNRRes_aa + msNRRes_d  + e'),file=filename,append=T)
write(paste('$ANALYSE 1 31 0 0'),file=filename,append=T)
write(paste('$DATA ASCII (9,',nTrait+1,',-9998.0)',paste0('dmudatInbredPop',4,'rep',iRepScheme,'.obs')),file=filename,append=T)
write(paste('$VARIABLE'),file=filename,append=T)
write(paste('#1  2   3    4      5     6        7        8     9    '),file=filename,append=T)
write(paste('id par1 par2 par21 par22 birthpop cupop fixFac1 fixFac2'),file=filename,append=T)
write(paste('# 1     2      3   4   5  '),file=filename,append=T)
write(paste(paste0(paste0('obs',c(1:nTrait)),collapse=' '),'rFac'),file=filename,append=T)
write(paste('$MODEL'),file=filename,append=T)
write(paste('1 1 0 0 0'),file=filename,append=T)
write(paste('0'),file=filename,append=T)
write(paste(nTrait,'0 7 7 8 4  5  2 1 1'),file=filename,append=T)
write(paste('5         1+ 1  2 3 4'),file=filename,append=T)
write(paste('0'),file=filename,append=T)
write(paste('0'),file=filename,append=T)
write(paste('$VAR_STR 1 GREL ASCII', paste0('gmatinvInbredPop',2,'.add')),file=filename,append=T)
write(paste('$VAR_STR 2 GREL ASCII', paste0('gmatinvInbredPop',1,'.add')),file=filename,append=T)
write(paste('$VAR_STR 3 GREL ASCII', paste0('gmatinvInbredPop4threewaymean.aa')),file=filename,append=T)
write(paste('$VAR_STR 4 GREL ASCII', paste0('gmatinvInbredPop',4,'threeway.PSKdom')),file=filename,append=T)
write(paste('$PRIOR'),file=filename,append=T)
write(paste('1  1  1',covarTraits[nTrait,nTrait]/8),file=filename,append=T)
write(paste('2  1  1',covarTraits[nTrait,nTrait]/2),file=filename,append=T)
write(paste('3  1  1',episcovarTraits[nTrait,nTrait]/4),file=filename,append=T)
write(paste('4  1  1',domcovarTraits[nTrait,nTrait]/4),file=filename,append=T)
write(paste('5  1  1',varResiduals[nTrait,nTrait]),file=filename,append=T)
write(paste('$DMUAI'),file=filename,append=T)
write(paste('10     Emstep       Number of steps before full weight on EM in imet = 2'),file=filename,append=T)
write(paste('1.0d-6 Conv_ndelta  Convergence criteria for norm of the update vector '),file=filename,append=T)
write(paste('1.0d-5 Conv_gnorm   Convergence criteria for norm of the gradient vector (AI)'),file=filename,append=T)
write(paste('0      Printout     1 -> Solution vector is printed/written to file SOL'),file=filename,append=T)
write(paste('0      Fspopt	     Time (0) or memory (1) optimised FSPAK'),file=filename,append=T)
write(paste('0      P_neval	     Restart an analysis from evaluation p_neval.'),file=filename,append=T)
dmuRun <- c(paste0('export jobname=',dmufile))
dmuRun <- c(dmuRun,paste('/usr/home/qgg/chuthinh/adam/rdmuai.bsh'))
write(dmuRun, file=paste0("dmuRun"))
system("bash dmuRun")
if (!file.exists('MODINF')) {	write(paste("DMU seems to have failed for",filename),file=outlog,append=T)	}

# clear junks
clean.junk.files()
temp1 <- c(paste0('gmatinvInbredPop',2,'.add'),paste0('gmatinvInbredPop',1,'.add'),
           paste0('gmatinvInbredPop4threewaymean.aa'),paste0('gmatinvInbredPop',4,'threeway.PSKdom'),
           'SOL',paste0(dmufile,'.DMU1.log'))
for (itemp in temp1) { 
  if (file.exists(itemp)) file.remove(itemp) } 

## reading DMU variance output
temp1 <- read.table(paste0(dmufile,'.PAROUT'),header = F)
filename='varDMUparout1.res'
if (!file.exists(filename)) write(paste('iRep',paste0(paste0('val',c(1:length(temp1[,4]))),collapse =' ')),file=filename,append=F)
write(paste(iRepScheme,paste0(temp1[,4],collapse =' ')),file=filename,append=file.exists(filename))


################################### Different GRM ###################################
## gmatrices in statistical models
# ms NR additive
idLst <- pop[pop$fixedFactor1>0 & pop$currentpop %in% c(2),'id']
#idFreq <- pop[pop$fixedFactor1==1 & pop$currentpop %in% c(2),'id']
gmatTemp <- calculateGmatbyRmths1(genodf=geno,idtobeused=idLst,igmatrixfile=paste0('gmatinvInbredPop',2,'.add'),
                                  methodgmat='VanRaden1.add',scaleCov=1,returngmat='none',scaleDiag=FALSE)
# Restorer
idLst <- pop[pop$fixedFactor1>0 & pop$currentpop %in% c(1),'id']
#idFreq <- pop[pop$fixedFactor1==1 & pop$currentpop %in% c(1),'id']
gmatTemp <- calculateGmatbyRmths1(genodf=geno,idtobeused=idLst,igmatrixfile=paste0('gmatinvInbredPop',1,'.add'),
                                  methodgmat='VanRaden1.add',scaleCov=1,returngmat='none',scaleDiag=FALSE)

# 3way hybrid ms NR Res: Epis
pedTemp <- pheno[,c('id','par1','par2','par21','par22')]
genoMean3way=(haplo_pat[pedTemp$par1,]+haplo_mat[pedTemp$par1,]-2)/2 +
  ((haplo_pat[pedTemp$par21,]+haplo_mat[pedTemp$par21,]-2)/2 + 
     (haplo_pat[pedTemp$par22,]+haplo_mat[pedTemp$par22,]-2)/2)/2
genoMean3way=cbind(pedTemp[,'id'],genoMean3way)
gmatTemp <- calculateGmatbyRmths1(genodf=genoMean3way,methodgmat='VanRaden1.add',scaleCov=1,returngmat='gmat',scaleDiag=FALSE)
calculateEpisGmatbyRmths1(gmat1=gmatTemp,igmatrixfile=paste0('gmatinvInbredPop4threewaymean.aa'),scaleDiag=FALSE)
gmatTemp <- NULL

# 3way hybrid ms NR Res: dominance by PSK
calculatedominancegmat3wayPSK(genodf=geno,pedigree=pedTemp,igmatrixfile=paste0('gmatinvInbredPop',4,'threeway.PSKdom'),scalemethod='frequencybased')

#################################### Running DMU #####################################
#Create DMU dir file
filename=paste0('inbredpop4Mod102','rep',iRepScheme,'.DIR')
dmufile=substring(filename, 1, last = (nchar(filename)-4))
write(paste('$COMMENT ADAM-DMU interface'),file=filename,append=F)
write(paste('y = mean + hys + ms_a + NR_a + Res_a + msNRRes_aa + msNRRes_d  + e'),file=filename,append=T)
write(paste('$ANALYSE 1 31 0 0'),file=filename,append=T)
write(paste('$DATA ASCII (9,',nTrait+1,',-9998.0)',paste0('dmudatInbredPop',4,'rep',iRepScheme,'.obs')),file=filename,append=T)
write(paste('$VARIABLE'),file=filename,append=T)
write(paste('#1  2   3    4      5     6        7        8     9    '),file=filename,append=T)
write(paste('id par1 par2 par21 par22 birthpop cupop fixFac1 fixFac2'),file=filename,append=T)
write(paste('# 1     2      3   4   5  '),file=filename,append=T)
write(paste(paste0(paste0('obs',c(1:nTrait)),collapse=' '),'rFac'),file=filename,append=T)
write(paste('$MODEL'),file=filename,append=T)
write(paste('1 1 0 0 0'),file=filename,append=T)
write(paste('0'),file=filename,append=T)
write(paste(nTrait,'0 7 7 8 4  5  2 1 1'),file=filename,append=T)
write(paste('5         1+ 1  2 3 4'),file=filename,append=T)
write(paste('0'),file=filename,append=T)
write(paste('0'),file=filename,append=T)
write(paste('$VAR_STR 1 GREL ASCII', paste0('gmatinvInbredPop',2,'.add')),file=filename,append=T)
write(paste('$VAR_STR 2 GREL ASCII', paste0('gmatinvInbredPop',1,'.add')),file=filename,append=T)
write(paste('$VAR_STR 3 GREL ASCII', paste0('gmatinvInbredPop4threewaymean.aa')),file=filename,append=T)
write(paste('$VAR_STR 4 GREL ASCII', paste0('gmatinvInbredPop',4,'threeway.PSKdom')),file=filename,append=T)
write(paste('$PRIOR'),file=filename,append=T)
write(paste('1  1  1',covarTraits[nTrait,nTrait]/8),file=filename,append=T)
write(paste('2  1  1',covarTraits[nTrait,nTrait]/2),file=filename,append=T)
write(paste('3  1  1',episcovarTraits[nTrait,nTrait]/4),file=filename,append=T)
write(paste('4  1  1',domcovarTraits[nTrait,nTrait]/4),file=filename,append=T)
write(paste('5  1  1',varResiduals[nTrait,nTrait]),file=filename,append=T)
write(paste('$DMUAI'),file=filename,append=T)
write(paste('10     Emstep       Number of steps before full weight on EM in imet = 2'),file=filename,append=T)
write(paste('1.0d-6 Conv_ndelta  Convergence criteria for norm of the update vector '),file=filename,append=T)
write(paste('1.0d-5 Conv_gnorm   Convergence criteria for norm of the gradient vector (AI)'),file=filename,append=T)
write(paste('0      Printout     1 -> Solution vector is printed/written to file SOL'),file=filename,append=T)
write(paste('0      Fspopt	     Time (0) or memory (1) optimised FSPAK'),file=filename,append=T)
write(paste('0      P_neval	     Restart an analysis from evaluation p_neval.'),file=filename,append=T)
dmuRun <- c(paste0('export jobname=',dmufile))
dmuRun <- c(dmuRun,paste('/usr/home/qgg/chuthinh/adam/rdmuai.bsh'))
write(dmuRun, file=paste0("dmuRun"))
system("bash dmuRun")
if (!file.exists('MODINF')) {	write(paste("DMU seems to have failed for",filename),file=outlog,append=T)	}

# clear junks
clean.junk.files()
temp1 <- c(paste0('gmatinvInbredPop',2,'.add'),paste0('gmatinvInbredPop',1,'.add'),
           paste0('gmatinvInbredPop4threewaymean.aa'),paste0('gmatinvInbredPop',4,'threeway.PSKdom'),
           'SOL',paste0(dmufile,'.DMU1.log'))
for (itemp in temp1) { 
  if (file.exists(itemp)) file.remove(itemp) } 

## reading DMU variance output
temp1 <- read.table(paste0(dmufile,'.PAROUT'),header = F)
filename='varDMUparout2.res'
if (!file.exists(filename)) write(paste('iRep',paste0(paste0('val',c(1:length(temp1[,4]))),collapse =' ')),file=filename,append=F)
write(paste(iRepScheme,paste0(temp1[,4],collapse =' ')),file=filename,append=file.exists(filename))

################################### Different GRM ###################################
# 3way hybrid ms NR Res: Epis
pedTemp <- pheno[,c('id','par1','par2','par21','par22')]
genoMean3way=(haplo_pat[pedTemp$par1,]+haplo_mat[pedTemp$par1,]-2)/2 +
  ((haplo_pat[pedTemp$par21,]+haplo_mat[pedTemp$par21,]-2)/2 + 
     (haplo_pat[pedTemp$par22,]+haplo_mat[pedTemp$par22,]-2)/2)/2
genoMean3way=cbind(pedTemp[,'id'],genoMean3way)
gmatTemp <- calculateGmatbyRmths1(genodf=genoMean3way,methodgmat='VanRaden1.add',scaleCov=1,returngmat='gmat',scaleDiag=TRUE,
                                  igmatrixfile=paste0('gmatinvInbredPop4threewaymean.add'))
calculateEpisGmatbyRmths1(gmat1=gmatTemp,igmatrixfile=paste0('gmatinvInbredPop4threewaymean.aa'),scaleDiag=TRUE)
gmatTemp <- NULL

# 3way hybrid ms NR Res: dominance by PSK
calculatedominancegmat3wayPSK(genodf=geno,pedigree=pedTemp,igmatrixfile=paste0('gmatinvInbredPop',4,'threeway.PSKdom'))

#################################### Running DMU #####################################
#Create DMU dir file
filename=paste0('inbredpop4Mod105','rep',iRepScheme,'.DIR')
dmufile=substring(filename, 1, last = (nchar(filename)-4))
write(paste('$COMMENT ADAM-DMU interface'),file=filename,append=F)
write(paste('y = mean + hys + ms_a + NR_a + Res_a + msNRRes_aa + msNRRes_d  + e'),file=filename,append=T)
write(paste('$ANALYSE 1 31 0 0'),file=filename,append=T)
write(paste('$DATA ASCII (9,',nTrait+1,',-9998.0)',paste0('dmudatInbredPop',4,'rep',iRepScheme,'.obs')),file=filename,append=T)
write(paste('$VARIABLE'),file=filename,append=T)
write(paste('#1  2   3    4      5     6        7        8     9    '),file=filename,append=T)
write(paste('id par1 par2 par21 par22 birthpop cupop fixFac1 fixFac2'),file=filename,append=T)
write(paste('# 1     2      3   4   5  '),file=filename,append=T)
write(paste(paste0(paste0('obs',c(1:nTrait)),collapse=' '),'rFac'),file=filename,append=T)
write(paste('$MODEL'),file=filename,append=T)
write(paste('1 1 0 0 0'),file=filename,append=T)
write(paste('0'),file=filename,append=T)
write(paste(nTrait,'0 5 7 8 1 1 1'),file=filename,append=T)
write(paste('3         1 2 3'),file=filename,append=T)
write(paste('0'),file=filename,append=T)
write(paste('0'),file=filename,append=T)
write(paste('$VAR_STR 1 GREL ASCII', paste0('gmatinvInbredPop4threewaymean.add')),file=filename,append=T)
write(paste('$VAR_STR 2 GREL ASCII', paste0('gmatinvInbredPop4threewaymean.aa')),file=filename,append=T)
write(paste('$VAR_STR 3 GREL ASCII', paste0('gmatinvInbredPop',4,'threeway.PSKdom')),file=filename,append=T)
write(paste('$PRIOR'),file=filename,append=T)
write(paste('1  1  1',covarTraits[nTrait,nTrait]),file=filename,append=T)
write(paste('2  1  1',episcovarTraits[nTrait,nTrait]),file=filename,append=T)
write(paste('3  1  1',domcovarTraits[nTrait,nTrait]),file=filename,append=T)
write(paste('4  1  1',varResiduals[nTrait,nTrait]),file=filename,append=T)
write(paste('$DMUAI'),file=filename,append=T)
write(paste('10     Emstep       Number of steps before full weight on EM in imet = 2'),file=filename,append=T)
write(paste('1.0d-6 Conv_ndelta  Convergence criteria for norm of the update vector '),file=filename,append=T)
write(paste('1.0d-5 Conv_gnorm   Convergence criteria for norm of the gradient vector (AI)'),file=filename,append=T)
write(paste('0      Printout     1 -> Solution vector is printed/written to file SOL'),file=filename,append=T)
write(paste('0      Fspopt	     Time (0) or memory (1) optimised FSPAK'),file=filename,append=T)
write(paste('0      P_neval	     Restart an analysis from evaluation p_neval.'),file=filename,append=T)
dmuRun <- c(paste0('export jobname=',dmufile))
dmuRun <- c(dmuRun,paste('/usr/home/qgg/chuthinh/adam/rdmuai.bsh'))
write(dmuRun, file=paste0("dmuRun"))
system("bash dmuRun")
if (!file.exists('MODINF')) {	write(paste("DMU seems to have failed for",filename),file=outlog,append=T)	}

# clear junks
clean.junk.files()
temp1 <- c(paste0('gmatinvInbredPop4threewaymean.add'),paste0('gmatinvInbredPop4threewaymean.aa'),paste0('gmatinvInbredPop',4,'threeway.PSKdom'),
           'SOL',paste0(dmufile,'.DMU1.log'))
for (itemp in temp1) { 
  if (file.exists(itemp)) file.remove(itemp) } 

## reading DMU variance output
temp1 <- read.table(paste0(dmufile,'.PAROUT'),header = F)
filename='varDMUparout5.res'
if (!file.exists(filename)) write(paste('iRep',paste0(paste0('val',c(1:length(temp1[,4]))),collapse =' ')),file=filename,append=F)
write(paste(iRepScheme,paste0(temp1[,4],collapse =' ')),file=filename,append=file.exists(filename))
