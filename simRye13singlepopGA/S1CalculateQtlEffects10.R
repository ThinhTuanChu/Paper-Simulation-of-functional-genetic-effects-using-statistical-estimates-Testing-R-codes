######### Starting population, calculate genetic effects ############
#####################################################################
# Simulate functional effects, or
# Simulate statistical effects, then transform them into functional effects: 2 different methods 

#list of library packages required
temp1 <- c('MASS') 			# c("readxl","plyr","ggplot2",'GA')
temp2 <- temp1 %in% rownames(installed.packages()) # Install packages not yet installed
if (any(temp2 == FALSE)) {
  install.packages(temp1[!temp2])}
#load packages
invisible(lapply(temp1, library, character.only = TRUE))

##### Functional matrices: y = ak + al + dk + dl +aakl; or plot y = a1k + a1l + a2k + a2l + dk + dl +aakl
setupWfklmat1 <- function (tx1=c(0,1,2),tx2=NULL) { # tx1 = c(0,1,2);tx2=NULL  #Any number within 0:2
  if (any(tx1 <0) | any(tx1 >2)) stop ('Wrong genotype codes used. The code should be within 0,1,2.')
  if (is.null(tx2)) {  # Genotype of indivdiuals
    nGCode=length(tx1)
    ta = matrix(tx1-1,nrow=nGCode,ncol=1)
    tx1[tx1>1]=abs(tx1[tx1>1]-2)
    td = matrix(tx1,nrow=nGCode,ncol=1) 
    Ivec = matrix(1,nrow=nGCode,ncol=1)
    Wfkl <- matrix(nrow=nGCode*nGCode,ncol=6)  # 6 here is 6 elements: u ak al dk dl aakl
    Wfkl[,1] <- 1.0
    Wfkl[,2] <- kronecker(Ivec,ta)
    Wfkl[,3] <- kronecker(ta,Ivec)
    Wfkl[,4] <- kronecker(Ivec,td)
    Wfkl[,5] <- kronecker(td,Ivec)
    Wfkl[,6] <- Wfkl[,2]*Wfkl[,3]
    
  } else {       # plot; Genotypes of parents; y = a1k + a1l + a2k + a2l + dk + dl +aakl
    #tx1=c(0,1,2);tx2=c(0,0.5,1,1.5,2)
    nGCode1=length(tx1)
    nGCode2=length(tx2)
    
    ta1 = matrix(kronecker(tx1-1,rep(1,length(tx2))),nrow=nGCode1*nGCode2,ncol=1)/2
    ta2 = matrix(kronecker(rep(1,length(tx1)),tx2-1),nrow=nGCode1*nGCode2,ncol=1)/2
    temp1 <- kronecker(tx1,rep(1,length(tx2)))
    temp2 <- kronecker(rep(1,length(tx1)),tx2)
    td = matrix((1-temp1/2)*(temp2/2) + (temp1/2)*(1-temp2/2),nrow=nGCode1*nGCode2,ncol=1)
    # if scale temp1 was -1,0,1 above equivalent to: 1/4*(1-temp1)*(temp2+1)+1/4*(temp1+1)*(1-temp2)
    #temp0 <- data.frame(ta1,ta2,temp1,temp2,td)
    Ivec = matrix(1,nrow=nGCode1*nGCode2,ncol=1)
    
    Wfkl <- matrix(nrow=nGCode1*nGCode2*nGCode1*nGCode2,ncol=8)  # 8 elements: u a1k a1l  a2k a2l dk dl aakl
    Wfkl[,1] <- 1.0
    Wfkl[,2] <- kronecker(Ivec,ta1)
    Wfkl[,3] <- kronecker(ta1,Ivec)
    Wfkl[,4] <- kronecker(Ivec,ta2)
    Wfkl[,5] <- kronecker(ta2,Ivec)
    Wfkl[,6] <- kronecker(Ivec,td)
    Wfkl[,7] <- kronecker(td,Ivec)
    temp1 = (Wfkl[,2:3] + Wfkl[,4:5])         # xAx=(aRes + aNRMS) for k and l # note: aRes was halfed up there. So no devided by 2
    Wfkl[,8] <- temp1[,1]*temp1[,2]           # xAk * xAl
  }
  return(Wfkl)
}

setupWfklmat2 <- function (tx1=c(0,1,2),p1_k,p1_l,tx2=NULL,p2_k=NULL,p2_l=NULL) { # tx1 = c(0,1,2);tx2=NULL  #Any number within 0:2
  if (!any(tx1 %in% c(0,1,2))) stop ('Wrong genotype codes used. The code should be 0,1,2.')
  if (is.null(tx2)) {  # Genotype of indivdiuals
    nGCode=length(tx1)
    ta_k = matrix((tx1-p1_k[2]-2.0*p1_k[3]),nrow=nGCode,ncol=1)
    ta_l = matrix((tx1-p1_l[2]-2.0*p1_l[3]),nrow=nGCode,ncol=1)
    tx1[tx1>1]=abs(tx1[tx1>1]-2)
    td_k = matrix(tx1-p1_k[2],nrow=nGCode,ncol=1) 
    td_l = matrix(tx1-p1_l[2],nrow=nGCode,ncol=1)
    Ivec = matrix(1,nrow=nGCode,ncol=1)
    Wfkl <- matrix(nrow=nGCode*nGCode,ncol=6)  # 6 here is 6 elements: u ak al dk dl aakl
    Wfkl[,1] <- 1.0
    Wfkl[,2] <- kronecker(Ivec,ta_k)
    Wfkl[,3] <- kronecker(ta_l,Ivec)
    Wfkl[,4] <- kronecker(Ivec,td_k)
    Wfkl[,5] <- kronecker(td_l,Ivec)
    Wfkl[,6] <- Wfkl[,2]*Wfkl[,3]
    
  } else {       # plot; Genotypes of parents; y = a1k + a1l + a2k + a2l + dk + dl +aakl
    if (!any(tx2 %in% c(0,1,2))) stop ('Wrong genotype codes used. The code should be 0,1,2.')
    #tx1=c(0,1,2);tx2=c(0,1,2)
    nGCode1=length(tx1)
    nGCode2=length(tx2)
    
    ta1_k = matrix(kronecker((tx1-p1_k[2]-2.0*p1_k[3]),rep(1,length(tx2))),nrow=nGCode1*nGCode2,ncol=1)/2
    ta1_l = matrix(kronecker((tx1-p1_l[2]-2.0*p1_l[3]),rep(1,length(tx2))),nrow=nGCode1*nGCode2,ncol=1)/2
    
    ta2_k = matrix(kronecker(rep(1,length(tx1)),(tx2-p2_k[2]-2.0*p2_k[3])),nrow=nGCode1*nGCode2,ncol=1)/2
    ta2_l = matrix(kronecker(rep(1,length(tx1)),(tx2-p2_l[2]-2.0*p2_l[3])),nrow=nGCode1*nGCode2,ncol=1)/2
    
    p3_k=((p1_k[1]+p1_k[2]/2)+(p2_k[1]+p2_k[2]/2))/2
    p3_l=((p1_l[1]+p1_l[2]/2)+(p2_l[1]+p2_l[2]/2))/2
    temp1 <- kronecker(tx1,rep(1,length(tx2)))
    temp2 <- kronecker(rep(1,length(tx1)),tx2)
    td_k = matrix((1-temp1/2)*(temp2/2) + (temp1/2)*(1-temp2/2)-2*p3_k*(1-p3_k),nrow=nGCode1*nGCode2,ncol=1)
    td_l = matrix((1-temp1/2)*(temp2/2) + (temp1/2)*(1-temp2/2)-2*p3_l*(1-p3_l),nrow=nGCode1*nGCode2,ncol=1)
    # if scale temp1 was -1,0,1 above equivalent to: 1/4*(1-temp1)*(temp2+1)+1/4*(temp1+1)*(1-temp2)
    #temp0 <- data.frame(ta1,ta2,temp1,temp2,td)
    Ivec = matrix(1,nrow=nGCode1*nGCode2,ncol=1)
    
    Wfkl <- matrix(nrow=nGCode1*nGCode2*nGCode1*nGCode2,ncol=8)  # 8 elements: u a1k a1l  a2k a2l dk dl aakl
    Wfkl[,1] <- 1.0
    Wfkl[,2] <- kronecker(Ivec,ta1_k)
    Wfkl[,3] <- kronecker(ta1_l,Ivec)
    Wfkl[,4] <- kronecker(Ivec,ta2_k)
    Wfkl[,5] <- kronecker(ta2_l,Ivec)
    Wfkl[,6] <- kronecker(Ivec,td_k)
    Wfkl[,7] <- kronecker(td_l,Ivec)
    temp1 = (Wfkl[,2:3] + Wfkl[,4:5])         # xAx=(aRes + aNRMS) for k and l # note: aRes was halfed up there. So no devided by 2
    Wfkl[,8] <- temp1[,1]*temp1[,2]           # xAk * xAl
  }
  return(Wfkl)
}

##### Functions for statistical matrices for indivdiuals
# locus k & l
# function to set up Wx matrix
setupWxsmat <- function(px) {
  Wx <- matrix(nrow=3,ncol=3)
  Wx[,1]=1.0
  Wx[,2]=c((0.0-px[2]-2.0*px[3]),
           (1.0-px[2]-2.0*px[3]),
           (2.0-px[2]-2.0*px[3]))
  #Wx[,3]=0  # no dominance simulated. So, using 0 instead of below column.
  temp1 <- (px[3]+px[1]-((px[1]-px[3])^2))
  if (temp1!=0) {
    Wx[,3]=c((-2.0*px[2]*px[3])/temp1,
             ( 4.0*px[1]*px[3])/temp1,
             (-2.0*px[1]*px[2])/temp1)
  } else {Wx[,3]=0}
  # tested: Duenk 2020 may be correct using his formula. But in this code -2.0 is correct.
  return(Wx)
}
# function to set up Wkl matrix
setupWsklmat <- function(pk,pl) {
  Wskl <- matrix(nrow=3*3,ncol=6)
  Wk=setupWxsmat(pk)
  Wl=setupWxsmat(pl)
  # mean
  # Wskl[,1]=1.0
  # a
  temp1=kronecker(Wl[,c(1,2)],Wk[,c(1,2)]) 
  Wskl[,1:3]=temp1[,1:3]
  # d
  temp1=kronecker(Wl[,c(1,3)],Wk[,c(1,3)]) 
  Wskl[,4:5]=temp1[,2:3]
  # aa 
  Wskl[,6]=kronecker(Wl[,c(2)],Wk[,c(2)]) 
  return(Wskl)
}

##### Functions for statistical matrices for plots
# function to set up Wkl matrix
setupWsklPlotmat <- function(tx1=c(0,1,2),ap1_k,ap2_k,tx2=c(0,0.5,1,1.5,2),ap1_l,ap2_l) {
  #tx1=c(0,1,2);tx2=c(0,0.5,1,1.5,2)
  nGCode1=length(tx1)
  nGCode2=length(tx2)
  # a1
  ta1_k = matrix(kronecker((tx1-2*(1-ap1_k)),rep(1,length(tx2))),nrow=nGCode1*nGCode2,ncol=1)/2
  ta1_l = matrix(kronecker((tx1-2*(1-ap1_l)),rep(1,length(tx2))),nrow=nGCode1*nGCode2,ncol=1)/2
  # a2
  ta2_k = matrix(kronecker(rep(1,length(tx1)),(tx2-2*(1-ap2_k))),nrow=nGCode1*nGCode2,ncol=1)/2
  ta2_l = matrix(kronecker(rep(1,length(tx1)),(tx2-2*(1-ap2_l))),nrow=nGCode1*nGCode2,ncol=1)/2
  
  #d
  temp1 <- kronecker(tx1,rep(1,length(tx2)))
  pa11=1-temp1/2; pa12=temp1/2
  temp2 <- kronecker(rep(1,length(tx1)),tx2)
  pa21=1-temp2/2; pa22=temp2/2 
  # if scale temp1 was -1,0,1 above equivalent to: 1/4*(1-temp1)*(temp2+1)+1/4*(temp1+1)*(1-temp2)
  temp0=pa11*pa22*2*(1-ap1_k)*ap2_k + pa12*pa21*2*ap1_k*(1-ap2_k) +
    pa11*pa21*(-2)*(1-ap1_k)*(1-ap2_k) + pa12*pa22*(-2)*ap1_k*ap2_k
  td_k = matrix(temp0,nrow=nGCode1*nGCode2,ncol=1)
  #temp0 <- data.frame(ta1_k,ta2_k,temp1,temp2,td_k)
  
  temp0=pa11*pa22*2*(1-ap1_l)*ap2_l + pa12*pa21*2*ap1_l*(1-ap2_l) +
    pa11*pa21*(-2)*(1-ap1_l)*(1-ap2_l) + pa12*pa22*(-2)*ap1_l*ap2_l
  td_l = matrix(temp0,nrow=nGCode1*nGCode2,ncol=1)
  #temp0 <- data.frame(ta1_l,ta2_l,temp1,temp2,td_l)
  
  Ivec = matrix(1,nrow=nGCode1*nGCode2,ncol=1)
  
  Wskl <- matrix(nrow=nGCode1*nGCode2*nGCode1*nGCode2,ncol=8)  # 8 elements: u a1k a1l  a2k a2l dk dl aakl
  Wskl[,1] <- 1.0
  Wskl[,2] <- kronecker(Ivec,ta1_k)
  Wskl[,3] <- kronecker(ta1_l,Ivec)
  Wskl[,4] <- kronecker(Ivec,ta2_k)
  Wskl[,5] <- kronecker(ta2_l,Ivec)
  Wskl[,6] <- kronecker(Ivec,td_k)
  Wskl[,7] <- kronecker(td_l,Ivec)
  temp1 = (Wskl[,2:3] + Wskl[,4:5])         # xAx=(aRes + aNRMS) for k and l # note: aRes was halfed up there. So no devided by 2
  Wskl[,8] <- temp1[,1]*temp1[,2]           # xAk * xAl
  
  return(Wskl)
}

# set up hai and hdi colum vector
findhai <- function(px,allelecount) {
  # if (!allelecount %in% c(0,1,2)) stop('error: wrong allele count input.')  # comment out to improve speed
  hai=(allelecount-px[2]-2.0*px[3]);  return(hai)
}
findhdi <- function(px,allelecount) {
  # if (!allelecount %in% c(0,1,2)) stop('error: wrong allele count input.')  # comment out to improve speed
  temp1 <- (px[3]+px[1]-((px[1]-px[3])^2))
  if (temp1!=0) {
    if (allelecount==0) {
      hdi=(-2.0*px[2]*px[3])/temp1; return(hdi)	}
    if (allelecount==1) {
      hdi=( 4.0*px[1]*px[3])/temp1; return(hdi)	}
    if (allelecount==2) {
      hdi=(-2.0*px[1]*px[2])/temp1; return(hdi)	}
  } else {
    hdi=0; return(hdi)
  }
}

findhdiPlot <- function(p1,p2,allelecount1,allelecount2) {  
  # allele counts: (0:2)
  pa11=1-allelecount1/2; pa12=allelecount1/2
  pa21=1-allelecount2/2; pa22=allelecount2/2 
  #if scale allelecount1 was -1,0,1:   pa11=1/2(1-a1);   pa12=1/2(a1+1)
  #if scale allelecount2 was -1,0,1:   pa21=1/2(1-a2);   pa22=1/2(a2+1)
  ## proBb = pa11*pa22 + pa12*pa21
  ## proBB = pa11*pa21
  ## probb = pa12*pa22
  hdi=pa11*pa22*2*(1-p1)*p2 + pa12*pa21*2*p1*(1-p2) +
    pa11*pa21*(-2)*(1-p1)*(1-p2) + pa12*pa22*(-2)*p1*p2
  return(hdi)
}

# set up kronecker(hai,hai)
findhaihai <- function(pk,ack,pl,acl) {
  # if ((!ack %in% c(0:2))|(!acl %in% c(0:2))) stop('error: wrong allele count input.')  # comment out to improve speed
  return(((ack-pk[2]-2.0*pk[3])) * ((acl-pl[2]-2.0*pl[3])))
}

# set up kronecker(hai,hai)
findhaihaiPlot <- function(ap1_k,ap2_k,ac1_k,ac2_k,    ap1_l,ap2_l,ac1_l,ac2_l) {
  hai_k = (ac1_k-2*(1-ap1_k))/2 + (ac2_k-2*(1-ap2_k))/2
  hai_l = (ac1_l-2*(1-ap1_l))/2 + (ac2_l-2*(1-ap2_l))/2
  return(hai_k*hai_l)
}

calculateVarbyLocus <- function(EffInput,px,hx=c(-1,0,1)) {  # if hx=c(-1,0,1), biological additive variances
  # nTrait=2;EffInput=c(100,8888);px=c(0.4,0.1,0.5);hx=c(-1,0,1)
  #if (length(px) != length(hx)) stop('error: wrong inputs. length of hx not equivalent to frequency.')  # comment out to improve speed
  #if (length(EffInput) != nTrait) stop('error: wrong inputs. Effects are not equal to number of traits.')     # comment out to improve speed
  Eff = array(data = 0,dim=c(length(hx),nTrait))
  Eff[,] =  kronecker(EffInput,hx)
  varLocusj = array(data = NA,dim=c(nTrait,nTrait))
  for (iTrait in 1:(nTrait)) {
    meanpopiTrait=sum(Eff[,iTrait]*px)
    for (jTrait in 1:(nTrait)) {
      meanpopjTrait=sum(Eff[,jTrait]*px)
      varLocusj[iTrait,jTrait]=sum(((Eff[,iTrait]-meanpopiTrait)*(Eff[,jTrait]-meanpopjTrait))*px)
    }	 
  } 
  return(varLocusj)
}

calculateVarbyLocusEpis <- function(EffInput,pk,pl,hk=c(-1,0,1),hl=c(-1,0,1)) {   # if hx=c(-1,0,1), biological additive variances
  # nTrait=2;EffInput=c(100,8888);pk=c(0.4,0.1,0.5);pl=c(0.2,0.7,0.1);hk=c(-1,0,1);hl=c(-1,0,1)
  #if (length(pk) != length(hk) | length(pl) != length(hl) ) stop('error: wrong inputs. length of hk not equivalent to frequency.')  # comment out to improve speed
  #if (length(EffInput) != nTrait) stop('error: wrong inputs. Effects are not equal to number of traits.')     # comment out to improve speed
  hkhl = kronecker(hl,hk)
  pkpl = kronecker(pl,pk)
  Eff = array(data = 0,dim=c(length(hkhl),nTrait))
  Eff[,] =  kronecker(EffInput,hkhl)
  varLocusj = array(data = NA,dim=c(nTrait,nTrait))
  for (iTrait in 1:(nTrait)) {
    meanpopiTrait=sum(Eff[,iTrait]*pkpl)
    for (jTrait in 1:(nTrait)) {
      meanpopjTrait=sum(Eff[,jTrait]*pkpl)
      varLocusj[iTrait,jTrait]=sum(((Eff[,iTrait]-meanpopiTrait)*(Eff[,jTrait]-meanpopjTrait))*pkpl)
    }	 
  } 
  return(varLocusj)
}

calculateStatVarDPlotbyLocus <- function(EffInput,ap1,ap2,p1,p2,h1=c(0,1,2),h2=c(0,1,2)) {
  #  nTrait=2;EffInput=c(100,8888);ap1=0.4;ap2=0.3;h1=c(0,1,2);h2=c(0,1,2)
  if (length(p1) != length(h1) | length(p2) != length(h2) ) stop('error: wrong inputs. length of hk not equivalent to frequency.')  # comment out to improve speed
  #  pa11=1-h1/2; pa12=h1/2
  #  pa21=1-h2/2; pa22=h2/2 
  pa11=kronecker(rep(1,length(h2)),1-h1/2)
  pa12=kronecker(rep(1,length(h2)),h1/2)
  pa21=kronecker(1-h2/2,rep(1,length(h1)))
  pa22=kronecker(h2/2,rep(1,length(h1)))
  
  h1h2 = pa11*pa22*2*(1-ap1)*ap2 + pa12*pa21*2*ap1*(1-ap2) +
    pa11*pa21*(-2)*(1-ap1)*(1-ap2) + pa12*pa22*(-2)*ap1*ap2
  p1p2 = kronecker(p2,p1)
  
  Eff = array(data = 0,dim=c(length(h1h2),nTrait))
  Eff[,] =  kronecker(EffInput,h1h2)
  varLocusj = array(data = NA,dim=c(nTrait,nTrait))
  for (iTrait in 1:(nTrait)) {
    meanpopiTrait=sum(Eff[,iTrait]*p1p2)
    for (jTrait in 1:(nTrait)) {
      meanpopjTrait=sum(Eff[,jTrait]*p1p2)
      varLocusj[iTrait,jTrait]=sum(((Eff[,iTrait]-meanpopiTrait)*(Eff[,jTrait]-meanpopjTrait))*p1p2)
    }	 
  } 
  return(varLocusj)
}

calculateVarPlotbyLocusEpis <- function(EffInput,apk1,apk2,pk1,pk2,hk1=c(0,1,2),hk2=c(0,1,2),
                                        apl1,apl2,pl1,pl2,hl1=c(0,1,2),hl2=c(0,1,2)) {  # if hx=c(-1,0,1), biological additive x additive variances
  #nTrait=2;EffInput=c(100,8888);apk1=0.4;apk2=0.3;hk1=c(-1,0,1);hk2=c(-1,0,1)
  #pk1=c(apk1^2,2*apk1*(1-apk1),(1-apk1)^2); pk2=c(apk2^2,2*apk2*(1-apk2),(1-apk2)^2);
  #apl1=0.1;apl2=0.7;hl1=c(-1,0,1);hl2=c(-1,0,1)
  #pl1=c(apl1^2,2*apl1*(1-apl1),(1-apl1)^2); pl2=c(apl2^2,2*apl2*(1-apl2),(1-apl2)^2);
  
  hk1 = (hk1-2*(1-apk1))/2
  hk2 = (hk2-2*(1-apk2))/2
  hl1 = (hl1-2*(1-apl1))/2
  hl2 = (hl2-2*(1-apl2))/2
  
  hk = kronecker(rep(1,length(hk2)),hk1)+kronecker(hk2,rep(1,length(hk1)))
  hl = kronecker(rep(1,length(hl2)),hl1)+kronecker(hl2,rep(1,length(hl1)))
  hkl= kronecker(hl,hk)
  
  pk = kronecker(pk2,pk1)
  pl = kronecker(pl2,pl1)
  pkl= kronecker(pl,pk)
  
  Eff = array(data = 0,dim=c(length(hkl),nTrait))
  Eff[,] =  kronecker(EffInput,hkl)
  varLocusj = array(data = NA,dim=c(nTrait,nTrait))
  for (iTrait in 1:(nTrait)) {
    meanpopiTrait=sum(Eff[,iTrait]*pkl)
    for (jTrait in 1:(nTrait)) {
      meanpopjTrait=sum(Eff[,jTrait]*pkl)
      varLocusj[iTrait,jTrait]=sum(((Eff[,iTrait]-meanpopiTrait)*(Eff[,jTrait]-meanpopjTrait))*pkl)
    }	 
  } 
  return(varLocusj)
}

#####################################################################
############# Function(s) to simulate genetic effects ###############
# required global variables: nQtl, nTrait, covarTraits, SNPinfo, episcovarTraits,domcovarTraits, domdegreemean, domdegreevar, nEpisPairLoci, episPairLoci
# prior 
covarTraits=aRes_covar
episcovarTraits=aaNRMSRes_covar
domcovarTraits=dNRMSRes_covar
# genodf=Pops[[2]][seq(1,200,by=2),]+Pops[[2]][seq(2,200,by=2),]-2
# covarTraits=matrix(c(400, 0.00,0.0, 400),nrow = nTrait,byrow = T)     # testing
# episcovarTraits=matrix(c(100, 0.00,0.0, 100),nrow = nTrait,byrow = T) # testing
# domcovarTraits=matrix(c(100, 0.00,0.0, 100),nrow = nTrait,byrow = T)  # testing
# idtobeused=NULL;idfrequencyCalculation=NULL;infuncParameters='biological';reportresult=TRUE;domdegreemean=NULL;domdegreevar=NULL

caculateQtlEffectsTwo <- function(genodf=NULL,idtobeused=NULL,idfrequencyCalculation=NULL,infuncParameters='biological',reportresult=TRUE,
                                  covarTraits=NULL,episcovarTraits=NULL,domcovarTraits=NULL,domdegreemean=NULL,domdegreevar=NULL,
                                  methodCalculateVar='byid') {
  #genodf=NULL;idtobeused=NULL;idfrequencyCalculation=NULL;infuncParameters='biological';reportresult=TRUE;
  #covarTraits=NULL;episcovarTraits=NULL;domcovarTraits=NULL;domdegreemean=NULL;domdegreevar=NULL;methodCalculateVar='byid'
  #methodCalculateVar='byid'# (options: byid or bylocus)
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
  geno <- geno[idtobeused,]
  rownames(geno) <- NULL   # No need rownames
  
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
  
  
  #### Calculate statistical effects given functional effects simulated, apply to any infuncParameters. Recalculated
  
  # setup bkl using ivares-Castro and Calborg2007 method, to calculate functional genetic effects from statistical effects of a, d and e
  # only work if (nQtlEachQtlInteractwith==1)
  if (nQtlEachQtlInteractwith!=1) stop ('this method works with nQtlEachQtlInteractwith=1 only.') 
  
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
    #Wfkl <- setupWfklmat1()
    
    Wskl <- setupWsklmat(pk,pl)
    temp1 <- t(Wskl)%*%Wskl
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
  
  ####################################
  ############ Report ################
  ####################################
  ## additive genetics 
  # biological by id
  idgv <- array(data = 0,dim=c(nAnimal,nTrait)) 
  for (idd in 1:nAnimal){
    genotypecode=geno[idd,SNPinfo[,4]=='T']
    for (iTrait in 1:(nTrait)) {
      idgv[idd,iTrait]=sum(aQtleffects[,iTrait] * (genotypecode-1))
    }
  }
  v_temp=var(idgv)   # calculate covariances for traits
  if (reportresult) print('Additive genetic biological effect. Variance calculated by id. Mean & variance:')
  temp1 <- c()
  for (iTrait in 1:(nTrait)) {	     temp1 <- c(temp1,mean(idgv[,iTrait]))	 }
  if (reportresult) print(temp1)
  if (reportresult) print(v_temp)
  
  # biological by locus
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
  if (reportresult) print('Additive genetic biological effect. Variance calculated by locus:')
  if (reportresult) print(v_temp)
  
  # statistical by id
  idgv <- array(data = 0,dim=c(nAnimal,nTrait)) 
  temp1 <- freqGenotypes[SNPinfo[,4]=='T',]
  px <- lapply(seq_len(nrow(temp1)), function(i) temp1[i,])
  for (idd in 1:nAnimal){
    genotypecode=geno[idd,SNPinfo[,4]=='T']
    haiVect <- mapply(findhai,px=px,allelecount=genotypecode)
    for (iTrait in 1:(nTrait)) {
      idgv[idd,iTrait] <- sum(haiVect*aQtleffects_Stat[,iTrait])
    }	
  }
  v_temp=var(idgv)   # calculate covariances for traits
  if (reportresult) print('Additive genetic statistical effect. Variance calculated by id. Mean & variance:')
  temp1 <- c()
  for (iTrait in 1:(nTrait)) {	     temp1 <- c(temp1,mean(idgv[,iTrait]))	 }
  if (reportresult) print(temp1)
  if (reportresult) print(v_temp)
  
  # statistical by locus
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
  if (reportresult) print('Additive genetic statistical effect. Variance calculated by locus:')
  if (reportresult) print(v_temp)
  
  ## dominance
  # biological by id
  idgv <- array(data = 0,dim=c(nAnimal,nTrait)) 
  for (idd in 1:nAnimal){
    genotypecode=geno[idd,SNPinfo[,4]=='T']
    for (iTrait in 1:(nTrait)) {
      idgv[idd,iTrait]=sum(domValues[genotypecode==1,iTrait])
    }
  }
  v_temp=var(idgv)   # calculate covariances for traits
  if (reportresult) print('Dominance biological effect. Variance calculated by id. Mean & variance:')
  temp1 <- c()
  for (iTrait in 1:(nTrait)) {	     temp1 <- c(temp1,mean(idgv[,iTrait]))	 }
  if (reportresult) print(temp1)
  if (reportresult) print(v_temp)
  
  # biological by locus
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
  if (reportresult) print('Dominance biological effect. Variance calculated by locus:')
  if (reportresult) print(v_temp)
  
  # statistical by id
  idgv <- array(data = 0,dim=c(nAnimal,nTrait)) 
  temp1 <- freqGenotypes[SNPinfo[,4]=='T',]
  px <- lapply(seq_len(nrow(temp1)), function(i) temp1[i,])
  for (idd in 1:nAnimal){
    genotypecode=geno[idd,SNPinfo[,4]=='T']
    hdiVect <- mapply(findhdi,px=px,allelecount=genotypecode)
    for (iTrait in 1:(nTrait)) {
      idgv[idd,iTrait] <- sum(hdiVect*domValues_Stat[,iTrait])
    }	
  }
  v_temp=var(idgv)   # calculate covariances for traits
  if (reportresult) print('Dominance statistical effect. Variance calculated by id. Mean & variance:')
  temp1 <- c()
  for (iTrait in 1:(nTrait)) {	     temp1 <- c(temp1,mean(idgv[,iTrait]))	 }
  if (reportresult) print(temp1)
  if (reportresult) print(v_temp)
  
  # statistical by locus
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
  if (reportresult) print('Dominance statistical effect. Variance calculated by locus:')
  if (reportresult) print(v_temp)
  
  ## Epistasis
  if (sum(diag(episcovarTraits))>0) {
    # biological 1 by id
    idgv <- array(data = 0,dim=c(nAnimal,nTrait))
    for (idd in 1:nAnimal){
      xAk <- geno[idd,episPairLoci[,1]]-1
      xAl <- geno[idd,episPairLoci[,2]]-1
      for (iTrait in 1:(nTrait)) {
        idgv[idd,iTrait]=sum(episPairLociEffect[1:nEpisPairLoci,iTrait]*xAk*xAl)
      }
    }
    v_temp=var(idgv)   # calculate covariances for traits
    if (reportresult) print('Epistasis biological effect. Variance calculated by id. Mean & variance:')
    temp1 <- c()
    for (iTrait in 1:(nTrait)) {	     temp1 <- c(temp1,mean(idgv[,iTrait]))	 }
    if (reportresult) print(temp1)
    if (reportresult) print(v_temp)
    
    # biological by locus
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
    if (reportresult) print('Biological epistatic effect. Variance calculated by locus:')
    if (reportresult) print(v_temp)	
    
    # statistical by id
    idgv <- array(data = 0,dim=c(nAnimal,nTrait))
    temp1 <- freqGenotypes[episPairLoci[,1],]
    pxAk <- lapply(seq_len(nrow(temp1)), function(i) temp1[i,])
    temp1 <- freqGenotypes[episPairLoci[,2],]
    pxAl <- lapply(seq_len(nrow(temp1)), function(i) temp1[i,])
    
    for (idd in 1:nAnimal){
      xAk <- geno[idd,episPairLoci[,1]]
      xAl <- geno[idd,episPairLoci[,2]]
      haihaiVect <- mapply(findhaihai,pk=pxAk,ack=xAk,pl=pxAl,acl=xAl)		
      for (iTrait in 1:(nTrait)) {
        idgv[idd,iTrait]=sum(haihaiVect*episPairLociEffect_Stat[,iTrait])    # episPairLociEffect = episPairLociEffect_Stat
      }
    }
    v_temp=var(idgv)   # calculate covariances for traits
    if (reportresult) print('Statistical epistatic effect. Variance calculated by id. Mean & variance:')
    temp1 <- c()
    for (iTrait in 1:(nTrait)) {	     temp1 <- c(temp1,mean(idgv[,iTrait]))	 }
    if (reportresult) print(temp1)
    if (reportresult) print(v_temp)
    
    # statistical by locus
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
    if (reportresult) print('Statistical epistatic effect. Variance calculated by locus:')
    if (reportresult) print(v_temp)	
  }
  
  # tgv
  idgv <- array(data = 0,dim=c(nAnimal,nTrait))
  for (idd in 1:nAnimal){
    
    genotypecode=geno[idd,SNPinfo[,4]=='T']
    for (iTrait in 1:(nTrait)) {
      idgv[idd,iTrait]=idgv[idd,iTrait]+sum(aQtleffects[,iTrait] * (genotypecode-1))
      idgv[idd,iTrait]=idgv[idd,iTrait]+sum(domValues[genotypecode==1,iTrait])
    }
    if (sum(diag(episcovarTraits))>0) {
      xAk <- geno[idd,episPairLoci[,1]]-1
      xAl <- geno[idd,episPairLoci[,2]]-1
      for (iTrait in 1:(nTrait)) {
        idgv[idd,iTrait]=idgv[idd,iTrait]+sum(episPairLociEffect[1:nEpisPairLoci,iTrait]*xAk*xAl)
      }	  
    }
  }
  v_temp=var(idgv)   # calculate covariances for traits
  if (reportresult) print('Total genetic effect. Mean & variance:')
  temp1 <- c()
  for (iTrait in 1:(nTrait)) {	     temp1 <- c(temp1,mean(idgv[,iTrait]))	 }
  if (reportresult) print(temp1)
  if (reportresult) print(v_temp)
  
  # return
  return(list(aQtleffects,domValues,domdegreeValues,episPairLociEffect))
}

# run calculating effects
temp1=Pops[[1]][seq(1,nFounderAnimal*2,by=2),]+Pops[[1]][seq(2,nFounderAnimal*2,by=2),]-2
temp2=caculateQtlEffectsTwo(genodf=temp1,infuncParameters='biological', reportresult=TRUE,
                            covarTraits=aRes_covar,episcovarTraits=aaNRMSRes_covar,domcovarTraits=dNRMSRes_covar)
#temp2=caculateQtlEffectsTwo(genodf=temp1,infuncParameters='statisticalOne', reportresult=TRUE,
#						covarTraits=aRes_covar,episcovarTraits=aaNRMSRes_covar,domcovarTraits=dNRMSRes_covar)
#temp2=caculateQtlEffectsTwo(genodf=temp1,infuncParameters='statisticalTwo', reportresult=TRUE,
#                            covarTraits=aRes_covar,episcovarTraits=aaNRMSRes_covar,domcovarTraits=dNRMSRes_covar,methodCalculateVar='byid')
#temp2=caculateQtlEffectsTwo(genodf=temp1,infuncParameters='statisticalTwo', reportresult=TRUE,
#                            covarTraits=aRes_covar,episcovarTraits=aaNRMSRes_covar,domcovarTraits=dNRMSRes_covar,methodCalculateVar='bylocus')
#temp2=caculateQtlEffectsTwo(genodf=temp1,infuncParameters='statisticalThree', reportresult=TRUE,
#                            covarTraits=aRes_covar,episcovarTraits=aaNRMSRes_covar,domcovarTraits=dNRMSRes_covar)

aQtleffects=temp2[[1]]
domValues=temp2[[2]]
domdegreeValues=temp2[[3]]
episPairLociEffect=temp2[[4]]

#
keepVariables <- ls()
# remove unneeded variables
#keepVariables <- c(keepVariables,'domValues','domdegreeValues','aQtleffects',
#                   'episPairLociEffect','nEpisPairLoci','episPairLoci')
rm(list=setdiff(ls(), keepVariables))
