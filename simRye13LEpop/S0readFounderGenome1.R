######## read genome structures, calculate genetic effects ##########
#####################################################################

######## Check some inputs and set defaults values ##################
nAlleles=2; idplotobs <- NULL
if (!inputParameters %in% c('biological','statisticalOne','statisticalTwo')) {
   print(paste('Invalid inputParameters provided',inputParameters,'. Thus, biological input paramters are assumed.'))
   inputParameters='biological'
}
if (is.null(domdegreevar))  {   domdegreevar <- diag(rep(sqrt(0.0973),nTrait))  }
if (is.null(domdegreemean)) {   domdegreemean <- rep(0.193,nTrait)  }
keepVariables <- ls()

###### read in geneticArchitecture ######
infile <- readLines(geneticArchitectureFile)
nChr0 <- as.numeric(infile[1])
SNP <- infile[2:(1+nChr0)]
write(SNP,file = 'temp.txt')
SNP <- read.table('temp.txt',header = F)
SNPinfo <- infile[(2+nChr0):length(infile)]
write(SNPinfo,file = 'temp.txt')
SNPinfo <- read.table('temp.txt',header = F)

if (is.null(nChr)) nChr=nChr0
if (is.null(nMarkerEachChr) | (length(nMarkerEachChr)!=nChr)) {
  nMarkerEachChr <- SNP[,2]
}
if (is.null(nQtlEachChr) | (length(nQtlEachChr)!=nChr)) {
  nQtlEachChr <- SNP[,1]
}
# extract genetic Architecture
SNPinfo$snpNo <- c(1:length(SNPinfo[,1]))
SNPinfo <- SNPinfo[SNPinfo[,2]<=nChr,]
SNPinfo[,4][SNPinfo[,4]==F] <- "F"
SNPinfo[,4][SNPinfo[,4]==T] <- "T"

# sample qtl and marker from above
snp_used <- c()
for (i in 1:nChr) {
  temp1  <- SNPinfo[SNPinfo[,2]==i,]
  temp11 <- temp1[temp1[,4]=='T',]
  temp12 <- temp1[temp1[,4]=='F',]
  if (length(temp11$snpNo)<nQtlEachChr[i]) nQtlEachChr[i]=length(temp11$snpNo)
  if (length(temp12$snpNo)<nMarkerEachChr[i]) nMarkerEachChr[i]=length(temp12$snpNo)
  temp21 <- sample(temp11$snpNo,nQtlEachChr[i])
  temp22 <- sample(temp12$snpNo,nMarkerEachChr[i])
  snp_used <- c(snp_used,temp21,temp22)
}
snp_used <- sort(snp_used)
SNP <- data.frame(nQtlEachChr,nMarkerEachChr)
SNPinfo <- SNPinfo[SNPinfo$snpNo %in% snp_used,]
nSnp = length(SNPinfo[,1])
nQtl=length(SNPinfo[,4][SNPinfo[,4]=="T"])
nMarkers=length(SNPinfo[,4][SNPinfo[,4]=="F"])

# change allele frequency # Redistribute SNP # only apply to LEpop
SNPinfo[,1] <- sample(SNPinfo[,1],nSnp)

# chromosome info
chrom <- vector(mode="list", length=nChr)  									# list hold chromosome info
#counter=0
for (ichrom in 1:nChr) {
  chrom[[ichrom]] <- list()
  chrom[[ichrom]]$chromlength= ceiling(max(SNPinfo[SNPinfo[,2]==ichrom,3]))  # Size of chromosome in cM
  chrom[[ichrom]]$nSnpEachChr = nMarkerEachChr[ichrom] + nQtlEachChr[ichrom]
  chrom[[ichrom]]$nMarkerEachChr = nMarkerEachChr[ichrom]
  chrom[[ichrom]]$nQtlEachChr = nQtlEachChr[ichrom]
  chrom[[ichrom]]$qtl = SNPinfo[SNPinfo[,2]==ichrom,4]
}

# read in founder populations

### declare important variables to be used frequently later
Pops <- vector(mode="list", length=nFounderPop+1)  									# list hold genotypes of purebred Populations (last pops is from allpop)
allPop <- matrix(NA,nrow=nFounderAnimal*nFounderPop*2,ncol=nSnp)
freqAllelesFounder <- array(data = 0,dim=c(nSnp,nAlleles,nFounderPop+1))			# 3 dimensions: nSnp,nAlleles,nFounderPop+1 (last pops is from allpop)

########### creat purebred Populations with LD from defined base haplotypes #############
for (iPop in 1:nFounderPop) Pops[[iPop]] <- matrix(NA, nrow = nFounderAnimal*2, ncol= nSnp) 
n=1		#;nBaseHap=200  # nBaseAnimal is number of animals*2 in baseHaplotypes for each basePop (should be equal number)
for (iPop in 1:nFounderPop) {
  
  # start: only apply to LEpop
  Poptemp <- matrix(NA, nrow=nBaseHap,ncol=nSnp)
  for (i in 1:nBaseHap) {
     Poptemp[i,] <- rbinom(nSnp,1,SNPinfo[,1])+1
  }
  # end: only apply to LEpop
  
  if (nFounderAnimal<(nBaseHap/2)) {
    # sample chro
    temp1 <- sample(nBaseHap,nFounderAnimal*2,replace = F)
    Pops[[iPop]][1:(nFounderAnimal*2),] <- Poptemp[temp1,]
  } else {
  for (i in 1:nFounderAnimal) {
    # sample chro
    temp1 <- sample(nBaseHap,2,replace = F)
    Pops[[iPop]][(2*i-1):(2*i),] <- Poptemp[temp1,]
  }
  }
  n=n+1+nBaseHap
}

# calculate frequency of alleles for Pops
MAF <- function(x,y){sum(x==y)/length(x)}
for (iPop in 1:nFounderPop) {
  temp1 <- as.matrix(Pops[[iPop]])
  for (iAllele in 1:nAlleles) {      
     freqAllelesFounder[,iAllele,iPop] <- apply(temp1, 2, MAF,y=iAllele)
  }
}

for (iPop in 1:nFounderPop) {
  allPop[((iPop-1)*(nFounderAnimal*2)+1):((nFounderAnimal*2)*iPop),]=Pops[[iPop]]
}
temp1 <- as.matrix(allPop)
for (iAllele in 1:nAlleles) {
    freqAllelesFounder[,iAllele,nFounderPop+1] <- apply(temp1, 2, MAF,y=iAllele)
}

# frequency of genotypes
freqGenotypeFounder <- array(data = 0,dim=c(nSnp,nAlleles*(nAlleles+1)/2,nFounderPop+1))				# 3 dimensions: nSnp,nAlleles*(nAlleles+1)/2,nFounderPop+1 (last pops is from allpop)
GenoFreq <- function(x,y){sum(x==y)/length(x)}
for (iPop in 1:nFounderPop) {
    temp1 <- as.matrix(Pops[[iPop]])
	temp2 <- temp1[seq(1,length(temp1[,1]),by=2),]+temp1[seq(2,length(temp1[,1]),by=2),]-2
	counter=0
    for (iAllele in 1:nAlleles) {
	  for (jAllele in iAllele:nAlleles) {
          counter=counter+1
		  freqGenotypeFounder[,counter,iPop] <- apply(temp2, 2, GenoFreq,y=(counter-1))
	  }
    }
}
temp1 <- as.matrix(allPop)
temp2 <- temp1[seq(1,length(temp1[,1]),by=2),]+temp1[seq(2,length(temp1[,1]),by=2),]-2
counter=0
for (iAllele in 1:nAlleles) {
	  for (jAllele in iAllele:nAlleles) {
          counter=counter+1
		  freqGenotypeFounder[,counter,nFounderPop+1] <- apply(temp2, 2, GenoFreq,y=(counter-1))
	  }
}

# sample for allpop to Pops with same number of nFounderAnimal
idaddingpop <- c() # ensure same number of individuals from each pop
for (ipop in 1:nFounderPop) {
  temp1 <- ((ipop-1)*nFounderAnimal+1):(ipop*nFounderAnimal)
  if (ipop==nFounderPop) {
      temp3 <- nFounderAnimal-length(idaddingpop)
  } else {temp3 <- floor(nFounderAnimal/nFounderPop)}
  temp2 <- sample(temp1,temp3)  
  idaddingpop <- c(idaddingpop,temp2)
}
idaddingpop <- idaddingpop[order(idaddingpop)]
idused=c()
for (idd in idaddingpop) {
   idused=c(idused,2*idd-1,2*idd)
}
Pops[[nFounderPop+1]]=allPop[idused,]

#important outputs: Pops, Lociw4alleles, allPop, freqAllelesFounder, freqGenotypeFounder
##################################--------****--------##################################

### simulate paired loci interactions 
  # identify number of interactions & sample interaction loci pairs
  nEpisPairLoci=floor(nQtlEachQtlInteractwith*nQtl/2)
  if (nQtlEachQtlInteractwith==1) {
    episPairLoci <- array(data = NA,dim=c(nEpisPairLoci,2))
    temp0 <- which(SNPinfo[,4]=="T")
    temp1 <- sample(temp0,nEpisPairLoci)
    temp2 <- sample(temp0[!temp0 %in% temp1],nEpisPairLoci)
    episPairLoci[,1] <- temp1
    episPairLoci[,2] <- temp2
    
  } else {
    counter=0; loopcondition=TRUE
    while (loopcondition | counter<100) {   # resample if sample loci
      episPairLoci <- array(data = NA,dim=c(nEpisPairLoci,2))
      episPairLoci[,1] <- sample(rep(which(SNPinfo[,4]=="T"),nQtlEachQtlInteractwith),nEpisPairLoci)
      episPairLoci[,2] <- sample(rep(which(SNPinfo[,4]=="T"),nQtlEachQtlInteractwith),nEpisPairLoci)
      loopcondition=any(episPairLoci[,1]==episPairLoci[,2])
      counter=counter+1
    }
    if (loopcondition) {   # resample not converged, accept solution, but exclude non-interaction loci
      nEpisPairLoci=nEpisPairLoci-sum(episPairLoci[,1]==episPairLoci[,2])
      episPairLoci=episPairLoci[!episPairLoci[,1]==episPairLoci[,2],]
    }
    # remove duplicated
    episPairLoci=unique(episPairLoci)
    nEpisPairLoci=dim(episPairLoci)[1]
  }  
  if (nEpisPairLoci<1) stop ('Pair up of episPairLoci is not possible.')
  print(paste('Number of paired loci for simulation of epistatic effects:',nEpisPairLoci))
#important outputs: 'nEpisPairLoci','episPairLoci'
##################################--------****--------##################################

# remove unneeded variables
keepVariables <- c(keepVariables,'Pops','Lociw4alleles','allPop','freqAllelesFounder',
				'freqGenotypeFounder','SNPinfo','nQtl','nMarkers','nSnp','chrom',
				'nEpisPairLoci','episPairLoci','keepVariables')
rm(list=setdiff(ls(), keepVariables))

##################################--------****--------##################################

