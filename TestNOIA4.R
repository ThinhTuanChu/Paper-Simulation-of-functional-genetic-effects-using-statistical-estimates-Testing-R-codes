############### Required functions ###################
##### Set up functional matrice Wfkl in NOIA formula of pPG scenarios.
# Required input: 
#    Vectors of tx1 & tx2 are possible genotype codes at locus x (k or l) 
#    for parental population 1 and 2. Genotype codes can be mean of genotypes.
#    Elements of the vectors can be real or integers within 0 to 2.
setupWfklplot <- function (tx1=c(0,1,2),tx2=c(0,1,2)) {
  ngeno1=length(tx1);ngeno2=length(tx2)
  if (ngeno1<2 | any(tx1<0) | any(tx1>2) | ngeno2<2 | any(tx2<0) | any(tx2>2)) {
  stop ('Wrong input to setupWfklplot function.')}
  # Calculate functional additive covariates
  ta1 = matrix(kronecker(tx1-1,rep(1,ngeno2)),nrow=ngeno1*ngeno2,ncol=1)/2
  ta2 = matrix(kronecker(rep(1,ngeno1),tx2-1),nrow=ngeno1*ngeno2,ncol=1)/2
  temp1 <- kronecker(tx1,rep(1,ngeno2)); temp2 <- kronecker(rep(1,ngeno1),tx2)
  # Use Eq. 9 to calculate functional dominance covariates
  td = matrix((1-temp1/2)*(temp2/2) + (temp1/2)*(1-temp2/2),nrow=ngeno1*ngeno2,ncol=1)
  Jvec = matrix(1,nrow=ngeno1*ngeno2,ncol=1)
  Wfkl <- matrix(nrow=ngeno1*ngeno2*ngeno1*ngeno2,ncol=8)
  Wfkl[,1] <- 1.0
  Wfkl[,2] <- kronecker(Jvec,ta1);  Wfkl[,3] <- kronecker(ta1,Jvec)
  Wfkl[,4] <- kronecker(Jvec,ta2);  Wfkl[,5] <- kronecker(ta2,Jvec)
  Wfkl[,6] <- kronecker(Jvec,td);   Wfkl[,7] <- kronecker(td,Jvec)
  temp1= (Wfkl[,2:3] + Wfkl[,4:5]); Wfkl[,8] <- temp1[,1]*temp1[,2]
  return(Wfkl) }
##### Set up statistical matrice Wskl in NOIA formula of pPG scenarios.
# Required input: 
#    Vectors of tx1 & tx2 are possible genotype codes at locus x (k or l) 
#    from parental population 1 and 2, respectively. Genotype codes can 
#    be mean of genotypes. Elements of the vectors can be real or integers 
#    within 0 to 2.
#    p1k and p2k are allele frequencies at locus k from parental population 1 and 2.
#    p1l and p2l are allele frequencies at locus l.
setupWsklPlot <- function(tx1=c(0,1,2),p1k,p2k,tx2=c(0,1,2),p1l,p2l) {
  ngeno1=length(tx1); ngeno2=length(tx2)
  if (ngeno1<2 | any(tx1<0) | any(tx1>2) | ngeno2<2 | any(tx2<0) | any(tx2>2)) {
  stop ('Wrong input to setupWsklPlot function.')}
  # Use Eq. 2 to calculate statistical additive covariates
  ta1k = matrix(kronecker((tx1-2*(1-p1k)),rep(1,ngeno2)),nrow=ngeno1*ngeno2,ncol=1)/2
  ta1l = matrix(kronecker((tx1-2*(1-p1l)),rep(1,ngeno2)),nrow=ngeno1*ngeno2,ncol=1)/2
  ta2k = matrix(kronecker(rep(1,ngeno1),(tx2-2*(1-p2k))),nrow=ngeno1*ngeno2,ncol=1)/2
  ta2l = matrix(kronecker(rep(1,ngeno1),(tx2-2*(1-p2l))),nrow=ngeno1*ngeno2,ncol=1)/2
  # Use Eq. 10 to calculate statistical dominance covariates
  temp1 <- kronecker(tx1,rep(1,ngeno2)); pa11=1-temp1/2; pa12=temp1/2
  temp2 <- kronecker(rep(1,ngeno1),tx2); pa21=1-temp2/2; pa22=temp2/2 
  temp0=pa11*pa22*2*(1-p1k)*p2k + pa12*pa21*2*p1k*(1-p2k) +
    pa11*pa21*(-2)*(1-p1k)*(1-p2k) + pa12*pa22*(-2)*p1k*p2k
  tdk = matrix(temp0,nrow=ngeno1*ngeno2,ncol=1)
  temp0=pa11*pa22*2*(1-p1l)*p2l + pa12*pa21*2*p1l*(1-p2l) +
    pa11*pa21*(-2)*(1-p1l)*(1-p2l) + pa12*pa22*(-2)*p1l*p2l
  tdl = matrix(temp0,nrow=ngeno1*ngeno2,ncol=1)
  Jvec = matrix(1,nrow=ngeno1*ngeno2,ncol=1)
  # Setting up Wskl matrix
  Wskl <- matrix(nrow=ngeno1*ngeno2*ngeno1*ngeno2,ncol=8)
  Wskl[,1] <- 1.0
  Wskl[,2] <- kronecker(Jvec,ta1k); Wskl[,3] <- kronecker(ta1l,Jvec)
  Wskl[,4] <- kronecker(Jvec,ta2k); Wskl[,5] <- kronecker(ta2l,Jvec)
  Wskl[,6] <- kronecker(Jvec,tdk);  Wskl[,7] <- kronecker(tdl,Jvec)
  temp1 = (Wskl[,2:3] + Wskl[,4:5]); Wskl[,8] <- temp1[,1]*temp1[,2]
  return(Wskl) }

#######################################################
############## Testing NOIA formulas ##################

#### Random inputs ####
Ef=matrix(1,nrow = 8,ncol=1)
Wxk=c(1,rep(runif(1),2),runif(1)); Wxl=c(1,rep(runif(1),2),runif(1))
temp1=kronecker(Wxl[c(1,2)],Wxk[c(1,2)]); Ef[2:3,]=temp1[2:3]
temp1=kronecker(Wxl[c(1,3)],Wxk[c(1,3)]); Ef[4:5,]=temp1[2:3]
temp1=kronecker(Wxl[c(1,4)],Wxk[c(1,4)]); Ef[6:7,]=temp1[2:3]
Ef[8,]= runif(1)
p1k=runif(1, min = 0.01, max = 0.99); p2k=runif(1, min = 0.01, max = 0.99)
p1l=runif(1, min = 0.01, max = 0.99); p2l=runif(1, min = 0.01, max = 0.99)

#### Default input ####
# Adjustment to ensure invertable matrix. It generally has no effects on results.
invertcrit = 1e-10 ;   diagadd = 1e-8

#### Setting up matrices
Wfkl=setupWfklplot(tx1=c(0,1,2),tx2=c(0,0.5,1,1.5,2))
Wskl=setupWsklPlot(tx1=c(0,1,2),p1k,p2k,tx2=c(0,0.5,1,1.5,2),p1l,p2l)

#### Calculate Es from Ef
temp1 <- t(Wskl)%*%Wskl
if (det(temp1)<invertcrit) diag(temp1)=diag(temp1)+diagadd
Es=solve(temp1)%*%t(Wskl)%*%Wfkl%*%Ef
print(data.frame(Ef,Es))
#### Calculate Ef from Es
Es=matrix(c(1,c(runif(7))),nrow = 8,ncol=1)
temp1 <- t(Wfkl)%*%Wfkl
if (det(temp1)<invertcrit) diag(temp1)=diag(temp1)+diagadd
Ef=solve(temp1)%*%t(Wfkl)%*%Wskl%*%Es
print(data.frame(Ef,Es))
