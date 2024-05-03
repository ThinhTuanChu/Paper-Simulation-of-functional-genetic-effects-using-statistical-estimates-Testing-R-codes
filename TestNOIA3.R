######################################################
########## Code for testing NOIA formula #############
######################################################

############### required functions ###################
##### Set up functional matrix Wfkl (9x6)
# Required input: Vector of tx as possible genotype code at locus x (k or l). 
#                 In this NOIA formula, tx = c(0,1,2)
setupWfkl <- function (tx=c(0,1,2)) {
    if (length(tx)!=3 | any(tx<0)| any(tx>2)) stop ('Wrong input to setupWfkl function.')
    ta = matrix(tx-1,nrow=3,ncol=1)  # additive covariate
    tx[tx>1]=abs(tx[tx>1]-2)
    td = matrix(tx,nrow=3,ncol=1)    # dominance covariate
    Jvec = matrix(1,nrow=3,ncol=1)
    Wfkl <- matrix(nrow=3*3,ncol=6)
    Wfkl[,1]= 1.0
    Wfkl[,2]= kronecker(Jvec,ta); Wfkl[,3]= kronecker(ta,Jvec)
    Wfkl[,4]= kronecker(Jvec,td); Wfkl[,5]= kronecker(td,Jvec)
    Wfkl[,6]= Wfkl[,2]*Wfkl[,3]
  return(Wfkl) }

##### Set up statistical matrix Wskl (9x6)
# Function to set up Wx matrix as a component to calculate Wskl matrix
#          Use Eq. 2 & 3 in the main paper
#          Required input: genotype frequencies at locus x.
setupWxs <- function(px) {
  if (round(sum(px),digits=1) !=1) stop ('Wrong input to setupWxs function.')
  Wx <- matrix(nrow=3,ncol=3)
  Wx[,1]=1.0
  Wx[,2]=c((0.0-px[2]-2.0*px[3]),(1.0-px[2]-2.0*px[3]),(2.0-px[2]-2.0*px[3]))
  temp1 <- (px[3]+px[1]-((px[1]-px[3])^2))
  if (temp1!=0) {
  Wx[,3]=c((-2.0*px[2]*px[3])/temp1,( 4.0*px[1]*px[3])/temp1,(-2.0*px[1]*px[2])/temp1)
  } else {Wx[,3]=0}
  return(Wx) }
# Function to set up Wskl matrix (9x6)
#          Required input: genotype frequencies at loci k & l
setupWskl <- function(pk,pl) {
  if (length(pk)!=3 | any(pk<0)| any(pk>1)) stop ('Wrong input to setupWskl function.')
  if (length(pl)!=3 | any(pl<0)| any(pl>1)) stop ('Wrong input to setupWskl function.')
  Wskl <- matrix(nrow=3*3,ncol=6)
  Wk=setupWxs(pk); Wl=setupWxs(pl)
  temp1=kronecker(Wl[,c(1,2)],Wk[,c(1,2)]) 
  Wskl[,1:3]=temp1[,1:3]   					# a
  temp1=kronecker(Wl[,c(1,3)],Wk[,c(1,3)]) 
  Wskl[,4:5]=temp1[,2:3]   					# d
  Wskl[,6]=kronecker(Wl[,c(2)],Wk[,c(2)]) 	# aa
  return(Wskl) }

###### Function to set up diagonal matrix Dkl (9x9)
#      Required input: genotype frequencies at loci k & l
setupDkl <- function(pk,pl) {
  if (length(pk)!=3 | any(pk<0)| any(pk>1)) stop ('Wrong input to setupDkl function.')
  if (length(pl)!=3 | any(pl<0)| any(pl>1)) stop ('Wrong input to setupDkl function.')
  pk=matrix(pk,nrow = 3); pl=matrix(pl,nrow = 3)
  temp1 <- as.vector(kronecker(pl,pk))
  Dkl <- diag(temp1)
  return(Dkl)
}

#######################################################
############## Testing NOIA formulas ##################

#### Random inputs ####
temp1 = runif(2, min = 0.01, max = 0.49); pk = c(temp1,1-sum(temp1))
temp1 = runif(2, min = 0.01, max = 0.49); pl = c(temp1,1-sum(temp1))
Eskl=matrix(c(1,c(runif(5))),nrow = 6,ncol=1)

#### Default input ####
# Adjustment to ensure invertable matrix. It generally has no effects on results.
invertcrit = 1e-10
diagadd = 1e-8

# Setting up matrices
Wfkl <- setupWfkl(); Wskl <- setupWskl(pk,pl); Dkl <- setupDkl(pk,pl)

                     #### Calculate Efkl from Eskl
temp1 <- t(Wfkl)%*%Wfkl        # Test formula without Dkl
if (det(temp1)<invertcrit) diag(temp1)=diag(temp1)+diagadd
Efkl = EfNoDkl = solve(temp1)%*%t(Wfkl)%*%Wskl%*%Eskl

temp1 <- t(Wfkl)%*%Dkl%*%Wfkl  # Test formula with Dkl
if (det(temp1)<invertcrit) diag(temp1)=diag(temp1)+diagadd
Efkl = EfYesDkl = solve(temp1)%*%t(Wfkl)%*%Dkl%*%Wskl%*%Eskl

# Show & compare results between formula
print(data.frame(Eskl,EfNoDkl,EfYesDkl))

                     #### Calculate Eskl from Efkl
Efkl=matrix(c(1,c(runif(5))),nrow = 6,ncol=1)
temp1 <- t(Wskl)%*%Wskl        # Test formula without Dkl
if (det(temp1)<invertcrit) diag(temp1)=diag(temp1)+diagadd
Eskl = EsNoDkl = solve(temp1)%*%t(Wskl)%*%Wfkl%*%Efkl

temp1 <- t(Wskl)%*%Dkl%*%Wskl  # Test formula with Dkl
if (det(temp1)<invertcrit) diag(temp1)=diag(temp1)+diagadd
Eskl = EsYesDkl = solve(temp1)%*%t(Wskl)%*%Dkl%*%Wfkl%*%Efkl

# Show & compare results between formula
print(data.frame(Efkl,EsNoDkl,EsYesDkl))

