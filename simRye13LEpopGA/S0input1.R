########################### Inputs ################################
# Note that not all inputs are used in a scenario. However, they ##
# are all needed due to laziness & minimum changes for scenarios ##
###################################################################

# Traits: Minimum of 2 traits are used. This paper uses only trait 2.
nTrait=nObs=nEbv=nRes=2
# Additive variance input. Variance in iPG; or Variance in pop 1 in pPG
aRes_covar <- matrix(c(400.0, 0.00,
                       0.0, 400.00),nrow = nTrait,byrow = T)
# Additive variance in pop 2 in pPG
aNRMS_covar <- matrix(c(320.0, 0.00,
                        0.0, 320.0),nrow = nTrait,byrow = T)
# AxA variance
aaNRMSRes_covar <- matrix(c(100.0, 0.00,
                            0.0, 100.00),nrow = nTrait,byrow = T)
# Dominance variance
dNRMSRes_covar <- matrix(c(100.0, 0.00,
                           0.0, 100.00),nrow = nTrait,byrow = T)

# criteria weights. Which values need to be met? & Other inputs
critWeight <- c(0.0,1.0,  0.0,1.0,  0.0,1.0)				# nTrait * 3
critWeightcorAbio <- matrix(rep(0.01,nTrait*nTrait), nrow=nTrait, byrow = T)	# nTrait*nTrait

# input parameters on genome
nChr <- nMarkerEachChr <- nQtlEachChr <- NULL     # give it as input or use default
nMarkerEachChr <- c(rep(2,7))                     # no markers used. Gmatrix is based on QTL
#nChr <- 7; nMarkerEachChr <- c(rep(2,nChr)); nQtlEachChr <- c(rep(50,nChr))   # nQtlEachChr <- NULL;  # 
nFounderPop=2
nFounderAnimal=100
nselfing=4
DMUavailable <- FALSE

# read Genome input 
geneticArchitectureFile="C:/ThinhInDenmark/Setup_ADAM_fortran 95/NewAdam/simulationDevelopment/epistasis/genome/geneticArchitecture.dat"
if (!file.exists(geneticArchitectureFile)) geneticArchitectureFile="/usr/home/qgg/chuthinh/Hybrimax/simRyeGenome1/genome2/geneticArchitecture.dat"
baseHaplotypesFile="C:/ThinhInDenmark/Setup_ADAM_fortran 95/NewAdam/simulationDevelopment/epistasis/genome/baseHaplotypes.dat"
if (!file.exists(baseHaplotypesFile)) baseHaplotypesFile="/usr/home/qgg/chuthinh/Hybrimax/simRyeGenome1/genome2/baseHaplotypes.dat"
nBaseHap=600

#other program parameter inputs
programdir="C:/ThinhPLANTS/Hybrimax/Paper1/code2GitHub/simRye13LEpopGA/"
if (!dir.exists(programdir)) programdir='/usr/home/qgg/chuthinh/Hybrimax/simulation4/simRye13LEpopGA/'
nQtlEachQtlInteractwith=1
methodEpistasis='alphaSim' 		# or: 'multiAlleles'
nfixedFactors=2
nranFactors=1
inputParameters <- 'biological'   #  'statisticalOne'
domdegreevar <- domdegreemean <- NULL
varResiduals <- diag(rep(100,nRes),nrow = nRes)
