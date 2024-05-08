This simulation program was used for paper: "Simulation of functional additive and 
non-additive genetic effects using statistical estimates from quantitative genetic 
models"
By Thinh Tuan Chu, QGG, Aarhus University, Denmark

For the interest of scientific transparency, I have made this R script available on
GitHub. However, this R code was meant for a preliminary exploration of the theory 
on the topic. Or, the code is only for this specific paper. It was not prepared for 
efficient computation or well-strutured program.

Please refer to paper for a description of the simulation program:
Folder 'simRye13LEpop'       Example 1, LE population, method SS_NOIA
Folder 'simRye13LEpopGA'     Example 1, LE population, method SF_GA
Folder 'simRye13singlepop'   Example 2, LD population, method SS_NOIA
Folder 'simRye13singlepopGA' Example 2, LD population, method SF_GA
Folder 'simRye13multipopGA'  Example 3, LD population, method SF_GA

To run a simulation program:
1. modify S0input1.R
     Change object geneticArchitectureFile, baseHaplotypesFile,programdir to a 
        correct location
2. (optional) modify S4run*.R for different method of calculating variance
3. execute S4run*.R

Note: 
 - Not all steps and inputs are used in a scenario. However, they are all needed due
to laziness & minimum changes needed for different scenarios.
 - To run the breeding scheme (generation 1-4), the code requires an external program 
DMU to estimate variance component & predict EBV. However, DMU is not provided together
with this code due to liciense etc. Therefore, it only generate functional QTL effects
given statistical variance. Please contact the author (chu.thinh @ qgg.au.dk remove 
space) if you'd like to run the complete breeding program.

This code is free and comes with ABSOLUTELY NO WARRANTY.
You are welcome to use and redistribute it.
