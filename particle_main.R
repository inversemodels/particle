rm(list = ls())
if(!is.null(dev.list())) dev.off()
cat("\014")

# Libraries---------------------------------

# Load necessary packages - smcDefaults() is a list of relevant packages
library(devtools)
library(cnpParticles)
library(Rage)
library(popbio)
library(pbapply)
library(popdemo)
library(seewave)
library(npreg)
library(gtools)
library(parallel)
library(usethis)
library(future)
library(future.apply)
library(purrr)
library(furrr)
library(stringr)
library(scales)
library(tidyr)


#** ALL FUNCTIONS----
# Functions for towardParal_15

newIPM <- function(survival, fertility, offspringSize, establishment, reproduction, growth){
  IPMplaceholder <- list(survival, fertility, offspringSize, establishment, reproduction, growth)
  rm(list = c("survival", "fertility", "offspringSize", "establishment", "reproduction", "growth"), envir = .GlobalEnv)
  return(IPMplaceholder)
}

z_Length <- function(IPMListElement){
  return(length(c(which(names(formals(IPMListElement))=="z"),
                  which(names(formals(IPMListElement))=="z_prime"))))
}


totalVarLength <- function(IPMListElement){
  return(length(formals(IPMListElement)))
}

variableLength <- function(IPMListElement){
  return(totalVarLength(IPMListElement)-z_Length(IPMListElement))
}

ParamFraming <- function(IPM_function_list){
  survLen <- variableLength(IPM_function_list$survival)
  growthLen <- variableLength(IPM_function_list$growth)
  reproLen <- variableLength(IPM_function_list$reproduction)
  return(c(survLen, growthLen, reproLen))
}

genFraming<- function(IPM_function_list){
  param <- ParamFraming(IPM_function_list)
  thisParaSeries <- data.frame(matrix(NA, nrow = sum(param), ncol = 1))
  colnames(thisParaSeries) <- "paramValue"
  
  rownames(thisParaSeries)[1:param[[1]]] <- names(formals(IPM_function_list$survival))[
    (length(names(formals(IPM_function_list$survival)))-variableLength(IPM_function_list$survival)+1):
      length(names(formals(IPM_function_list$survival)))
  ]
  
  
  rownames(thisParaSeries)[(param[[1]]+1):(param[[1]]+param[[2]])] <- names(formals(IPM_function_list$growth))[
    (length(names(formals(IPM_function_list$growth)))-variableLength(IPM_function_list$growth)+1):
      length(names(formals(IPM_function_list$growth)))
  ]
  
  rownames(thisParaSeries)[(param[[1]]+param[[2]]+1):(param[[1]]+param[[2]]+param[[3]])] <- names(formals(IPM_function_list$reproduction))[
    (length(names(formals(IPM_function_list$reproduction)))-variableLength(IPM_function_list$reproduction)+1):
      length(names(formals(IPM_function_list$reproduction)))
  ]
  
  return(thisParaSeries)
}

insertParam <- function(parameter, value, paramFrame){
  newFrame <- paramFrame
  newFrame[which(parameter == rownames(newFrame)),] <- value
  return(newFrame)
}

unknownToList <- function(newUnknParams){
  unknParams <- list()
  for(i in 1:length(newUnknParams)){
    unknParams[[i]] <- solveParameter(newUnknParams[[i]])
  }
  return(unknParams)
}


funcFraming <- function(IPM_function_list){
  param <- ParamFraming(IPM_function_list)
  thisParaSeries <- data.frame(matrix(NA, nrow = sum(param), ncol = 1))
  colnames(thisParaSeries) <- "paramValue"
  
  return(list(1:param[[1]],
              (param[[1]]+1):(param[[1]]+param[[2]]),
              (param[[1]]+param[[2]]+1):(param[[1]]+param[[2]]+param[[3]])))
}


matrix.power <- function(A, n) {   # only works for diagonalizable matrices
  e <- eigen(A)
  M <- e$vectors   # matrix for changing basis
  d <- e$values    # eigen values
  return(M %*% diag(d^n) %*% solve(M))
}

simpleProjection <- function(mat,
                             timepoint,
                             freq){
  holdProject <- matrix.power(mat, timepoint)
  preNormalised <- holdProject %*% freq
  projection <- (preNormalised/sum(preNormalised))
  return(Re(projection))
}



IPMify <- function(knownParameters, ...){
  holdNewPars <- list(...)
  knownParameters[which(is.na(knownParameters))] <- unlist(holdNewPars)
  span <- z[2] - z[1]
  
  s_kernel <- IPMlist$survival(z,
                               knownParameters$survInt,
                               knownParameters$survSlope)
  
  g_kernel <- span * outer(z,
                           z_prime,
                           IPMlist$growth,
                           growthInt = knownParameters$growthInt, 
                           growthSlope = knownParameters$growthSlope,
                           growthSD = knownParameters$growthSD)
  
  f_kernel <- span * outer(z,
                           z_prime,
                           IPMlist$reproduction,
                           establishment = knownParameters$establishment, 
                           recruitSize = knownParameters$recruitSize,
                           fertInt = knownParameters$fertInt,
                           fertSlope = knownParameters$fertSlope)
  
  p_kernel <- g_kernel
  
  for (i in 1:length(z)) {
    p_kernel[, i] = g_kernel[, i] * s_kernel[i]
  }
  
  a_kernel <- p_kernel+f_kernel
  
  return(a_kernel)
}

IPMdistance_3d <- function(IPM.pre,
                           observed.stage.dist,
                           evaluate.params,
                           initial.stage.dist,
                           time.projection,
                           i){
  
  fitted_IPM <- IPMify(IPM.pre, c(evaluate.params[i, 1], 
                                  evaluate.params[i, 2]))
  
  particle_StageDist <- NA
  particle_StageDist <- try(simpleProjection(fitted_IPM,
                                             time.projection,
                                             initial.stage.dist))
  distance <- ShannonDistance(observed.stage.dist,
                              particle_StageDist)
  if (is.na(distance)) {
    distance <- 0
  }
  return(distance)
}



singleChain <- function(k){
  iRegister <- list()
  iRegister[[1]]  <- cnpParticles::initiateSMC(model.parameters = PriorsList[[k]],
                                               particle.number = particleNumber)
  for(i in 1:depth){
    pre_iRegister <- list()
    for(r in time){
      startingTime <- Sys.time()
      pre_iRegister[[i]] <- list()
      try(pre_iRegister[[i]][[r]] <- pbmapply(IPMdistance_3d,
                                              1:particleNumber,
                                              MoreArgs = list(IPM.pre = IPM_pre,
                                                              observed.stage.dist = ref_observed[[r]],
                                                              initial.stage.dist = initial_distribution,
                                                              time.projection = r,
                                                              evaluate.params = iRegister[[i]])
      ))
    }
    
    try(iRegister[[i]][,length(true_IPMparams[,1])+1] <- Reduce(`+`,
                                                                Filter(Negate(is.null),
                                                                       pre_iRegister[[i]])))
    
    try(prob_Surface <- multiSurface(perf.dataframe = iRegister[[i]],
                                     div = discretisationPar))
    
    try(new_Domain <- subDomain(surface.object = prob_Surface,
                                cut.threshold = tabuThreshold))
    
    try(new_Sample <- sirSampling(new_Domain,
                                  temperature[[i]],
                                  particleNumber))
    
    try(iRegister[[i+1]] <- new_Sample)
    try(print(paste("local", i, "|", "iteration", k)))
    endTime <- Sys.time()
    print(endTime-startingTime)
  }
  return(iRegister)
}


#************************************* ------------
# Model Definition (IPM STARTS HERE) ========================================

#--------------- Start of IPM

#* Survival ~ z --------------------------------------------------------------
survival <- function(z,
                     survInt,
                     survSlope){
  rho <- survInt + survSlope * z
  return(rho)
}

#* growth ~ z, z' -------------------------------------------------------------
growth <- function(z_prime,
                   z,
                   growthInt,
                   growthSlope,
                   growthSD){
  growthProb <- dnorm(z_prime, mean = growthInt + growthSlope * z, sd = growthSD)
  return(growthProb)
}

#* reproduction ~ z, z' --------------------------------------------------------
# fertility ~ z
fertility <- function(z,
                      fertInt,
                      fertSlope){ # Could be linear or logistic
  initialProfile <- fertInt + fertSlope * z
  truncationLine <- dunif(z, min = 3, max = 15) * (1/dunif(z, min = 3, max = 15))
  truncateProfile <- initialProfile*truncationLine
  return(truncateProfile)
}

# offspring size ~ z'
offspringSize <- function(z_prime,
                          recruitMean,
                          recruitSD){ # Could be a Gamma
  estabSize <- dnorm(z_prime, mean = recruitMean, sd = recruitSD)
  return(estabSize)
}

establishment <- NULL

# total reproductive contributions ~ z, z'
reproduction <- function(z_prime,
                         z,
                         fertInt,
                         fertSlope,
                         establishment,
                         recruitSize){
  return(fertility(z, fertInt, fertSlope)*
           establishment*
           dnorm(z_prime, conditionalSize(recruitSize), 0.2))
}

#--------------- End of IPM


# Vital Rate Functions - To - IPM List -----------------------------------------------------------------

vitalRates <- list(survival, fertility, offspringSize,
                   establishment, reproduction, growth)
names(vitalRates) <- c("survival", "fertility", "offspringSize",
                       "establishment", "reproduction", "growth")

IPMlist <- newIPM(survival, fertility, offspringSize,
                  establishment, reproduction, growth)
names(IPMlist) <- c("survival", "fertility", "offspringSize",
                    "establishment", "reproduction", "growth")

param <- ParamFraming(IPMlist)
thisParaSeries <- data.frame(matrix(NA, nrow = sum(param), ncol = 1))
colnames(thisParaSeries) <- "paramValue"

new_IPMframe <- genFraming(IPMlist)



#************************************* ------------
#*** Specify All IPM parameters -------------------------------------

IPM_LH1 <- insertParam("survInt", 0.8301, new_IPMframe) # -3.9650
IPM_LH1 <- insertParam("survSlope", -0.0715, IPM_LH1) # 0.5000
IPM_LH1 <- insertParam("growthInt", 3.641, IPM_LH1) # 2.6810
IPM_LH1 <- insertParam("growthSlope", 0.2563, IPM_LH1) # 0.5790
IPM_LH1 <- insertParam("growthSD", 0.8850, IPM_LH1) # 0.8850
IPM_LH1 <- insertParam("fertInt", -1.7328, IPM_LH1) #1.3320
IPM_LH1 <- insertParam("fertSlope", 0.5945, IPM_LH1) #0.200
IPM_LH1 <- insertParam("establishment", 1, IPM_LH1) # 0.0198
IPM_LH1 <- insertParam("recruitSize", 1.9, IPM_LH1) # 1.5240


lifeHistory <- list(IPM_LH1)
names(lifeHistory) <- c("slow")


# Creating a linked prior channel

pPriors <- IPM_LH1

pPriors$nar_Par1 <- rep(-5, length(pPriors$paramValue))
pPriors$nar_Par2 <- rep(5, length(pPriors$paramValue))

priorTabs <- list()

priorTabs[1] <- list(pPriors)

pPriors <- IPM_LH1
pPriors$broad_Par1 <- rep(-5, length(pPriors$paramValue))
pPriors$broad_Par2 <- rep(5, length(pPriors$paramValue))

priorTabs[2] <- list(pPriors)

names(priorTabs) <- c("narrow", "broad")


#


Stages <- list(c(c(rep(0, 99), 1)))
names(Stages) <- c("right")


#


metaCross <- crossing(c("slow"),
                      c("right"),
                      c("narrow", "broad"),
                      c("narrow", "broad"))



# Sequential monte carlo system parameters ---------------------------------


iKey <- IndexKey #needs changing back to index key

chains <- 1 # replicates / number of independent SMC iterations

depth <- 35 # replicates / number of independent SMC iterations

particleNumber <- 1000 # number of samples / particles

discretisationPar <- 25 # number of samples / particles

temperature <- rep(0.1, 50) # Sequence of parameters weighing the performance measures for resampling

tabuThreshold <- 0.05 # Minimum percentage of parameter domains restricted from resampling


# Data - Initial stage structure + st/age vectors ---------------------------

time <- c(1000) # time step(s) at which data is evaluated



# Manual Sim Keys --------------------------------------------------------

initCross <- crossing(rownames(IPM_LH1), rownames(IPM_LH1))
uniqCross <- unique(t(apply(initCross, 1, sort)))

nonUniq <- list()
for(i in 1:length(uniqCross[,1])){
  nonUniq[[i]] <- sum(duplicated(uniqCross[i,]))>0
}

uniqueParamCross <- data.frame(uniqCross)
uniqueParamCross$unique <- unlist(nonUniq)

uniqueParamCross <- uniqueParamCross[-which(uniqueParamCross$unique==TRUE),]
uniqueParamCross <- uniqueParamCross[,1:2]
colnames(uniqueParamCross) <- c("par1", "par2")
rownames(uniqueParamCross) <- 1:length(uniqueParamCross[,1])

triTryp <- uniqueParamCross

mainUnknownDatFrame <- crossing(triTryp, metaCross)


# Populating the unknowns 

initialScaling <- Stages[which(names(Stages)==mainUnknownDatFrame[iKey, 4][[1]])][[1]]

IPM_complete <- lifeHistory[which(names(lifeHistory)==mainUnknownDatFrame[iKey, 3][[1]])][[1]]

newUnknParams <- c(mainUnknownDatFrame[iKey, 1][[1]],
                   mainUnknownDatFrame[iKey, 2][[1]])

wPriors_1 <- priorTabs[which(names(priorTabs)==mainUnknownDatFrame[iKey, 5][[1]])][[1]]

wPriors_2 <- priorTabs[which(names(priorTabs)==mainUnknownDatFrame[iKey, 6][[1]])][[1]]


#*** Specify Unknown IPM Parameters --------------------------------------



# The below requires a pull from cnpParticles
model_Parameters <- unknownToList(newUnknParams) # Create a list of parameters w/ paired prior

model_Parameters <- model_Parameters[rank(match(newUnknParams, rownames(IPM_complete)))]

names(model_Parameters) <- newUnknParams[rank(match(newUnknParams, rownames(IPM_complete)))] # Assign list elements the names of the vectors

list2env(model_Parameters, envir = globalenv()) # Push each element of the list to an object in R

true_IPMparams <- trueParameter(param.name = names(model_Parameters), # Format unknowns to table
                                param.dataframe = IPM_complete)

IPM_paramaters <- rmParameter(param.name = names(model_Parameters), # Format complete param table
                              param.dataframe = IPM_complete)

IPM_pre <- IPMparameters(param.dataframe = IPM_paramaters) # Format parameter dataframe as list


# IPM mesh definition ---------------------------------------------------------

z <- IPMmesh(IPMboundary(minValue = 2, maxValue = 8, jointBuffer = 0.2),
             mesh.resoution = 100)

z_prime <- IPMmesh(IPMboundary(minValue = 2, maxValue = 8, jointBuffer = 0.2),
                   mesh.resoution = 100)


#**** IPM Generation ---------------------------------------------------

# Re-order the unknowns to ensure they slot into the right places in the gen function
true_IPMparams[rank(match(rownames(true_IPMparams), names(IPM_pre))),1] <- true_IPMparams[,1]
rownames(true_IPMparams)[rank(match(rownames(true_IPMparams), names(IPM_pre)))] <- rownames(true_IPMparams)

# Reference parameterisation -------------------------
# Generate the "true" (known/a priori)  Integral Population Model
IPM_observed <- IPMify(knownParameters = IPM_pre,
                       true_IPMparams[[1]])


lambda(IPM_observed)
net.reproductive.rate(IPM_observed)
generation.time(IPM_observed)


print(mainUnknownDatFrame[11,])


#-------



# Generate the "true"  st/age distribution - sub w/ projection
initial_distribution <- initialScaling

# Priors -------------------------
# Setting up a prior series

# Rank order is working correctly. When the focal parameter changes from position 1 to 2, we need to track with priors.

mean_1 <- true_IPMparams[[1]][1]
mean_2 <- true_IPMparams[[1]][2]


PriorsList <- list()

model_Parameters[[1]]$prior <- function(x){return(mean_1 + runif(x,
                                                                 min = wPriors_1[which(names(model_Parameters)[[2]]==rownames(wPriors_1)), 2],
                                                                 max = wPriors_1[which(names(model_Parameters)[[2]]==rownames(wPriors_1)), 3]))}


model_Parameters[[2]]$prior <- function(x){return(mean_2 + runif(x,
                                                                 min = wPriors_2[which(names(model_Parameters)[[2]]==rownames(wPriors_2)), 2],
                                                                 max = wPriors_2[which(names(model_Parameters)[[2]]==rownames(wPriors_2)), 3]))}


# Modify the order of parameters

for(i in 1:chains){
  PriorsList[[i]] <- model_Parameters
}


# Projection + Distance -------------------------

simpleRef <- list()

for(i in time){
  simpleRef[[i]] <- simpleProjection(IPMify(IPM_pre,
                                            true_IPMparams[[1]][1],
                                            true_IPMparams[[1]][2]), 
                                     i,
                                     initial_distribution)
}


ref_observed <- simpleRef


# Evaluate Param Performance -------------------------

iNestedRegister <- list()

#---

outobject <- lapply(1:chains,
                    singleChain)


save(outobject, file=paste("/data/zool-salgo-team/pemb5819/", "3155", "remote",
                           "/",
                           "id",
                           IndexKey,
                           ".",
                           str_pad(round(runif(1, 0, 10000), 0), 5, pad = "0"),
                           ".RData",
                           sep = ""))





















