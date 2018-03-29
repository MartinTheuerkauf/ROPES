#' ROPES
#' 
#' The ROPES approach (Reveals withOut PpES; Theuerkauf and Couwenberg 2018; 
#' doi: 10.3389/feart.2018.00014 ) aims to translate pollen data from large 
#' sites into regional vegetation composition. ROPES builds on the REVEALS model
#' (Sugita 2007). However, other than REVEALS, the ROPES approach allows 
#' reconstruction without given pollen productivity estimates (known as PPEs or 
#' RPPES).  To that end, ROPES addinally employs pollen accumulation rates 
#' (PAR). PARs provide independend information about the changes in abundance 
#' for each species. ROPES uses this additional information, i.e. it fits a 
#' REVEALS reconstruction to the changes known from the PAR data. So, all 
#' changes shown in the PAR data should be similarly present in a REVEALS 
#' reconstruction. ROPES applies optimization to find the set of PPEs that 
#' produces the best matching reconstruction. To that end, ROPES minmizes 
#' variation in the ratio of PAR values over abundance reconstructed with 
#' REVEALS (PoR ratio).
#' 
#' ROPES requires both pollen counts and PAR value from a pollen records For 
#' succesful application, pollen records with sufficient variation in most 
#' pollen taxa along the record. Furthermore, succesful application is 
#' deteremined by the quality of the PAR value.
#' 
#' The ROPESinR function returns a list with VTR (the Value to be reached in
#' optimization), iter (the number of iterations needed), the resulting PPEs,
#' the resulting cover of each taxon, the initial PAR values and the PoR ratio.
#' 
#' 
#' @param pc data frame with raw pollen counts. Counts are arranged as pollen 
#'   taxa in columns and samples in rows. The first colum includes sample ages 
#'   or depth, i.e. only numeric values! The first row includes taxon names.
#'   
#' @param par data frame with pollen accumulation rates. Values are arranged 
#'   like the pollen counts, i.e. as pollen taxa in columns and samples in rows.
#'   pc and par must have the same dimensions. The first colum includes sample 
#'   ages or depth, i.e. only numeric values! The first row includes taxon 
#'   names.
#'   
#' @param params data frame of parameters. Includes 'fallspeed' of pollen (in m 
#'   s-1), relative pollen productivity estimates ('PPE') and their standard 
#'   error. Pollen taxa are given in rows, parameters in colums.  Column names 
#'   must be "species", "fallspeed", "PPE", and "PPE.error". See 'paramsTS' in 
#'   the 'Tiefer See' data set as example. For each pollen taxon in the pollen 
#'   data, one record of parameters with exactly the same name is required. 
#'   Params may include more taxa, also in different order, so that a standard 
#'   list may be established.
#'   
#' @param ppelim data frame defining the lower and upper limit of each PPE. As a
#'   suggestion, set the lower limit to 0.01 and the upper to 20. This range 
#'   covers common PPE values. For the reference taxon, both lower and upper 
#'   limit are set to 1.
#'   
#' @param parError if TRUE, noise is added in PAR data in each model run. 
#'   Magnutide of error is defined by parErrorValue. This option may be useful 
#'   in simulation studies. By default FALSE.
#'   
#' @param pcError if TRUE, noise is added in pollen counts through re-drawing. 
#'   Magnutide of noise is defined by pollen sum, here set to 1000 (may be 
#'   changed in the code). This option may be useful in simulation studies. By 
#'   default FALSE.
#'   
#' @param parErrorValue defines the magnitude of noise added in PAR data. To add
#'   noise, each PAR value is multiplied by a random value from a normal 
#'   distribution, devided by parErrorValue and than added to the original PAR 
#'   values. So higher parErrorValues produce lower noise.
#'   
#' @param dBasin basin diameter in meter
#'   
#' @param dwm distance weighting method. The following methods are implemented: 
#'   'lsm unstable', 'gpm unstable', 'gpm neutral'. 'lsm' refers to the the 
#'   Lagrangian stochastic model presented by Kuparinen et al. 2007 (Ecological 
#'   Modelling 208:177-188 ). 'gpm' refers to the Gaussian plume model, based on
#'   Sutton's equations, used by Prentice 1985 (Quaternary Research 23: 76-86). 
#'   Further methods can be added.
#'   
#' @param tBasin type of basin, either 'lake' oder 'peatland'
#'   
#' @param regionCutoff	diameter of the reconstruction region in meters, by 
#'   default 100000
#'   
#' @param strategy	control paramter of DEoptim
#'   
#' @param itermax control paramter of DEoptim, setting the maximum number of 
#'   iterations. By default 20000.
#'   
#' @param reltol control paramter of DEoptim, setting a stop criterion for 
#'   operations. By default 0.001.
#'   
#' @param steptol control paramter of DEoptim, setting a stop criterion for 
#'   operations. By default 500.
#'   
#' @param parallelType control paramter of DEoptim for multicore calculations. 
#'   With the default value of 1, all processor cores are used.
#'   


#############################################################################################
###  ROPESinR with DEoptim
#############################################################################################


ROPESinR <- function (pc, par, params, ppelim, parError = FALSE, parErrorValue, pcError = FALSE, dBasin, 
                      dwm="lsm unstable", tBasin="lake", regionCutoff=1e5, strategy=4, 
                      itermax = 20000, reltol=0.001, steptol=500, parallelType=1){
  
  vg <- unlist(params[,2]) # read fall speed (vg) and ppes
  
  # calculate K factors (parameters should provided by users)
  K <- do.call(rbind,lapply(vg, DispersalFactorK, tBasin=tBasin, dBasin=dBasin, dwm=dwm, regionCutoff=regionCutoff))
  
  # remove depth/age column
  # par <- par[,-1]
  # pc <- pc[,-1]
  # extract number of taxa and samples
  nTaxa <- nrow(ppelim)
  nSamples <- nrow(par)
  
  # define lower and upper PPE limit with data provided
  lo <- as.numeric(ppelim[,1])
  up <- as.numeric(ppelim[,2])
  
  # optionally add PAR error, range defined by parErrorValue
  if (parError) {
    noise <- rnorm(nSamples)/parErrorValue
    par <- par + par*noise
  }
  
  if (pcError) {
    pcp <- prop.table(pc,1)
    pc <- t(apply(pcp, 1, function(x) rmultinom(n=1, size=1000, prob = x))) # may in the future use disqover function
  }

  # optimization
  ropes <- DEoptim(dist.dplr.c, lo, up, pac=par, pc=pc, K=K,
                   DEoptim.control(strategy=strategy, itermax = itermax, reltol=reltol, steptol=steptol, VTR = 0, 
                                   parallelType=parallelType))
  
  # extract bestvalue, number of iterations and resulting PPEs
  vtr <- ropes$optim$bestval
  it <- ropes$optim$iter
  rppes <- ropes$optim$bestmem
  
  # apply REVEALS ro calculate resulting cover values (REVEALS step II) and PoR ratio
  f <- as.vector(K*rppes)
  RI <- t(t(pc) / f) # Reveals step I
  RII <- 100*prop.table(RI, 1) # Reveals step II
  
  return(invisible(list(VTR = vtr, iter = it, PPEs = rppes, cover = RII, counts = pc, PARs = par, PoRratio = par/RII)))
}

#############################################################################################
###  target functions to calculate first REVEALS cover and then PoR ratio
#############################################################################################

library(compiler)

# target function, calculates variation in the PoR ratio
# to that end, an overall PoR ratio is calculated from the PoR ratio of each pollen taxa using
# a tree ring method from the dplr package

dist.dplr <- function(o, pac, pc, K){
  library(dplR)
  f <- as.vector(K*o)
  RI <- t(t(pc) / f) # Reveals step I
  RII <- prop.table(RI, 1) # Reveals step II 
  por <-  pac/RII
  por[is.infinite(por)] <- NA
  porm <-apply(por, 2, function(x) x/mean(x, na.rm = TRUE))
  porm.crn <- chron(porm, prefix="HUR") #create a common PoR ratio
  dist <- var(porm.crn[,1])
  return(dist)
}
dist.dplr.c <- cmpfun(dist.dplr)

