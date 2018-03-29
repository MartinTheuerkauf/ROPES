
# ROPES requires some functions from the "discover" package (avaialble at http://disqover.botanik.uni-greifswald.de/) 
# and from the "DEoptim" package (available from the CRAN)
# importing the .ods example files requires the "readODS" package

install.packages("readODS")
library(readODS)
library(disqover)
library(DEoptim)

setwd(".../ExampleData/")

# load data (needs readODS package)
pc <- read_ods("TSK 250 Counts.ods")
par <- read_ods("TSK 250 PAR.ods")
params <- read_ods("params.ods")
ppelim <- read_ods("ppe_limits.ods")
dBasin <- 600

# ROPES with standard parameters
ROPES_TSK250 <- ROPESinR(pc=pc, par=par, params=params, ppelim=ppelim, dBasin = dBasin)

# ROPES applied with noise in the pollen counts
ROPES_TSK250_perc_noise <- ROPESinR(pc=pc, par=par, params=params, ppelim=ppelim, dBasin = dBasin, pcError = TRUE)

# ROPES applied with noise in the pollen counts and PARs
ROPES_TSK250_perc_and_par_noise <- ROPESinR(pc=pc, par=par, params=params, ppelim=ppelim, dBasin = dBasin, parError = TRUE, parErrorValue = 10, pcError = TRUE)

