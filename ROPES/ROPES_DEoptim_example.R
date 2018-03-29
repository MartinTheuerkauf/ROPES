

install.packages("readODS")
library(readODS)

setwd("D:/Arbeit/DISQOVER/GIT/ROPES/ROPES/")

# load data (needs readODS package)
pc <- read_ods("TSK 250 Counts.ods")
par <- read_ods("TSK 250 PAR.ods")
params <- read_ods("params.ods")
ppelim <- read_ods("ppe_limits.ods")
dBasin <- 600

ROPES_TSK250_PROerr <- ROPESinR(pc=pc, par=par, params=params, ppelim=ppelim, dBasin = dBasin)
params
vg <- params[,2] # read fall speed (vg) and ppes
as.double(vg)
vg
