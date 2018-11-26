#---------------------------------------------
#   JunctionSeq run on computing cluster
#

## Load relevant libraries
library(JunctionSeq)
library(BiocParallel)

print("Estimate effect sizes")
jscs <- get(load("jscs.5Groups.tests.RData"))
jscs <- estimateEffectSizes(jscs)
gc()

print("Writing object")
save(jscs,
     file = "jscs.5Groups_FC.RData")




