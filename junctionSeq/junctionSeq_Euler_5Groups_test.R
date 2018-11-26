#---------------------------------------------
#   JunctionSeq run on computing cluster
#


## Load relevant libraries
library(JunctionSeq)
library(BiocParallel)


jscs = get(load("jscs.5Groups.Disps.RData"))

multicoreParam <- MulticoreParam(workers = 23)

print("Test diff usage")
jscs <- testForDiffUsage(jscs, nCores = multicoreParam)
gc()

print("Writing object")
save(jscs,
     file = "jscs.5Groups.tests.RData")

