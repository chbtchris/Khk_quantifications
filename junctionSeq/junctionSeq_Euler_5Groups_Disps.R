#---------------------------------------------
#   JunctionSeq run on computing cluster
#


## Load relevant libraries
library(JunctionSeq)
library(BiocParallel)
print("Libraries loaded")



## Directory where reports are stored
directory.qorts = "/qorts/"

#---------------------------------------------
## Create JunctionSeq annotations

colData = read.csv("tissue_annotations.csv",
                   stringsAsFactors = FALSE,
                   header = TRUE)

print("Decoder loaded")

rownames(colData) = as.character(colData$Sample)
colData$sample = factor(colData$Sample,
			levels = colData$Sample)


## We only focus on the 5 tissues being extracted: heart, spleen, kidney, lung and liver
colData = colData[colData$Group %in% c("kidney", "heart", "liver",  "lung", "spleen"),]


## Factor re-ordering
colData$condition = factor(colData$Group,
                           levels = c("kidney", "heart", "liver",  "lung", "spleen"))

print("colData loaded")

multicoreParam <- MulticoreParam(workers = 23)

print("Reading data")
count.files = paste(directory.qorts, 
                    paste0(colData$Sample, "_qorts"), 
                    "QC.spliceJunctionAndExonCounts.forJunctionSeq.txt.gz", sep = "/")
print(count.files)
print(colData[, c("sample", "condition")])

jscs = readJunctionSeqCounts(countfiles = count.files,
                             samplenames = colData$Sample,
                             design = colData[, c("sample", "condition")],
                             flat.gff.file = "/mouse/M14_GRCm38.p5/gtf/gencode.vM14.primary_assembly.annotation.flattened.QoRTs.junctions.gtf",
                             test.formula1 = formula(~ sample + countbin + condition : countbin),
                             analysis.type = c("junctionsAndExons"),
                             nCores = multicoreParam)

print("Object loaded")
gc()

print("Estimate size factors")
jscs <- estimateJunctionSeqSizeFactors(jscs)
gc()

print("Estimate Dispersions")
jscs <- estimateJunctionSeqDispersions(jscs, nCores = multicoreParam)
gc()

print("Fit dispersion functions")
jscs <- fitJunctionSeqDispersionFunction(jscs)
gc()

#print("Test diff usage")
#jscs <- testForDiffUsage(jscs, nCores = multicoreParam)
#gc()

print("Writing object")
save(jscs,
     file = "jscs.5Groups.Disps.RData")


