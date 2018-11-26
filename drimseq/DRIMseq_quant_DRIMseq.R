###################################################################
### Load relevant libraries
library(tximport)
library(DRIMSeq)
library(BiocParallel)


###################################################################
### Get arguments passed from command line
args <- commandArgs(TRUE)

dir.in = args[1]
dir.out = args[2]
sample_annot = args[3]
cond = args[4]


###################################################################
### Function to reformat transcript annotation
splitTxNames <- function(txNames){
  txSp = strsplit(txNames, split = "\\|")
  
  ## Extract information
  tx_id = lapply(txSp, function(x){x[1]})
  gene_id = lapply(txSp, function(x){x[2]})
  gene_name = lapply(txSp, function(x){x[6]})
  gene_type = lapply(txSp, function(x){x[8]})
  entrez_id = lapply(txSp, function(x){x[7]})
  
  ## Create DF
  df = data.frame(transcript_id = unlist(tx_id),
                  gene_id = unlist(gene_id),
                  gene_name = unlist(gene_name),
                  gene_type = unlist(gene_type),
                  entrez_id = unlist(entrez_id))
  return(df)
}


###################################################################
### Read quantification and counts from Salmon

sample.annot = read.csv(sample_annot, stringsAsFactors = F)
files <- file.path(dir.in, sample.annot$Sample, "quant.sf/quant.sf")
names(files) <- sample.annot$Sample

tx.quant <- tximport(files, type = "salmon", txOut = TRUE)
names(tx.quant)


tx.abundance = tx.quant$abundance
nonZero = !(rowSums(tx.abundance) == 0)
tx.abundance.noZero = tx.abundance[nonZero,]
tx.names = rownames(tx.abundance.noZero)


###################################################################
### Reformat transcript names

df.tx = splitTxNames(tx.names)
gc()


print("Formatting transcripts")
tx.abundance.df = data.frame(tx.abundance.noZero)
tx.abundance.df$feature_id = df.tx$transcript_id
tx.abundance.df$gene_id = df.tx$gene_id
rownames(tx.abundance.df) <- df.tx$transcript_id
head(tx.abundance.df)
gc()


###################################################################
### Format sample annotations

print("Formatting sample annotations")
drimseq_samples = data.frame(sample_id = sample.annot$Sample,
                             group = sample.annot$Group)
head(drimseq_samples)
gc()


###################################################################
### Create DRIMSeq object and filter

print("Create object")
d <- DRIMSeq::dmDSdata(counts = tx.abundance.df, samples = drimseq_samples)
d <- dmFilter(d, 
              min_samps_gene_expr = 3,
              min_samps_feature_expr = 3)


###################################################################
### DRIMSeq analysis

print("Estimate precisions")
design_full <- model.matrix(~ group, data = samples(d))
set.seed(123)
multicoreParam <- MulticoreParam(workers = 23)
d <- dmPrecision(d, design = design_full,
                 BPPARAM = multicoreParam,
                 prec_grid_range = c(-15, 15))
gc()


print("Fitting")
multicoreParam <- MulticoreParam(workers = 23)
d <- dmFit(d, design = design_full, verbose = 1,
           BPPARAM = multicoreParam)

gc()
print("Testing")
multicoreParam <- MulticoreParam(workers = 23)
d.coef <- dmTest(d, coef = "groupkidney", verbose = 1,
                 BPPARAM = multicoreParam)

## Save output
coef.out = paste0(dir.out, cond, ".drimseq_coef.salmon.RData")
print(coef.out)
print("Done")

save(d.coef,
     file = coef.out)
