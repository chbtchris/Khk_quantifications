import glob
import os
#import Config

#### Configuration
localrules: all
configfile: "config_DRIMseq.json"

###########
## Import sample list
CONDITIONS = list(config["conditions"].keys())

###########
## List final outputs
dir_out = config["directories"]["out"]
coefs = expand(dir_out + "{c}.drimseq_coef.salmon.no_testis.RData", c = CONDITIONS)

rule all:
	input:
		coefs#,
#		contrasts

## Run the test on DEU
rule quantify_isoforms:
	input:
		dir = lambda wildcards: config["conditions"][wildcards.condition],
		annot = lambda wildcards: config["annotations"][wildcards.condition]
	output:
		coef = dir_out + "{condition}.drimseq_coef.salmon.no_testis.RData"
	params:
		Rscript = config["software"]["drimseq"],
		out = dir_out
	run:
		shell('module load new gcc/4.8.2 r/3.5.1 \nRscript {params.Rscript} {input.dir} {params.out} {input.annot} {wildcards.condition}')

