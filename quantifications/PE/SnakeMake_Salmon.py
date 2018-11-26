import glob
import os

#### Configuration
localrules: all
configfile: "config.json"

###########
## Import sample list

TISSUES = list(config["samples"].keys())
print(TISSUES)

READS = ["1", "2"]

TISSUE_SAMPLE = []

for t in TISSUES:
	samples = list(config["samples"][t].keys())
	TISSUE_SAMPLE = TISSUE_SAMPLE + ["{t}_{s}".format(t = t, s = s) for s in samples]
print(TISSUE_SAMPLE)


###########
## List final outputs
SALMON = expand(config["directories"]["salmon"]["work"] + "{sample}/quant.sf", sample = TISSUE_SAMPLE)

#### Rules

rule all:
	input:
		SALMON


rule salmon_quantify_wIR:
	input:
		expand(config["directories"]["fastq"]["work"] + "{{tissue}}/{{sample}}_{read}.fastq.gz", read = READS)
	output:
		config["directories"]["salmon"]["work"] + "{tissue}_{sample}/quant.sf"
	params:
		index = config["reference"]["salmon_quasi"],
		mean = config["salmon"]["mean"],
		sd = config["salmon"]["sd"]
	run:
 		shell(' '.join(['salmon',
 		'quant -i',
 		'{params.index}',
 		'-l A -r {input}',
 		'-p 24 --dumpEq --seqBias --gcBias --numBootstraps 100 -o {output}']))

