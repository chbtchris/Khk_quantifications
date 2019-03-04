import glob
import os
#import Config

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

for i in TISSUE_SAMPLE:
	sample_dir = config["directories"]["stringtie"]["naive"] + i
	if not os.path.isdir(sample_dir):
		shell('mkdir {d}'.format(d = sample_dir))


###########
## List final outputs
NAIVE = expand(config["directories"]["stringtie"]["naive"] + "{sample}/{sample}_stringTie.gtf", sample = TISSUE_SAMPLE)


#### Rules
rule all:
	input:
		NAIVE


## StringTie
rule stringtie:
	input:
		config["directories"]["alignment"]["scratch"] + "{tissue}_{sample}.sortedByCoord.bam"
	output:
		config["directories"]["stringtie"]["naive"] + "{tissue}_{sample}/{tissue}_{sample}_stringTie.gtf"
	params:
		annotation = config["reference"]["gtf"]
	run:
 		shell(' '.join(['module load gcc/4.8.2 gdc stringtie/1.3.3b \n',
 		'echo {input} \n',
 		'stringtie {input} -G {params.annotation} -e -B -p 8 -o {output}']))


