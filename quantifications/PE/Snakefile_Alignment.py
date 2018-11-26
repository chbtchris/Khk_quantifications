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


for s in TISSUE_SAMPLE:
    qorts_dir = config["directories"]["qorts"]["work"] + s + "_qorts"
    if not os.path.isdir(qorts_dir):
        shell('mkdir {d}'.format(d = qorts_dir))


###########
## List final outputs
ALIGNMENTS = expand(config["directories"]["alignment"]["work"] + "{sample}Aligned.out.bam", sample = TISSUE_SAMPLE)
JUNCTIONS  = expand(config["directories"]["alignment"]["work"] + "{sample}SJ.out.tab", sample = TISSUE_SAMPLE)

SORTED_BAM = expand(config["directories"]["alignment"]["work"] + "{sample}.sortedByCoord.bam", sample = TISSUE_SAMPLE)
INDEXES = expand(config["directories"]["alignment"]["work"] + "{sample}.sortedByCoord.bam.bai", sample = TISSUE_SAMPLE)

DESEQ = expand(config["directories"]["qorts"]["work"] + "{sample}_qorts/QC.geneCounts.formatted.for.DESeq.txt.gz", sample = TISSUE_SAMPLE)
DEXSEQ = expand(config["directories"]["qorts"]["work"] + "{sample}_qorts/QC.exonCounts.formatted.for.DEXSeq.txt.gz", sample = TISSUE_SAMPLE)
JUNC = expand(config["directories"]["qorts"]["work"] + "{sample}_qorts/QC.spliceJunctionAndExonCounts.forJunctionSeq.txt.gz", sample = TISSUE_SAMPLE)
NOVEL = expand(config["directories"]["qorts"]["work"] + "{sample}_qorts/QC.spliceJunctionCounts.novelSplices.txt.gz", sample = TISSUE_SAMPLE)
KNOWN = expand(config["directories"]["qorts"]["work"] + "{sample}_qorts/QC.spliceJunctionCounts.knownSplices.txt.gz", sample = TISSUE_SAMPLE)


#### Rules
rule all:
	input:
		ALIGNMENTS, 
		JUNCTIONS, 
		SORTED_BAM, 
		INDEXES,
		DESEQ, 
		DEXSEQ, 
		JUNC, 
		NOVEL, 
		KNOWN

## Star alignment
rule star_align:
	input:
		expand(config["directories"]["fastq"]["work"] + "{{tissue}}/{{sample}}_{read}.fastq.gz", read = READS)
	output:
		align = config["directories"]["alignment"]["work"] + "{tissue}_{sample}Aligned.out.bam",
		junc = config["directories"]["alignment"]["work"] + "{tissue}_{sample}SJ.out.tab"
	params:
		genome = config["reference"]["star"],
		barcode_length = config["barcode"]["length"],
		outdir = config["directories"]["alignment"]["work"]
	run:
 		shell(' '.join(['module load star \n',
 		'echo {input} \n',
 		'STAR --runMode alignReads --runThreadN 24 --genomeDir {params.genome}',
 		'--outSAMunmapped Within --outSAMtype BAM Unsorted --readFilesCommand gunzip -c',
 		'--clip5pNbases {params.barcode_length} --genomeLoad LoadAndKeep',
 		'--readFilesIn  {input} --outFileNamePrefix  {params.outdir}/{wildcards.tissue}_{wildcards.sample} \n']))



## Bam sorting by position
rule sort:
	input:
		config["directories"]["alignment"]["work"] + "{tissue}_{sample}Aligned.out.bam"
	output:
		config["directories"]["alignment"]["work"] + "{tissue}_{sample}.sortedByCoord.bam"
	params:
		tmp = config["directories"]["bam_sort"]["scratch"]
	run:
		shell(' '.join(['module load samtools\n',
		'mkdir -p {params.tmp}/{wildcards.sample}/ \n',
		'echo {input} \n',
		'samtools sort -m 8G -T {params.tmp}/{wildcards.sample}/ -o {output} {input} \n']))


## Bam indexing
rule bam_index:
	input:
		config["directories"]["alignment"]["work"] + "{tissue}_{sample}.sortedByCoord.bam"
	output:
		config["directories"]["alignment"]["work"] + "{tissue}_{sample}.sortedByCoord.bam.bai"
	run:
		shell(' '.join(['module load samtools\n',
		'samtools index {input} \n']))


## Generation of count tables and of extensive QC metrics
rule QoRT:
	input:
		config["directories"]["alignment"]["work"] + "{tissue}_{sample}.sortedByCoord.bam"
	output:
		deseq  = config["directories"]["qorts"]["work"] + "{tissue}_{sample}_qorts/QC.geneCounts.formatted.for.DESeq.txt.gz",
		dexseq = config["directories"]["qorts"]["work"] + "{tissue}_{sample}_qorts/QC.exonCounts.formatted.for.DEXSeq.txt.gz",
		junc = config["directories"]["qorts"]["work"] + "{tissue}_{sample}_qorts/QC.spliceJunctionAndExonCounts.forJunctionSeq.txt.gz",
		novel = config["directories"]["qorts"]["work"] + "{tissue}_{sample}_qorts/QC.spliceJunctionCounts.novelSplices.txt.gz",
		known = config["directories"]["qorts"]["work"] + "{tissue}_{sample}_qorts/QC.spliceJunctionCounts.knownSplices.txt.gz"
	params:
		qorts = config["software"]["qorts"],
		gtf = config["reference"]["gtf"],
		outdir = config["directories"]["qorts"]["work"]
	run:
		shell(' '.join(['module load java \n',
		'java -Xmx16G -jar {params.qorts} QC',
		'{input} {params.gtf} {params.outdir}{wildcards.tissue}_{wildcards.sample}_qorts']))






