# Correction of gene model annotations improves isoform abundance estimates: the example of ketohexokinase (Khk)

[![DOI](https://zenodo.org/badge/159229285.svg)](https://zenodo.org/badge/latestdoi/159229285)

This Git repository contains the code used to generate the analysis and figures presented 
in the manuscript. Please make sure to adjust the path to the scripts and annotations before use.

The repository is organised as follows:

..* The **quantification** folder contains all scripts necessary to generate Salmon estimates and alignments
to the reference genomes (together with the associated count tables)

..* The **gene_expression** folder contains the scripts used to generate the analysis relating to the tissue
specific expression of Khk.

..* The **junctionSeq** folder contains all scripts used to estimate differential exon and junction usage

..* The **drimseq** folder should enable users to run the DRIMSeq analysis on an high performance cluster

⋅⋅* Finally, the annotation files design.txt, tissue_annotations.csv and sollner_sample_annotation.csv
contain all the sample annotations for all analysis presented in the study
