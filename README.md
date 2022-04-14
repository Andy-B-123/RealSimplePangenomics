# RealSimplePangenomics
Analysis script for reference-based pangenome analysis

### Introduction
This script is aimed at quick assessment of a group of samples in a pan-genome context.

Each sample is expected to be aligned to the same reference and is expected to be in BAM format, sorted and indexed. 
The reference is required to have an annotation of genes, including exon and transcript descriptors.
The script runs on a single thread per and for even very large files (Mouse alignment, sorted BAM > 60GB, over 100,000 features) takes a max of 30 minutes.

### Requirements
The script requires three non-standard Python3 libraries - pysam, icecream and tqdm
Please install these locally:
https://github.com/tqdm/tqdm
`pip install tqdm --user`

https://github.com/gruns/icecream
`pip install icecream --user`

https://github.com/pysam-developers/pysam
`pip install pysam --user`

Also required are:
os, pathlib, sys, os, glob, csv, argparse, ntpath, pysam, numpy, pandas 
but these are mostly standard.

### Usage
A typical command line looks like this:
```
python3 RealSimplePangenome_v1.2.Single.py \
--input_bam <path/to/input.bam> \
--input_annotation <path/to/input.gff> \
--delim [either ' ' or '=' in case of gtf or gff] \
--feature exon \
--annot_gene_str <eg. gene_id> \
--annot_transcript_str <eg transcript_id> \
  --output_filepath ./
```
#### From top to bottom:

--input_bam      is expected to be sorted and indexed.

--input_annotation      can be either in gff or gtf format. Depending on the format you will have to specify how the key:value pairs in the annotation file are parsed. gff files have key:value pairs seperated by '=', while gtf format has these values seperated by ' ' (eg. a space).

--delim     depending on your input file edit to this to fit either gtf or gff convention

--feature     this will be the string looked for in the annotation for which to generate the proportion of coverage. Most common is 'exon' but 'CDS' can also be used

--annot_gene_str    This is the highest level feature that the exons will be collapsed to. The collapsing averages the coverage of ALL of the exons for a gene.

--annot_transcript_str     This is the second highest level that the exons will be collapsed to. The collapsing will average the coverage of the exons to the parent transcript.

--output_filepath     The output filepath to write files to.

### Outputs
The output of the script are three csv files:
  [BAM filename].raw.csv
  [BAM filename].processed.transcriptLevel.csv
  [BAM filename].processed.transcriptLevel.csv
  
The '.raw.csv' file contains the raw 1x, 2x and average coverage of every exon in the input annotation. This is a very large file.
  
The '.processed.transcriptLevel.csv' file contains the exons coverages collapsed to transcript leve. The collapsing is done by averaging the coverage of all of the exons for a parent transcript.
  
The '.processed.geneLevel.csv' file contains the exons coverages collapsed to the gene level. The collapsing is done by averaging the coverage of all of the exons for a parent gene. 

The processed files can be merged from multiple files for analyis in R or other programs. Just take the relevant column from each sample and paste them together. This can ONLY be done if you have used the same annotation file for each sample. If you have used different annotation files the order will be different and you will need to do a more complex merge!  


### Specific for users
If using this on the CSIRO HPC system (eg. Petrichor) please load the latest python module:<br />
`module load python/3.9.4`<br />


