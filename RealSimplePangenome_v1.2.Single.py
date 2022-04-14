###Get the coverage of exons in a gff for a BAM alignment file
###Feels insane that this doesn't exist
###The problem with featureCounts is the counts rather than coverage - it makes it hard to set a threshold.
###Set a threshold for all exons in a given transcript
###Usage: python3 RealSimplePangenome_v1.2.Single.py --input_bam <path/to/input.bam> --input_annotation <path/to/input.gff> --delim [either ' ' or '=' in case of gtf or gff] --feature exon --annot_gene_str <eg. gene_id> --annot_transcript_str <eg transcript_id> --output_filepath ./
###
###v1.2
###Previous version didn't update value for coverage if an exon was incorrectly formatted (eg. start and end are the same), so these exons were still printed into the output but with the preceeding coverage value
###This isn't a maaaajor problem, it would mess up until a properly formatted exon was found, but still, not good.
###Updated here, these exons are now skipped properly.
###Also updated that the output tables have the sample name in the header as well as outputting two summary files to take the overhead off from R. These are the 'processed' csv files, at transcript and gene level

import os
from pathlib import Path
import sys
import os
import glob
import csv
import argparse
from icecream import ic
import ntpath
import pysam
import numpy as np
from tqdm import tqdm
import pandas as pd

pd.set_option('display.max_columns', None)
#pd.set_option('display.max_rows', None)

#First, take in a list of BAM files and also an annotation file

parser = argparse.ArgumentParser()

parser.add_argument("--input_bam", help="A text list with the full path of the BAM files to be analysed")
parser.add_argument("--input_annotation", help="An annotation file to be used to check the BAM files. Can be either GFF or GTF")
parser.add_argument("--delim", help="Delimiter for use in key:value extraction of metadata. '=' in the NCBI ref, ' ' in the Stringtie ref")
parser.add_argument("--feature", help="The target feature to extract coverage info for.",default=str("exon"))
parser.add_argument("--annot_gene_str", help="The string used in the annotation to denote the 'gene' level of an exon.",default=str("gene_id"))
parser.add_argument("--annot_transcript_str", help="The string used in the annotation to denote the 'transcript' level of an exon.",default=str("transcript_id"))
parser.add_argument("--output_filepath", help="Output filepath",default=str(""))

args = parser.parse_args()

input_bam = args.input_bam
input_annotation = args.input_annotation
output_filepath = args.output_filepath
delim = args.delim
feature = args.feature
annot_gene_str = args.annot_gene_str
annot_transcript_str = args.annot_transcript_str

anno_data = [line.strip() for line in open(input_annotation, 'r')]
#ic(anno_data)
ic(anno_data[2].split('\t'))
ic(input_bam)
ic(input_annotation)
ic(output_filepath)
ic(delim)
ic(annot_gene_str)
ic(annot_transcript_str)
ic(feature)

#os.makedirs('SingleResults', exist_ok=True)	
#os.makedirs('CombinedResults', exist_ok=True)	

#Function which takes a BAM file and returns the coverage info for the features in the annotation file

def CheckFeatureCoverage(filename):
	global base_cov0
	global base_cov1
	global base_covAvg
	samfile = pysam.AlignmentFile(filename, "rb")
	ic(samfile)

	sample_id = os.path.basename(filename)[:-4]
	ic(sample_id)

	output_data = []
	sample_cov1 = []
	sample_cov2 = []
	sample_covAvg = []

	for line in tqdm(anno_data):
		if line.startswith("#"):
			continue
		chr, start, stop, geneID, source = line.split('\t')[0], int(line.split('\t')[3]) ,int(line.split('\t')[4]), line.split('\t')[8],line.split('\t')[2]
		if source == 'transcript':
			continue
		elif source == feature:
			metadata = dict()
			for item in geneID.split(';'):
				data_split = dict()
				if len(item.split(delim)) > 1:
					data_split[item.lstrip().split(delim)[0]] = item.lstrip().split(delim)[1]
					metadata.update(data_split)
			for key,value in metadata.items():
				metadata[key] = value.replace('"',"")

			try:
				coverage = samfile.count_coverage(chr,start,stop)
			except:
				print("Error in this line:")
				print(line)
				print('SKIPPING THIS FEATURE')
				coverage = []
				continue
			gene_coverage_list = []
			for letter in coverage:
				gene_coverage_list.append(letter.tolist())
			gene_coverage_combined = zip(gene_coverage_list[0],gene_coverage_list[1],gene_coverage_list[2],gene_coverage_list[3])
			gene_coverage_combined_sum = [x + y + z + a for (x,y,z,a) in gene_coverage_combined]
			exon_cov_avg = sum(gene_coverage_combined_sum) / len(gene_coverage_combined_sum)
			#ic("-------")
			#ic(exon_cov_avg)
			exon_cov_avg = round(exon_cov_avg,3)
			#ic(exon_cov_avg)
			exon_cov_over0 = sum(x > 0 for x in gene_coverage_combined_sum) / len(gene_coverage_combined_sum)
			#ic(exon_cov_over0)
			exon_cov_over0 = round(exon_cov_over0,3)
			#ic(exon_cov_over0)
			exon_cov_over1 = sum(x > 1 for x in gene_coverage_combined_sum) / len(gene_coverage_combined_sum)
			#ic(exon_cov_over1)
			exon_cov_over1 = round(exon_cov_over1,3)
			#ic(exon_cov_over1)
			data_out = {'chr' : chr, 'start' : start, 'stop' : stop, f'{sample_id}.exon_cov_avg' : exon_cov_avg, f'{sample_id}.exon_cov_over0' : exon_cov_over0, f'{sample_id}.exon_cov_over1' : exon_cov_over1}
			data_out.update(metadata)
			output_data.append(data_out)

		else:
			continue
			#ic("Not an exon or a transcript?")

	raw = pd.DataFrame.from_dict(output_data)
	ic(raw)
	raw.to_csv(f'{output_filepath}/{sample_id}.raw.csv',index=False)

	ic(annot_gene_str)
	ic(annot_transcript_str)
	
	processed_df_transcriptLevel = raw.groupby([annot_gene_str,annot_transcript_str])[f'{sample_id}.exon_cov_avg',f'{sample_id}.exon_cov_over0',f'{sample_id}.exon_cov_over1'].mean()
	processed_df_transcriptLevel = processed_df_transcriptLevel.round(4)
	ic(processed_df_transcriptLevel)
	
	basic_data = raw.groupby([annot_gene_str,annot_transcript_str]).agg({'chr' : [pd.Series.mode], 'start' : ['min'], 'stop': ['max']})
	ic(basic_data)
	ic(basic_data.columns.get_level_values(0))
	basic_data.columns = basic_data.columns.get_level_values(0)
	ic(basic_data)

	processed_transcript_output = basic_data.merge(processed_df_transcriptLevel, on=[annot_gene_str,annot_transcript_str])
	processed_transcript_output.to_csv(f'{output_filepath}/{sample_id}.processed.transcriptLevel.csv')

	processed_df_geneLevel = processed_df_transcriptLevel.groupby([annot_gene_str])[f'{sample_id}.exon_cov_avg',f'{sample_id}.exon_cov_over0',f'{sample_id}.exon_cov_over1'].mean()
	processed_df_geneLevel = processed_df_geneLevel.round(4)
	ic(processed_df_geneLevel)

	basic_data_gene = basic_data.groupby([annot_gene_str]).agg({'chr' : [pd.Series.mode], 'start' : ['min'], 'stop': ['max']})
	ic(basic_data_gene)
	ic(basic_data_gene.columns.get_level_values(0))
	basic_data_gene.columns = basic_data_gene.columns.get_level_values(0)
	ic(basic_data_gene)

	processed_gene_output = basic_data_gene.merge(processed_df_geneLevel, on=[annot_gene_str])
	processed_gene_output.to_csv(f'{output_filepath}/{sample_id}.processed.geneLevel.csv')


samples = []
if input_bam.endswith(".bam"):
	print("Single input file detected. Moving to single analysis")
	samples.append(input_bam)
else:
	with open(input_bam) as file:		#Take in the filepaths from the BAM file

		for line in file: 
			line = line.strip()
			if not line.lower().endswith(".bam"):		#check that they at least end in .bam
				print(f"One of the .bam filepaths is invalid. \nSpecifically:\n {line} \nPlease make sure it ends in .bam")
				exit()
			samples.append(line)
ic(samples)

base_cov0 = pd.DataFrame()
base_cov1 = pd.DataFrame()
base_covAvg = pd.DataFrame()
base_data = pd.DataFrame()

for sample in samples:
	ic(sample)
	CheckFeatureCoverage(sample)

ic(base_cov0)

#base_cov0.to_csv(f'CombinedResults/ExonCovOver0.csv',index=False)
#base_cov1.to_csv(f'CombinedResults/ExonCovOver1.csv',index=False)
#base_covAvg.to_csv(f'CombinedResults/ExonCovAvg.csv',index=False)

"""




output_data = []
#ic(output_data)
df = pd.DataFrame(output_data)
ic(df)
ic(output_filepath)
output_file_name = ""
if output_filepath == "":
	print("No output filepath detected, using input filename with .bam replaced by suffix.tsv")
	output_file_name = input_file.replace(".bam",".") + os.path.basename(input_gff)[-35:]+ ".tsv"
else:
	output_file_name = output_filepath
ic(output_filepath)
df.to_csv(output_file_name,index=False,sep='\t')

exit()
"""
