#!/usr/local/opt/python/bin/python
#author: Tobias Hofmann, tobiashofmann@gmx.net
#workflow inspired by Yann Bertrand

#_____________________________________________________________________________________
#%%% Imports %%%
import os
import sys
import glob
import shutil
import random
import argparse
import fileinput
from cogent import LoadSeqs, DNA
from Bio import AlignIO
from Bio.Alphabet import IUPAC, Gapped

# Complete path function
class CompletePath(argparse.Action):
    """give the full path of an input file/folder"""
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, os.path.abspath(os.path.expanduser(values)))
		
# Get arguments
def get_args():
	parser = argparse.ArgumentParser(
		description="Clean and trim raw Illumina read files",
		formatter_class=argparse.ArgumentDefaultsHelpFormatter
	)
	parser.add_argument(
		'--input',
		required=True,
		action=CompletePath,
		default=None,
		help='The directory containing all fasta-alignment files'
	)
	parser.add_argument(
		'--config',
		action=CompletePath,
		help='A configuration file containing all taxa names to be used for SNP extraction'
	)
	parser.add_argument(
		'--phased',
		action='store_true',
		default=False,
		help='Use flag if alignments contain phased sequences'
	)
	parser.add_argument(
		'--output',
		required=True,
		action=CompletePath,
		default=None,
		help='The output directory where results will be saved'
	)
	return parser.parse_args()

# Get arguments
args = get_args()
# Set working directory
work_dir = args.input
out_dir = args.output
# Create the output directory
if not os.path.exists(out_dir):
    os.makedirs(out_dir)
config = args.config
name_file = open (config, 'r')
taxa_list = name_file.readlines()
taxa_names = []
for element in taxa_list:
	taxa_names.append(element.replace("\n", ""))

print "\n\n"
print " ____________________________________________________"
print "|Launching SNP extraction script...                  |"
print "|                                                    |"
print "|Written by Tobias Hofmann, inspired by Yann Bertrand|"
print "|Version 1.1, August 2015                            |"
print "|____________________________________________________|"
if not args.phased:
	print "\n\nScript is treating data as unphased alignment (add flag --phased to command if your data is phased)\n\n"
else:
	delimiter = raw_input("\n\n**************************************************\nUSER INPUT REQUIRED\n**************************************************\n\nYou indicated in your command that the alignment contains phased data, i.e. multiple alleles per sample.\nWhat is the delimiter that separates the different alleles from the sample name in the alignment?\nExample: If your alignment contains two alleles for sample1 which are named sample1_allele1 and sample1_allele2, the delimiter would be _\n\nPlease enter your delimiter: ")

#_______________________________________________________________________________
#%%% Functions %%%
def find_names(sequence_names):
	names = []
	for element in sequence_names:
		name, allele = element.rsplit(delimiter, 1)
		names.append(name)
	return names


def variable_positions(alignment):
	var_col = []
	for x in range(len(alignment)):
		if alignment[x].filtered(lambda x: len(set(x)) == 2 and "n" not in x and "N" not in x and "-" not in x):
			var_col.append(x)
	return var_col

	
def unphased_snps(list_var):
	if len(list_var)>0:
		list_positive = list_var
	else:
		print("no snp extraction performed due too a lack of polymorphic sites")
		return None
	#chooses randomly one snp position and saves position-coordinate
	snp = random.sample(list_positive, 1)[0]
	print "sampling position", snp
	
	#creates an alignment with only the extracted position
	temp_snp_align = edited_alignment[snp]
	#creates dictionary from the extracted snp position
	seq_dict = temp_snp_align.todict()
	#creates a set from the dictionary, with ordered characters (in order to replace them properly)
	set_values = list(set(seq_dict.values()))
	#code the base-letters into 0 or 1
	if "A" or "C" or "T" or "G" or "a" or "c" or "t" or "g" in set_values:
		zero = set_values[0]
		one = set_values[1]
		snp_dict = {}
		for name_seq, nucleotide in seq_dict.items():
			if nucleotide == zero:
				snp_dict[name_seq] = "0"
			else:
				snp_dict[name_seq] = "1"
	return snp_dict	
	
	
def phased_snps(list_var):
	if len(list_var)>0:
		list_positive = list_var
	else:
		print("no SNP extraction performed due too a lack of polymorphic sites")
		return None
	#chooses randomly one snp position and saves position-coordinate
	snp = random.sample(list_positive, 1)[0]
	print "sampling position", snp
	#creates an alignment with only the extracted position
	temp_snp_align = edited_alignment[snp]
	#creates dictionary from the extracted snp position
	seq_dict = temp_snp_align.todict()
	#finds elements in the dictionary that belong to the same sample (in case of multiple alleles)
	clean_names = sorted(find_names(taxa_names))
	no_duplicates = list(set(clean_names))
	new_dict = {taxon : [value for key, value in seq_dict.items() if key.startswith(taxon)] for taxon in no_duplicates}
	#returns which two bases are present in the dictionary (ordered, to replace them properly)
	set_values = list(set(seq_dict.values()))
	# Code the base-letters into 0, 1 or 2
	if "A" or "C" or "T" or "G" or "a" or "c" or "t" or "g" in set_values:
		zero = set_values[0]
		two = set_values[1]
		snp_dict = {}
		for name_seq, nucleotide in new_dict.items():
			if nucleotide[0] == nucleotide[1] == zero:
				snp_dict[name_seq] = "0"
			elif nucleotide[0] == nucleotide[1] == two:
				snp_dict[name_seq] = "2"
			elif nucleotide[0] != nucleotide[1]:
				snp_dict[name_seq] = "1"
	return snp_dict

	
#_____________________________________________________________________________________
#%%% Workflow %%%
#get all fasta files in working directory
os.chdir(work_dir)
fasta_files = [x for x in os.listdir(os.getcwd()) if ".fasta" in x and "snp.fasta" not in x]
final_dict = {}
list_of_genes=[]
for fasta in fasta_files[:]:
	print "_" * 50
	print "processing file:", fasta
	aln = LoadSeqs(fasta, moltype=DNA)
	# Apply filter of user-set taxa names to be used for snp-extraction
	edited_alignment = aln.takeSeqs(sorted(taxa_names))
	# Get the variable positions for each fasta file
	var_pos_list = variable_positions(edited_alignment)
	print var_pos_list
	D = ""
	if not args.phased:
		D = unphased_snps(var_pos_list)
	else:
		D = phased_snps(var_pos_list)
	if D != None:
		list_of_genes.append(fasta)
		for key, values in D.items():
			if key not in final_dict.keys():
				final_dict[key] = [D[key]]
			else:
				final_dict[key].append(D[key])

#join all snps into one dictionary
final_snp_alignment = {}
for key, value in final_dict.items():
	final_snp_alignment[key] = "".join(value)

# Create the output file in output directory
output_file_fasta = os.path.join(out_dir,'snp.fasta')
#print the snp dictionary into a fasta-file
with open(output_file_fasta, "wb") as f:
	for k, v in final_snp_alignment.items():
		f.write(">" + k+ "\n")
		f.write(v+ "\n")

# Create output file for SNAPP
output_file_nexus = os.path.join(out_dir,'snp.nexus')
aln = AlignIO.read(open(output_file_fasta), "fasta", alphabet=Gapped(IUPAC.ambiguous_dna))
with open(output_file_nexus, "wb") as n:
	n.write(aln.format("nexus"))
if not args.phased:
	for line in fileinput.input(output_file_nexus, inplace = 1):
		print line.replace("format datatype=dna missing=? gap=-;", "format datatype=binary symbols=01;").rstrip()
else:
	for line in fileinput.input(output_file_nexus, inplace = 1):
		print line.replace("format datatype=dna missing=? gap=-;", "format datatype=integerdata symbols=\"012\";").rstrip()
