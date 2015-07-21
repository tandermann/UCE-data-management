#!/usr/local/opt/python/bin/python
#author: Tobias Hofmann, tobiashofmann@gmx.net


#_____________________________________________________________________________________
#%%% Imports %%%
import os, sys, shutil, time, random, getopt
from cogent import LoadSeqs, DNA


#_____________________________________________________________________________________
#%%% Load files %%%
# read the argument after the script name (has to be given by user per commandline) and assign it to string "path"
path = sys.argv[1]

# Get the full path of the target file
full_path = os.path.abspath(path)

#define the working directory
work_dir = full_path[:full_path.rfind("/"):]

# read in the fasta file as a DNA sequence
DNA_seq = LoadSeqs(full_path, moltype=DNA)

#here the user can specify which sequences/taxa should be included in the SNP extraction
taxa_names = [x for x in DNA_seq.Names if x == "T_pyra1_0" or x == "T_pyra2_0" or x == "T_pyra3_0" or x == "T_pyra4_0" or x == "T_pella5_0" or x == "T_pella6_0" or x == "T_pella7_0" or x == "T_pella8_0" or x == "T_pella9_0" or x == "T_pyra1_1" or x == "T_pyra2_1" or x == "T_pyra3_1" or x == "T_pyra4_1" or x == "T_pella5_1" or x == "T_pella6_1" or x == "T_pella7_1" or x == "T_pella8_1" or x == "T_pella9_1"]


#_____________________________________________________________________________________
#%%% Functions %%%
def variable_positions(alignment):
	var_col = []
	for x in range(len(alignment)):
		if alignment[x].filtered(lambda x: len(set(x)) == 2 and "n" not in x and "N" not in x and "-" not in x):
			var_col.append(x)
	return var_col

	
#This function can be used to extract one variable position per alignment, need list of variable positions as input
def extract(list_var):
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
	# read in the fasta file as a DNA sequence
	aln = LoadSeqs(fasta, moltype=DNA)
	#applying the user-named taxa names that shall be used for extraction
	edited_alignment = aln.takeSeqs(taxa_names)
	var_pos_list = variable_positions(edited_alignment)
	print var_pos_list
	D = extract(var_pos_list)
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

#print the snp dictionary into a fasta-file
with open("snp.fasta", "wb") as f:
	for k, v in final_snp_alignment.items():
		f.write(">" + k+ "\n")
		f.write(v+ "\n")
