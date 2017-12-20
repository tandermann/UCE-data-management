#!/usr/local/opt/python/bin/python
#author: Tobias Hofmann, tobiashofmann@gmx.net

import os
import sys
import re
import glob
import shutil
import argparse
import csv
import random


#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# Input

# Complete path function
class CompletePath(argparse.Action):
    """give the full path of an input file/folder"""
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, os.path.abspath(os.path.expanduser(values)))

# Get arguments
def get_args():
    parser = argparse.ArgumentParser(
        description="This script will create chinmeric allele sequences from pairs of allele sequences (randomely shuffeling SNPs).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        '--input',
        required=True,
        action=CompletePath,
        default=None,
        help='The directory containing fasta alignments'
    )
    parser.add_argument(
        '--output',
        required=True,
        action=CompletePath,
        default=None,
        help='The output directory where results will be safed'
    )

    return parser.parse_args()

    
# Get arguments
args = get_args()
# Set working directory
#work_dir = args.input
#out_dir = args.output

<<<<<<< HEAD
work_dir = '/Users/tobias/GitHub/topaza_uce/additional_analyses_2017/data/empirical/all_alignments/allele_alignments_topaza'
out_dir = '/Users/tobias/GitHub/topaza_uce/additional_analyses_2017/data/empirical/all_alignments/chimeric_allele_alignments_topaza'
=======
work_dir = '/Users/tobias/GitHub/topaza_uce/additional_analyses_2017/data/empirical/allele_alignments_topaza'
out_dir = '/Users/tobias/GitHub/topaza_uce/additional_analyses_2017/data/empirical/chimeric_allele_alignments_topaza'
>>>>>>> 037dae023d75677fcd763fbc0af39156e2b836c2

if not os.path.exists(out_dir):
    os.makedirs(out_dir)


#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# Functions
    
    
def read_fasta(fasta):
    name, seq = None, []
    for line in fasta:
        line = line.rstrip()
        if line.startswith(">"):
            if name: yield (name, ''.join(seq))
            name, seq = line, []
        else:
            seq.append(line)
    if name: yield (name, ''.join(seq))

        
        
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# Workflow  


# Create a list of all fasta files
fasta_files = []    
for fasta in os.listdir(work_dir):
    if fasta.endswith(".fasta") or fasta.endswith(".fa"):
        fasta_files.append(fasta)

#fasta = fasta_files[1]

for fasta in fasta_files:
    # Create an output consensus fasta file for each allele alignment
    fasta_cons_name = re.sub("phased","chimeric_allele",fasta)
    out_fasta = open(os.path.join(out_dir, fasta_cons_name), 'w')
    
    # Create a dictionary for each fasta file, where both allele sequences are assigned to the same key for each sample
    seq_dict = {}
    with open("%s/%s" %(work_dir,fasta)) as f:
       for name, seq in read_fasta(f):
         name = re.sub('>', '', name)
         species = re.sub('_[0,1]', '', name)
         #name = key
         seq_dict.setdefault(species,[])
         seq_dict[species].append(seq)
    # Create a chimeric allele dict for each fasta file with the correct new header name as key and the consensus sequence of the two alleles as value
    chimeric_dict = {}
    for header in seq_dict:
        sequence = seq_dict[header]
        allele0 = sequence[0]
        allele1 = sequence[1]
        new_chimeric_allele_0_header = '%s_0' %header
        new_chimeric_allele_1_header = '%s_1' %header
        chimeric_dict.setdefault(new_chimeric_allele_0_header,[])
        chimeric_dict.setdefault(new_chimeric_allele_1_header,[])
        # Find those positions where the two alleles differ from each other and make a random pick of one of the versions, simulationg a consensus sequence
        for id, base in enumerate(allele0):
<<<<<<< HEAD
            if base == 'N' or allele1[id] == 'N':
                chimeric_dict[new_chimeric_allele_0_header].append(base)
                chimeric_dict[new_chimeric_allele_1_header].append(allele1[id])
            elif base != allele1[id]:
=======
            if base != allele1[id] and base != 'N':
>>>>>>> 037dae023d75677fcd763fbc0af39156e2b836c2
                variation = [base,allele1[id]]
                snp_new_chimeric_allele_0 = random.choice(variation)
                variation.remove(snp_new_chimeric_allele_0)
                snp_new_chimeric_allele_1 = variation[0]
                chimeric_dict[new_chimeric_allele_0_header].append(snp_new_chimeric_allele_0)
                chimeric_dict[new_chimeric_allele_1_header].append(snp_new_chimeric_allele_1)                
            else:
                chimeric_dict[new_chimeric_allele_0_header].append(base)
                chimeric_dict[new_chimeric_allele_1_header].append(allele1[id])
            
    # Write the consensus dictionary into a fasta output file
    for chim_header in chimeric_dict:
        chim_sequence =  "".join(chimeric_dict[chim_header])
        chim_header = ">%s" %chim_header
        out_fasta.write(chim_header+"\n")
        out_fasta.write(chim_sequence+"\n")
    
    out_fasta.close()
        
        

    
