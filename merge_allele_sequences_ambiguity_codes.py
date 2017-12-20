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
work_dir = args.input
out_dir = args.output

work_dir = '/Users/tobias/GitHub/topaza_uce/additional_analyses_2017/data/empirical/allele_alignments_topaza'
out_dir = '/Users/tobias/GitHub/topaza_uce/additional_analyses_2017/data/empirical/iupac_consensus_alignments_topaza'
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

iupac_dict = dict([("['A', 'C']",'M'),("['C', 'A']",'M'),("['C', 'T']",'Y'),("['T', 'C']",'Y'),("['C', 'G']",'S'),("['G', 'C']",'S'),("['A', 'G']",'R'),("['G', 'A']",'R'),("['G', 'T']",'K'),("['T', 'G']",'K'),("['A', 'T']",'W'),("['T', 'A']",'W')])


for fasta in fasta_files:
    # Create an output consensus fasta file for each allele alignment
    fasta_cons_name = re.sub("allele","iupac_consensus",fasta)
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
    iupac_cons_dict = {}
    for header in seq_dict:
        sequence = seq_dict[header]
        allele0 = sequence[0]
        allele1 = sequence[1]
        iupac_cons_dict.setdefault(header,[])
        # Find those positions where the two alleles differ from each other and make a random pick of one of the versions, simulationg a consensus sequence
        for id, base in enumerate(allele0):
            if base == 'N' or allele1[id] == 'N':
                iupac_cons_dict[header].append(base)
            elif base == '-' or allele1[id] == '-':
                iupac_cons_dict[header].append('N')
            elif base != allele1[id]:
                variation = [base,allele1[id]]
                iupac_code  = iupac_dict[str(variation)]
                iupac_cons_dict[header].append(iupac_code)
            else:
                iupac_cons_dict[header].append(base)
            
    # Write the consensus dictionary into a fasta output file
    for cons_header in iupac_cons_dict:
        iupac_cons_sequence =  "".join(iupac_cons_dict[cons_header])
        cons_header = ">%s" %cons_header
        out_fasta.write(cons_header+"\n")
        out_fasta.write(iupac_cons_sequence+"\n")
    
    out_fasta.close()
        
        

    
