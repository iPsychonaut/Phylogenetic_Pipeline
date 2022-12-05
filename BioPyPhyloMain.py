# -*- coding: utf-8 -*-
"""
Created on Wed Nov 30 11:19:41 2022

Main Phylogenetic Program

@author: ian.michael.bollinger@gmail.com
"""
###############################################################################
# IMPORT LIBRARIES
###############################################################################
from BioPyNCBISeqs import search_ncbi
from BioPyTreeGen import run_iqtree, gen_nj_tree, gen_boostrap_consensus_tree
from BioPySeqAlign import MUSCLE_alignment, trimAI_alignment
from BioPyFolderSeqs import combine_fastas
import os
from datetime import datetime

###############################################################################
# DEBUG WORKSPACE
###############################################################################
# Pick one of the following:
    
# 1) Combined Fasta from Folder
# Set fasta directory
save_path = 'TEST_IMPORT'
# Generate Combined Sequence List
combined_path = combine_fastas(save_path)

# # 2) Combined Fasta from NCBI Search
# # Set search term
# user_email = 'ian.michael.bollinger@gmail.com'
# #PsiM, PsiD, PsiK, PsiT2, PsiT1, 
# gene_name = 'Internal transcribed spacer'.replace(' ','_')
# organism = 'Fungi'.replace(' ','_')
# return_number = 25
# search_minbp = 100
# search_maxbp = 5000
# search_term = f'"{gene_name}"[All Fields] AND "{organism}"[Organism] AND ("{search_minbp}"[SLEN] : "{search_maxbp}"[SLEN]) NOT "whole genome shotgun"[All Fields]'

# search_db = 'nucleotide'
# now = datetime.now() # datetime object containing current date and time
# dt_string = now.strftime("%Y%m%d_%H%M")
# save_path = f'{gene_name}_{organism}_{return_number}_{dt_string}'

# # Check whether the specified path exists or not
# isExist = os.path.exists(save_path)
# if not isExist:
#     # Create a new directory because it does not exist
#     os.makedirs(save_path)

# # Download Sequences and Generate Combined Sequence List
# combined_path = search_ncbi(user_email, search_term, return_number, search_db, save_path)

###############################################################################
# MAIN FUNCTIONS CALL
###############################################################################
# Convert Combile Sequence List into MUSCLE Alignment (output fasta)
MUSCLE_aln_path = MUSCLE_alignment(combined_path, save_path)

# Trim the MUSCLE Alignment and convert to Phylip (relaxed or strict)
trimmed_aln_path = trimAI_alignment(MUSCLE_aln_path)

# Creating Tree Files and Figures based on MUSCLE Alignment
gen_nj_tree(trimmed_aln_path)
gen_boostrap_consensus_tree(trimmed_aln_path, 200)
run_iqtree(trimmed_aln_path, save_path)
