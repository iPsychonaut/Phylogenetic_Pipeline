# -*- coding: utf-8 -*-
"""
Created on Wed Nov 30 11:19:41 2022

Main Phylogenetic GUI

@author: ian.michael.bollinger@gmail.com
"""
###############################################################################
# IMPORT LIBRARIES
###############################################################################
from BioPyNCBISeqs import search_ncbi
from BioPyTreeGen import run_iqtree, tree_from_alignment, gen_boostrap_consensus_tree
from BioPySeqAlignments import MUSCLE_alignment, trimAI_alignment, fasta_to_strict_phylip, fasta_to_relaxed_phylip
from BioPyFolderSeqs import combine_fastas
import os
from datetime import datetime

###############################################################################
# DEBUG WORKSPACE
###############################################################################
# Set max number of records to return
return_number = 100
# Set user Email for Entrez
user_email = 'ian.michael.bollinger@gmail.com'

# # Set fasta directory
# save_path = 'TEST IMPORT'
# # Generate Combined Sequence List
# combined_path = combine_fastas(save_path)

# Set search term
#PsiM, PsiD, PsiK, PsiT2, PsiT1, 
gene_name = 'Internal transcribed spacer'.replace(' ','_')
organism = 'Psilocybe'.replace(' ','_')
search_term = f'"{gene_name}"[All Fields] AND "{organism}"[Organism]'
search_db = 'nucleotide'
# datetime object containing current date and time
now = datetime.now()
# dd/mm/YY H:M:S
dt_string = now.strftime("%Y%m%d_%H%M")
save_path = f'{gene_name}_{organism}_{return_number}_{dt_string}'

# Check whether the specified path exists or not
isExist = os.path.exists(save_path)
if not isExist:
    # Create a new directory because it does not exist
    os.makedirs(save_path)

# Download Sequences and Generate Combined Sequence List
combined_path = search_ncbi(user_email, search_term, return_number, search_db, save_path)



###############################################################################
# MAIN FUNCTIONS CALL
###############################################################################
# Convert Combile Sequence List into MUSCLE Alignment (output fasta)
MUSCLE_aln_path = MUSCLE_alignment(combined_path, save_path)

# Trim the MUSCLE Alignment and convert to Phylip (relaxed or strict)
trimmed_aln_path = trimAI_alignment(MUSCLE_aln_path)
relaxed_phylip_aln_path = fasta_to_relaxed_phylip(trimmed_aln_path)
#strict_phylip_aln_path = fasta_to_strict_phylip(trimmed_aln_path)

# CONSIDER OTHER ALIGNMENTS(?): MAFFT, CLUSTALW

# Creating Tree Files and Figures based on MUSCLE Alignment
run_iqtree(relaxed_phylip_aln_path, save_path)
tree_from_alignment(relaxed_phylip_aln_path, 'upgma')
tree_from_alignment(relaxed_phylip_aln_path, 'nj')
gen_boostrap_consensus_tree(relaxed_phylip_aln_path, 500)
