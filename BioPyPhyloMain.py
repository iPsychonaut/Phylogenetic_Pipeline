# -*- coding: utf-8 -*-
"""
Created on Wed Nov 30 11:19:41 2022

@author: theda
"""

# Import Necessary Libraries
from BioPyNCBISeqs import search_ncbi
from BioPyTreeFigures import trees_from_alignment
from BioPySeqMUSCLE import MUSCLE_alignment
from BioPyFolderSeqs import combine_fastas
import os
from datetime import datetime

# Set max number of records to return
return_number = 25
# Set user Email for Entrez
user_email = 'ian.michael.bollinger@gmail.com'

# Set fasta directory
save_path = 'TEST IMPORT'
combined_path = combine_fastas(save_path)

# # Set search term
# gene_name = 'Internal transcribed spacer'
# organism = 'Psilocybe'
# search_term = f'"{gene_name}"[All Fields] AND "{organism}"[Organism]'
# search_db = 'nucleotide'
# # datetime object containing current date and time
# now = datetime.now()
# # dd/mm/YY H:M:S
# dt_string = now.strftime("%d.%m.%Y_%H.%M")
# save_path = f'{gene_name}_{organism}_{return_number}_{dt_string}'

# # Check whether the specified path exists or not
# isExist = os.path.exists(save_path)
# if not isExist:
#    # Create a new directory because it does not exist
#    os.makedirs(save_path)

# # Download  Sequence List
# combined_path = search_ncbi(user_email, search_term, return_number, search_db, save_path)

# Combile Sequence List and convert into MUSCLE Alignment (output FASTA format)
muscle_alignment = MUSCLE_alignment(combined_path, save_path)

# Creating Tree Files and Figures based on MUSCLE Alignment
trees_from_alignment(muscle_alignment, save_path)