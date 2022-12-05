# -*- coding: utf-8 -*-
"""
Created on Fri Dec  2 09:51:01 2022

Generate a Multiple Sequence Alignment (MUSLCE) from a series of Sequences

@author: ian.michael.bollinger@gmail.com
"""
###############################################################################
# IMPORT LIBRARIES
###############################################################################
import os
from Bio import AlignIO, SeqIO

# Specify the location of your MUSCLE exe file
muscle_exe = r"C:/Users/theda/.spyder-py3/BioPy/muscle5.1.win64.exe" 
# Specify the location of your TrimAI exe file
trimAI_exe = r"C:/Users/theda/.spyder-py3/BioPy/trimal.exe"

###############################################################################
# ALIGNMENT FUNCTIONS
###############################################################################
# Functiont to genreate a MUSCLE Alignment
def MUSCLE_alignment(input_alignment, save_path):   
    
    # Build out MUSCLE Command Line
    MUSCLE_alignment_path = f"{save_path}/combined_MUSCLE.fasta"
    muscle_command = f'{muscle_exe} -align "{input_alignment}" -output "{MUSCLE_alignment_path}"'
    print(muscle_command)
    os.system(muscle_command)
    
    # Display Output Alignment
    # MUSCLE_alignment = AlignIO.read(MUSCLE_alignment_path, "fasta") 
    # print(MUSCLE_alignment)
    return(MUSCLE_alignment_path)

###############################################################################
# ALIGNMENT CLEAN-UP FUNCTIONS
###############################################################################
# Function to trim an alignment 
def trimAI_alignment(input_alignment):   
    
    # Build out trimAI Command Line
    output_trimAI_path = f"{input_alignment.replace('.fasta', '_trimAI.fasta')}"
    trimAI_command = f'{trimAI_exe} -in "{input_alignment}" -out "{output_trimAI_path}" -automated1'
    print(trimAI_command)
    os.system(trimAI_command)
    return(output_trimAI_path)

# Function to Convert a combined fasta alignment into a relaxed phylip format
def fasta_to_relaxed_phylip(input_alignment):
    output_phylip_path = f"{input_alignment.replace('.fasta', '_relaxed.phylip')}"
    AlignIO.convert(input_alignment, 'fasta', output_phylip_path, 'phylip-relaxed')
    return(output_phylip_path)

# Function to Convert a combined fasta alignment into a strict phylip format
def fasta_to_strict_phylip(input_alignment):
    import_list = []
    with open(input_alignment) as handle:
        for seq_record in SeqIO.FastaIO.FastaIterator(handle):        
            temp_id = seq_record.id.split('_')[2]
            seq_record.id = temp_id
            print(f'STRICT PHYLIP ID: {seq_record.id}')
            print(f'STRICT PHYLIP DESCRIPTION: {seq_record.description}')
            import_list.append(seq_record)
    output_phylip_path = f"{input_alignment.replace('.fasta', '_strict.phylip')}"
    SeqIO.write(import_list, output_phylip_path, "phylip")
    return(output_phylip_path)

###############################################################################
# DEBUG WORKSPACE
###############################################################################
# file_name = 'combined.fasta'
# save_path = 'C:/Users/theda/.spyder-py3/BioPy/TEST_IMPORT'
# # # save_path= 'C:/Users/theda/.spyder-py3/BioPy'
# # # file_name = 'neoxalapensis.fasta'
# input_alignment = f'{save_path}/{file_name}'

# MUSCLE_aln_path = MUSCLE_alignment(input_alignment, save_path)
# trimmed_aln_path = trimAI_alignment(MUSCLE_aln_path)
# # strict_phylip_aln_path = fasta_to_strict_phylip(trimmed_aln_path)
# relaxed_phylip_aln_path = fasta_to_relaxed_phylip(trimmed_aln_path)
# # #print(relaxed_phylip_aln_path)
