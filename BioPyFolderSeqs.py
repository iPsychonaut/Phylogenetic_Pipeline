# -*- coding: utf-8 -*-
"""
Created on Thu Dec  1 11:53:14 2022

Combine Fasta Files into a single file for alignment

@author: ian.michael.bollinger@gmail.com
"""
###############################################################################
# IMPORT LIBRARIES
###############################################################################
import os
from Bio import SeqIO 

###############################################################################
# FASTA FILE FROM FOLDER OF FASTA FILES COMBINER
###############################################################################
# Function to combine a folder of fasta files into a single fasta file
def combine_fastas(input_directory):
    # Iterate over files in given directory
    sequence_list = []
    description_list = []
    import_list = []
    for file_name in os.listdir(input_directory):
        file_path = f'{input_directory}/{file_name}'
        # Checking if it is a file
        if os.path.isfile(file_path):
            #print(file_path)
            with open(file_path) as handle:
                for seq_record in SeqIO.FastaIO.FastaIterator(handle):        
                    temp_desc = seq_record.id
                    #temp_id = seq_record.description.replace(f'{temp_desc} ', '').split('voucher')[0:2]
                    temp_id = f'{temp_desc}'.replace('  ', ' ').replace(' ', '_')
                    seq_record.id = temp_id
                    seq_record.description = temp_desc
                    print(f'FINAL ID: {seq_record.id}')
                    print(f'FINAL DESCRIPTION: {seq_record.description}')
                    sequence_list.append(seq_record.id)
                    description_list.append(seq_record.description)
                    import_list.append(seq_record)
    # Combine all of the individual sequences into a new file 
    combined_path = f"{input_directory}/combined.fasta"
    print(combined_path)
    SeqIO.write(import_list, combined_path, "fasta")
    return(combined_path)

###############################################################################
# DEBUG WORKSPACE
###############################################################################
# combine_fastas('C:/Users/theda/.spyder-py3/BioPy/TEST IMPORT')
