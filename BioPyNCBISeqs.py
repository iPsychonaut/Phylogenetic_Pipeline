# -*- coding: utf-8 -*-
"""
Created on Wed Nov 30 11:44:25 2022

@author: theda
"""

# Download Sequences from NCBI and combine into a single file for alingment

from Bio import Entrez
from Bio import SeqIO

# What Databases can Entrez view?
# handle = Entrez.einfo()
# rec = Entrez.read(handle)
# handle.close()
# print(rec.keys())
# rec['DbList']

def search_ncbi(user_email, search_term, return_number, search_db, save_path):
    # Set user Email for Entrez
    Entrez.email = user_email
    # Search a Database (db) for a specific term
    print(search_db)
    handle = Entrez.esearch(db = search_db, term = search_term, retmax = return_number)
    rec_list = Entrez.read(handle)
    handle.close()
    #print(rec_list['Count']) # Print the Number of total items returned
    #print(len(rec_list['IdList'])) # Print the Number Returned by Search
    #print(rec_list['IdList']) # Print the Id's of each item in the Returned Search
    
    id_list = rec_list['IdList']
    handle = Entrez.efetch(db = search_db, id=id_list, rettype='gb') # Returns as Genbank Format that we need to parse with SeqIO
    recs = list(SeqIO.parse(handle, 'gb'))
    sequence_list = []
    description_list = []
    import_list = []
    for seq_record in recs:
        print(f"Dealing with GenBank record {seq_record.description}") 
        faa_filename = f'{save_path}/{seq_record.id}.fasta'
        import_list.append(seq_record)
        with open(faa_filename, 'w') as output_handle:
            output_handle.write(f"{seq_record.id} {seq_record.description}\n{seq_record.seq}\n")
            sequence_list.append(seq_record.id)
            description_list.append(seq_record.description)
            sequence_list.append(seq_record.seq)        
    handle.close()
    
        
    # Combine all of the individual sequences into a new file 
    combined_path = f"{save_path}/combined.fasta"
    print(combined_path)
    combined_seqs = SeqIO.write(import_list, combined_path, "fasta")
    
    return(combined_path)

# # Set user Email for Entrez
# user_email = 'ian.michael.bollinger@gmail.com'
# # Set search term
# search_term = 'CRT[Gene Name] AND "Plasmodium falciparum"[Organism]'

# returned_seqs = search_ncbi(user_email, search_term)
# print(returned_seqs)
