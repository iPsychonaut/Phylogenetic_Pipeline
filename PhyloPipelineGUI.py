# -*- coding: utf-8 -*-
"""
Created on Wed Nov 30 11:19:41 2022

GUI Phylogenetic Pipeline

@author: ian.michael.bollinger@gmail.com
"""
###############################################################################
# IMPORT LIBRARIES
###############################################################################
from Bio import Phylo, AlignIO, SeqIO, Entrez
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.Phylo.Consensus import majority_consensus, bootstrap_consensus, get_support
import tkinter as tk
from tkinter import filedialog, StringVar, IntVar
import tkinter.font as font
from PIL import ImageTk, Image
import os
from datetime import datetime
import subprocess

# Get the path to the user's home directory
home_path = os.path.expanduser('~')
base_path = home_path.split('Users')[0] 
install_path = f"{base_path}PhyloPipeline"
bin_path = f"{install_path}/bin"

# Set Executable paths 
iqtree_path  = f'{bin_path}/iqtree-1.6.12-Windows/bin/iqtree.exe'
muscle_exe = f'{bin_path}/muscle.exe'
trimAI_exe = f'{bin_path}/trimal.exe'
figtree_exe = f'{bin_path}/FigTree.exe'

# Set Distance Calculator to 'identity' for Neucleotide Processing
distance_calculator = DistanceCalculator('identity')

###############################################################################
# GUI Functions Workspace
###############################################################################

def display_update(msg):
    search_entry_display.config(text=msg)
    folder_entry_display.config(text=msg)
    alignment_entry_display.config(text=msg)  

# Function to have user to select a directory and store it in global var called folder_path
def update_folder_param():
    global folder_path
    global save_path
    filename = filedialog.askdirectory(title = 'Select a folder')
    folder_path.set(filename)
    save_path = folder_path.get() # NEEDS UPDATE
    msg = f'Saving alignments and trees to\n{save_path}'
    display_update(msg)
    return(save_path, msg)

# Function to run the Folder Compiler Pipeline
def run_folder_pipeline():
    global save_path
    display_update('Processing...\nLarger Data Sets take longer')
    save_path, msg = update_folder_param()
    
    # Generate Combined Sequence List
    combined_path = combine_fastas(save_path)
    run_main_pipeline(combined_path, msg)

# Function to run the NCBI Search/Compiler Pipeline
def run_search_pipeline():
    global save_path
    global search_email
    global search_gene_abrv
    global search_organism
    global return_count
    
    display_update('Processing...\nLarger Data Sets take longer')
    
    search_email = email_entry.get()
    search_gene_name = gene_name_entry.get()
    search_gene_abrv = gene_abrv_entry.get()
    search_organism = organism_entry.get()
    return_count = return_count_entry.get()
    search_minbp = 100
    search_maxbp = 5000
    
    now = datetime.now() # datetime object containing current date and time
    dt_string = now.strftime("%Y%m%d_%H%M")
    save_path = f'{search_gene_name.replace("*","(WC)").replace(" ","_")}_{search_gene_abrv.replace("*","(WC)").replace(" ","_")}_{search_organism.replace("*","(WC)").replace(" ","_")}_{return_count}_{dt_string}'

    search_term = f"""
(({search_gene_abrv}[Gene Name] OR "{search_gene_name}"[All Fields])) AND ("{search_organism}"[Organism] OR ("{search_organism}"[Organism] OR {search_organism}[All Fields])) AND ({search_minbp}[SLEN] : {search_maxbp}[SLEN]) NOT "whole genome shotgun"[All Fields]
"""
    search_db = 'nucleotide'
    
    if search_email == '':
        search_entry_display.config(text='Please add an Email')
    elif return_count == '':
        search_entry_display.config(text='Please set the Max Return Records Count')
    elif search_gene_abrv == '' and search_organism == '':
        search_entry_display.config(text='Please enter a Gene AND/OR an Organism')
    else:
        msg = f"""
EMAIL: {search_email}\nGENE: {search_gene_abrv}
ORGNAISM: {search_organism}
Saving up to {return_count}
NCBI records (between {search_minbp}-{search_maxbp}bp),
alignments, and trees to
{save_path}
"""
        # Check whether the specified path exists or not
        isExist = os.path.exists(save_path)
        if not isExist:
            # Create a new directory because it does not exist
            os.makedirs(save_path)
    
        # Download Sequences and Generate Combined Sequence List
        combined_path = search_ncbi(search_email, search_term, return_count, search_db, save_path)
        
        # Check if combined Fasta is empty
        print(len(combined_path))
        if combined_path[0] == '':
            msg = f'No sequences were found in {combined_path}\nRetry search with new parameters'
            display_update(msg)
        else:
            run_main_pipeline(combined_path, msg)           

# Function to run the Alignment Pipeline on a trimmed sequence
def run_alignment_pipeline():
    global save_path
    display_update('Processing...\nLarger Data Sets take longer')
    file_types = (('FASTA files', '*.fasta'),
                 ('Phylip files', '*.phylip'),
                 ('All files', '*.*'))
    combined_path = filedialog.askopenfilename(
                 title='Please load a Multiple Sequence Alignment OR Sequence List',
                 initialdir='./',
                 filetypes=file_types)
    file_name = combined_path.split('/')[-1]
    save_path = combined_path.strip(file_name)
    msg = f"""
Alignments and/or trees were generated from
{combined_path}
and saved to
{save_path}
"""    
    run_main_pipeline(combined_path, msg)

# Function to run the Tree Generators on a trimmed alignment
def run_tree_generators(trimmed_aln_path, save_path, msg):
    display_update(msg)
    # Creating Tree Files and Figures based on MUSCLE Alignment
    gen_nj_tree(trimmed_aln_path)
    gen_boostrap_consensus_tree(trimmed_aln_path, save_path, 200)
    run_iqtree(trimmed_aln_path)

# Function to take the funneled combined fastas and run then through main pipeline
def run_main_pipeline(combined_path, msg):
    # Convert Compile Sequence List into MUSCLE Alignment (output fasta)
    MUSCLE_aln_path = MUSCLE_alignment(combined_path, save_path)

    # Trim the MUSCLE Alignment # and convert to Phylip (relaxed or strict)
    trimmed_aln_path = trimAI_alignment(MUSCLE_aln_path)
    
    # Run Tree Generators
    run_tree_generators(trimmed_aln_path, save_path, msg)
    
# Function to swap from the search frame to the folder frame.
def change_to_folder():
    search_frame.forget()
    alignment_frame.forget()
    folder_frame.pack(fill='both', expand=1)

# Function to swap from the folder frame to the search frame
def change_to_search():
    folder_frame.forget()
    alignment_frame.forget()
    search_frame.pack(fill='both', expand=1)

# Function to swap from the folder frame to the search frame
def change_to_alignment():
    folder_frame.forget()
    search_frame.forget()
    alignment_frame.pack(fill='both', expand=1)
    
###############################################################################
# COMPILER FUNCTIONS
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
            with open(file_path) as handle:
                for seq_record in SeqIO.FastaIO.FastaIterator(handle):        
                    temp_desc = seq_record.id
                    temp_id = seq_record.description[0:40] # Longest organism name is 36 characters
                    seq_record.id = temp_id
                    seq_record.description = temp_desc
                    seq_record.seq = seq_record.seq.replace('.','-')
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

# Function to search NCBI for a given Gene and Organism and return the maximum number of records
def search_ncbi(user_email, search_term, return_number, search_db, save_path):
    
    # Set user Email for Entrez
    Entrez.email = user_email
    
    # Search a Database (db) for a specific term
    print(search_db)
    handle = Entrez.esearch(db = search_db, term = search_term, retmax = return_number)
    rec_list = Entrez.read(handle)
    handle.close()
    print(rec_list['Count']) # Print the Number of total items returned
    print(len(rec_list['IdList'])) # Print the Number Returned by Search
    print(rec_list['IdList']) # Print the Id's of each item in the Returned Search  
    id_list = rec_list['IdList']
    # Returns as Genbank Format that we need to parse with SeqIO
    handle = Entrez.efetch(db = search_db, id=id_list, rettype='gb') 
    recs = list(SeqIO.parse(handle, 'gb'))
    sequence_list = []
    description_list = []
    import_list = []
    for seq_record in recs:
        print(f"Dealing with GenBank record {seq_record.description}") 
        faa_filename = f'{save_path}/{seq_record.id}.fasta'
        import_list.append(seq_record)
        with open(faa_filename, 'w') as output_handle:
            try:
                output_handle.write(f'>{seq_record.id} {seq_record.description}\n{seq_record.seq}\n')
                sequence_list.append(seq_record.id)
                description_list.append(seq_record.description)
                seq_record.seq = seq_record.seq.replace('.','-')
                sequence_list.append(seq_record.seq)        
            except:
                print(f'PASS DUE TO UNDEFINED SEQUENCE ERROR: {faa_filename}')
    print(import_list)
    handle.close()    
        
    # Combine all of the individual sequences into a new file 
    combined_path = f"{save_path}/combined.fasta"
    print(combined_path)
    SeqIO.write(import_list, combined_path, "fasta")
    return(combined_path)

###############################################################################
# ALIGNMENT FUNCTIONS
###############################################################################

# Functiont to genreate a MUSCLE Alignment
def MUSCLE_alignment(input_alignment, save_path):   
    
    # Build out MUSCLE Command Line
    MUSCLE_alignment_path = f"{save_path}/combined_MUSCLE.fasta"
    subprocess.run([muscle_exe, '-align', input_alignment, '-output', MUSCLE_alignment_path])
    # Display Output Alignment
    # MUSCLE_alignment = AlignIO.read(MUSCLE_alignment_path, "fasta") 
    # print(MUSCLE_alignment)
    return(MUSCLE_alignment_path)

# Function to trim an alignment 
def trimAI_alignment(input_alignment):   
    
    # Build out trimAI Command Line
    output_trimAI_path = f"{input_alignment.replace('.fasta', '_trimAI.fasta')}"
    subprocess.run([trimAI_exe, '-in', input_alignment, '-out', output_trimAI_path, '-automated1'])
    return(output_trimAI_path)

###############################################################################
# TREE FUNCTIONS
###############################################################################

# Function to Generate Trees from an input alignment
def run_iqtree(input_alignment_path):
    
    # Generate IQ-Tree Files
    #os.system(f'{iqtree_path} -s "{input_alignment_path}" -nt AUTO')
    subprocess.run([iqtree_path, '-s', input_alignment_path, '-nt', 'AUTO'])
    
    # Generate and attach Branch Support to the Fina IQ-Tree
    output_tree_path = f'{input_alignment_path}.treefile'
    supported_tree = get_branch_support(output_tree_path, 'newick')
    output_supported_tree_path = f'{output_tree_path}_supported.tre'
    print(f'Saving: {output_supported_tree_path}')
    Phylo.write(supported_tree, output_supported_tree_path, 'newick')    

# Function to generate a Neighbor Joining from an input alignment
def gen_nj_tree(input_alignment_path):
    global distance_calculator
    
    # Open the alignment file as a MultipleSeqAlignment object
    with open(input_alignment_path,'r') as aln:
        working_alignment = AlignIO.read(aln, 'fasta')
    distance_matrix = distance_calculator.get_distance(working_alignment)
    
    # Open and initiate the appropriate Tree Constructor
    constructor = DistanceTreeConstructor()      
    output_tree = constructor.nj(distance_matrix)
    output_tree.rooted = True
    output_tree_path = input_alignment_path.replace('.fasta','_nj.tre')
    print(f'Saving: {output_tree_path}')
    Phylo.write(output_tree, output_tree_path, 'newick')
    
    # Generate and attach Branch Support to the Tree
    supported_tree = get_branch_support(output_tree_path, 'newick')
    output_supported_tree_path = output_tree_path.replace('_nj.tre','_nj_supported.tre')
    print(f'Saving: {output_supported_tree_path}')
    Phylo.write(supported_tree, output_supported_tree_path, 'newick')    

# Function to generate Bootstrap Trees from an input alignment
def gen_boostrap_consensus_tree(input_alignment_path, save_path, replicate_count):
    print(f'Processing {input_alignment_path} with {replicate_count}x Replicates')
    print('NOTE: THIS CAN TAKE A WHILE')
    
    # Open the alignment file as a MultipleSeqAlignment object 
    with open(input_alignment_path,'r') as aln:
        working_alignment = AlignIO.read(aln, 'fasta')
    global distance_calculator
    boostrap_constructor = DistanceTreeConstructor(distance_calculator)
    output_tree_path = input_alignment_path.replace('.fasta','_bootstrap.tre')
    bootstrap_consensus_tree = bootstrap_consensus(working_alignment, replicate_count,
                                                   boostrap_constructor, majority_consensus)
    print(f'Saving: {output_tree_path}')
    Phylo.write(bootstrap_consensus_tree, output_tree_path, 'newick')
    
    # Generate and attach Branch Support to the Consensus Tree
    supported_bootstrap_tree = get_branch_support(output_tree_path, 'newick')
    output_supported_tree_path = output_tree_path.replace('.tre','_supported.tre')
    print(f'Saving: {output_supported_tree_path}')
    Phylo.write(supported_bootstrap_tree, output_supported_tree_path, 'newick')    

# Function to Get Branch Support For Specific Tree
def get_branch_support(input_tree_path, input_tree_format):
    trees = list(Phylo.parse(input_tree_path, input_tree_format))
    target_tree = trees[0]
    support_tree = get_support(target_tree, trees)
    return(support_tree)
    
###############################################################################
# MAIN GUI CONSTRUCTION
###############################################################################

# Run main guard loop
if __name__ == '__main__':
    root= tk.Tk()
    root.title('PhyloPipeline')
    root.iconbitmap(f'{bin_path}/DNA_icon.ico')
    
    # Set placeholders for user variables
    search_email = StringVar() 
    search_gene_abrv = StringVar()
    search_organism = StringVar()
    return_count = IntVar()
    folder_path = StringVar()
    save_path = StringVar()

    # Generate Main root Window and load Frames
    width, height = 1000, 500
    screen_width = root.winfo_screenwidth()
    screen_height = root.winfo_screenheight()
    search_frame = tk.Frame(root)
    folder_frame = tk.Frame(root)
    alignment_frame = tk.Frame(root)
    img_logo = ImageTk.PhotoImage(Image.open(f'{bin_path}/PhyloPipelineLogo.png'))
    
    # Let's create the fonts that we need.
    font_large = font.Font(family = 'Avenir', size = '12', weight = 'bold')
    font_small = font.Font(family = 'Avenir', size = '12')
    font_button = font.Font(family = 'Avenir', size = '16', weight = 'bold')
    
    # Build out Folder Frame
    folder_frame_logo = tk.Label(folder_frame, image=img_logo).pack(pady = 5)
    folder_frame_heading = tk.Label(folder_frame, text='FOLDER COMPILER',
                                font=font_large).pack(pady = 5)    
    change_to_search_button = tk.Button(folder_frame, font=font_small,
                                   text='Change to NCBI Search Compiler',
                                   command=change_to_search).pack(pady = 5)
    change_to_alignment_button = tk.Button(folder_frame,
                                    text='Change to Sequence/MSA->Tree Builder',
                                    font = font_small,
                                    command = change_to_alignment).pack(pady = 5)
    folder_label = tk.Label(folder_frame,
                           text='You were prompted to select a folder containing\nONLY Fasta files for alignment',
                           font=font_small).pack(pady = 5)
    
    folder_entry_display = tk.Label(folder_frame,
                           text='No search submitted yet',
                           font=font_large,
                           bg='brown', fg='lightyellow')
    folder_entry_display.pack(pady = 5)
    
    folder_pipeline_button = tk.Button(folder_frame,
                             text='SELECT FOLDER & RUN PIPELINE',
                             fg='darkred', bg='darkgray', font = font_button,
                             command= run_folder_pipeline).pack(pady = 5)
    
    # Build out Search Frame
    search_frame_logo = tk.Label(search_frame, image = img_logo).pack(pady = 0)    
    change_to_folder_button = tk.Button(search_frame, text='Change to Folder Compiler',
                                   font = font_small, command = change_to_folder).pack(pady = 5)
    search_heading_label = tk.Label(search_frame, text = 'NCBI SEARCH COMPILER',
                                    font = font_large).pack(pady = 20)
    change_to_alignment_button = tk.Button(search_frame,
                                    text='Change to Sequence/MSA->Tree Builder',
                                    font = font_small,
                                    command = change_to_alignment).pack(pady = 5)
    email_label = tk.Label(search_frame,
                           text='Enter an email associated with an NCBI login\nEX: ian.michael.bollinger@gmail.com',
                           font=font_small).pack(pady=5)
    email_entry = tk.Entry(search_frame, font=font_small)
    email_entry.pack(pady=5)    
    gene_name_label = tk.Label(search_frame,
                           text='Enter the NAME of a gene to search\nEX: norbaeocystin methyltransferase, or Internal Transcribed Spacer',
                           font=font_small).pack(pady=5)
    gene_name_entry = tk.Entry(search_frame, font=font_small)
    gene_name_entry.pack(pady=5)
    gene_abrv_label = tk.Label(search_frame,
                           text='Enter the ABBREVIATION of the gene to search\nEX: PsiM, ITS*',
                           font=font_small).pack(pady=5)
    gene_abrv_entry = tk.Entry(search_frame, font=font_small)
    gene_abrv_entry.pack(pady=5)
    organism_label = tk.Label(search_frame,
                           text='Enter the scientific name of the organisms to search\nEX: Fungi, or Drosophila melanogaster',
                           font=font_small).pack(pady=5)
    organism_entry = tk.Entry(search_frame, font=font_small)
    organism_entry.pack(pady=5)
    return_count_label = tk.Label(search_frame,
                           text='Enter the maximum number of records to return',
                           font=font_small).pack(pady=5)
    return_count_entry = tk.Entry(search_frame, font=font_small)
    return_count_entry.pack(pady=5)
    search_entry_display = tk.Label(search_frame,
                           text='No search submitted yet',
                           font=font_large,
                           bg='brown', fg='lightyellow')
    search_entry_display.pack(pady=5)
    search_pipeline_button = tk.Button(search_frame,
                             text='RUN NCBI SEARCH & PIPELINE',
                             fg='darkred', bg='darkgray', font = font_button,
                             command=run_search_pipeline).pack(pady = 5)
    
    # Build out Alignment Frame
    search_frame_logo = tk.Label(alignment_frame, image = img_logo).pack(pady = 0)    
    change_to_folder_button = tk.Button(alignment_frame, text='Change to Folder Compiler',
                                   font = font_small, command = change_to_folder).pack(pady = 5)
    change_to_search_button = tk.Button(alignment_frame, font=font_small,
                                   text='Change to NCBI Search Compiler',
                                   command=change_to_search).pack(pady = 5)
    alignment_heading_label = tk.Label(alignment_frame, text = 'SEQUENCE/MSA->TREE BUILDER',
                                    font = font_large).pack(pady = 20)    
    alignment_entry_display = tk.Label(alignment_frame,
                           text='No sequences submitted yet',
                           font=font_large,
                           bg='brown', fg='lightyellow')
    alignment_entry_display.pack(pady=5)
    alignment_pipeline_button = tk.Button(alignment_frame,
                             text='SELECT SEQUENCES & RUN PIPELINE',
                             fg='darkred', bg='darkgray', font = font_button,
                             command=run_alignment_pipeline).pack(pady = 5)
    
    # Build the Folder Frame out and run it and Main Loop
    folder_frame.pack(fill='both', expand=1)
    root.mainloop()
