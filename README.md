# Phylogenetic_Pipeline
### PURPOSE

Generate a MUSCLE Alignment based on either a Folder with Fasta files or NCBI Search Criteria and then output three (3) statistically supported trees: Neighbor Joining, Consensus Bootstrap, and the IQ-Tree output.


### REQUIREMENTS

IQ-Tree for Tree Building
1) Download the command line version it at: http://www.iqtree.org/
2) Update line 16 in BioPyTreeGen (or line 23 in PhyloPipelineGUI.py)  with the path of the IQ-Tree executable

MUSCLE for alingment
1) Download it at: https://2018-03-06-ibioic.readthedocs.io/en/latest/install_muscle.html
2) Update line 14 in BioPySeqAlign (or line 24 in PhyloPipelineGUI.py) with the path of the MUSCLE executable

TrimAI for alignment
1) Download the command line version it at: http://trimal.cgenomics.org/getting_started_with_trimal_v1.2
2) Update line 16 in BioPySeqAlign (or line 25 in PhyloPipelineGUI.py)  with the path of the TrimAI executable

### GUI DESCRIPTION

This will load a window where you can toggle between the three (3) options:
1) Compile as Folder with ONLY FASTA files in it
2) Search NCBI for a GENE AND ORGANISM and RETURN MAX RECORDS for compilation
3) Load an already generated sequence file or multiple sequence alignment (MSA) for tree generation

### FILE DESCRIPTIONS

-BioPyFolderSeqs contains the function combine_fastas(input_directory) which will combine all Fasta files in a given directory into a single fasta for alignment

-BioPyNCBISeqs contains the function search_ncbi(user_email, search_term, return_number, search_db, save_path) which will search NCBI for, download, and combine a desired number of Fasta files for alignment

-BioPySeqAlign contains the functions:
MUSCLE_alignment(combined_path, save_path) which will take a Combined Fasta and output a MUSCLE alignment;
trimAI_alignment(input_alignment) to call the AI trimmer

-BioPyTreeGen contains the functions:
tree_from_alignment(input_alignment_path, tree_format) which will take an alignment and generate a Neighbor Joining Tree with Branch Support;
gen_boostrap_consensus_tree(input_alignment_path, replicate_count) which will generate a consensus Bootstrap Tree based on a set number of replicates with Branch Support;
run_iqtree(input_alignment_path, save_path) which will run the IQ-Tree generator and Support the final Tree's Branches


### OUTPUT FILES

The outputs from this will either go into the provided Fasta Folder or a new one based on search criteria will be generated.
The following files will be present:

-ALL Relevent FASTAS Files (if downloaded)

-Combined FASTA File

-MUSCLE Alignment FASTA of the Combined Sequences

-A Neighbor Joining Tree File with Branch Support

-A Bootstrap Consensus Tree File with Branch Support

-An IQ-Tree Series of Files, indlucing a tree with Branch Support


### UPCOMING UPDATES

CSV FILE OUTPUT (Harte S.)
LOCATION EXTRACTION AND TRACKING (Harte S.)
