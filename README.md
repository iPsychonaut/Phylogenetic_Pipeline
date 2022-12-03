# Phylogenetic_Pipeline
#####
DESCRIPTION
#####

Generate a MUSCLE Alignment based on either a Folder with Fasta files or NCBI Search Criteria and then outputs four (4) pertinent trees: Neighbor Joining, UPGMA, Consensus Bootstrap, and the IQ-Tree output.


#####
REQUIREMENTS
#####

MUSCLE for alingment
1) Download it at: https://2018-03-06-ibioic.readthedocs.io/en/latest/install_muscle.html
2) Update line 14 in BioPySeqAlign with the path of the MUSCLE executable

TrimAI for alignment
1) Download the command line version it at: http://trimal.cgenomics.org/getting_started_with_trimal_v1.2
2) Update line 16 in BioPySeqAlign with the path of the TrimAI executable

IQ-Tree for Tree Building
1) Download the command line version it at: http://www.iqtree.org/
2) Update line 16 in BioPyTreeGen with the path of the IQ-Tree executable


#####
FILE DESCRIPTIONS
#####

-BioPyFolderSeqs contains the function combine_fastas(input_directory) which will combine all Fasta files in a given directory into a single fasta for alignment

-BioPyNCBISeqs contains the function search_ncbi(user_email, search_term, return_number, search_db, save_path) which will search NCBI for, download, and combine a desired number of Fasta files for alignment

-BioPySeqAlign contains the functions:
MUSCLE_alignment(combined_path, save_path) which will take a Combined Fasta and output a MUSCLE alignment;
trimAI_alignment(input_alignment) to call the AI trimmer;
fasta_to_relaxed_phylip(input_alignment)/fasta_to_strict_phylip(input_alignment) for final Phylip output

-BioPyTreeGen contains the functions:
tree_from_alignment(input_alignment_path, tree_format) which will take an alignment and generate a Phyloxml Tree, Nexus Tree, and IQ-Tree2 series;
gen_boostrap_consensus_tree(input_alignment_path, replicate_count) which will generate a consensus bootstrap tree based on a set number of replicates;
run_iqtree(input_alignment_path, save_path) which will run the IQ-Tree generator


#####
OUTPUT FILES
#####

The outputs from this will either go into the provided Fasta Folder or a new one based on search criteria will be generated.
The following files will be present:

-ALL Relevent FASTAS Files (if downloaded)

-Combined FASTA File

-MUSCLE Alignment of the Combined FASTA

-A Neighbor Joining Tree File

-A UPGMA Tree File

-A Bootstrap Consensus Tree File

-An IQ-Tree Series of Files


#####
UPCOMING UPDATES
#####

Simple GUI for ease of user input, general data viewing, and aesthetics
