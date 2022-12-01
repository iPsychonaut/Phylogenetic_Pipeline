# Phylogenetic_Pipeline_

REQUIRES MUSCLE FOR ALIGNMENT
1) Download it at: https://2018-03-06-ibioic.readthedocs.io/en/latest/install_muscle.html
2) Replace the 'muscle5.1.win64.exe' file with the downloaded version's executable

This pipeline is designed with two (2) different input pathways:
1) A folder already containing Fasta files
OR
2) An NCBI Search criteria, currently formatted as "{Gene name}[All Fields] AND {Species name}[Organism]"

What Each File Does:

-BioPyFolderSeqs contains the function combine_fastas(input_directory) which will combine all Fasta files in a given directory into a single fasta for alignment

-BioPyNCBISeqs contains the function search_ncbi(user_email, search_term, return_number, search_db, save_path) which will search NCBI for, download, and combine a desired number of Fasta files for alignment

-BioPySeqMUSCLE contains the function MUSCLE_alignment(combined_path, save_path) which will take a Combined Fasta and output a MUSCLE alignment

-BioPyTreeFigures contains the function trees_from_alignment(input_alignment, save_path) which will take an alignment and generate a Phyloxml Tree, Nexus Tree, and IQ-Tree2 series

The output of this will either go into the provided Fasta Folder or a new one based on search criteria will be generated.
The following files will be present:

-ALL Relevent FASTAS Files (if downloaded)

-Combined FASTA File

-MUSCLE Alignment of the Combined FASTA

-A Phyloxml Tree File

-A Nexus Tree File

-An SVG of the Nexus Tree

-An IQ-Tree2 Series of Files

UPCOMING UPDATES:
Simple GUI for ease of user input, general data viewing, and aesthetics
