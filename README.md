# Phylogenetic_Pipeline_

REQUIRES MUSCLE FOR ALIGNMENT
1) Download it at: https://2018-03-06-ibioic.readthedocs.io/en/latest/install_muscle.html
2) Replace the 'muscle5.1.win64.exe' file with the downloaded version's executable

This pipeline is designed with two (2) different input pathways:
1) A folder already containing Fasta files
OR
2) An NCBI Search criteria, currently formatted as "{Gene name}[All Fields] AND {Species name}[Organism]"

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
