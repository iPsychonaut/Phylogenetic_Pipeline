# -*- coding: utf-8 -*-
"""
Created on Fri Dec  2 08:42:29 2022

Generate Phylogenetics Trees from Alignments

@author: ian.michael.bollinger@gmail.com
"""
###############################################################################
# IMPORT LIBRARIES
###############################################################################
from Bio import Phylo, AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.Phylo.Consensus import majority_consensus, bootstrap_consensus, get_support
import os

# Set IQ-Tree2 Executable path
iqtree_path  = 'E:/iqtree-1.6.12-Windows/iqtree-1.6.12-Windows/bin/iqtree.exe'

# Open and initiate the Distance Calculator using the Identity model
distance_calculator = DistanceCalculator('identity')

###############################################################################
# ALIGNMENT INPUT FUNCTIONS
###############################################################################
# Function to Generate Trees from an input alignment
def run_iqtree(input_alignment_path, save_path):
    # Generate IqTree2 Files
    os.system(f'{iqtree_path} -s "{input_alignment_path}" -nt AUTO')

# Function to generate a Neighbor Joining or UPGMA tree from an input alignment
def tree_from_alignment(input_alignment_path, tree_format):
    global distance_calculator
    # Open the alignment file as a MultipleSeqAlignment object
    input_alignment_extension = input_alignment_path.split('.')[1]
    input_alignment_format = input_alignment_extension
    with open(input_alignment_path,'r') as aln:
        working_alignment = AlignIO.read(aln, input_alignment_format)
    distance_matrix = distance_calculator.get_distance(working_alignment)
    # Open and initiate the appropriate Tree Constructor
    constructor = DistanceTreeConstructor()
    if tree_format == 'upgma' or tree_format == 'UPGMA':        
        output_tree = constructor.upgma(distance_matrix)
    elif tree_format == 'nj' or tree_format == 'Neighbor Joining' or tree_format == 'neighbor joining':        
        output_tree = constructor.nj(distance_matrix)
    else:
        output_tree = constructor.nj(distance_matrix)
    output_tree.rooted = True
    output_tree_path = input_alignment_path.replace(input_alignment_extension,'')
    output_tree_path = f'{output_tree_path}_{tree_format}.tre'
    print(f'Saving: {output_tree_path}')
    Phylo.write(output_tree, output_tree_path, 'newick')   

# Generate Bootstrap Trees
def gen_boostrap_consensus_tree(input_alignment_path, replicate_count):
    print(f'Processing {input_alignment_path} with {replicate_count}x Replicates')
    print('NOTE: THIS CAN TAKE A WHILE')
    # Open the alignment file as a MultipleSeqAlignment object
    input_alignment_extension = input_alignment_path.split('.')[1]
    input_alignment_format = input_alignment_extension
    with open(input_alignment_path,'r') as aln:
        working_alignment = AlignIO.read(aln, input_alignment_format)
    global distance_calculator
    boostrap_constructor = DistanceTreeConstructor(distance_calculator)
    output_alignment_path = input_alignment_path.strip(input_alignment_extension)
    bootstrap_consensus_tree = bootstrap_consensus(working_alignment, replicate_count, boostrap_constructor, majority_consensus)
    output_tree_path = f'{output_alignment_path}_boostrap.tre.tre'
    print(f'Saving: {output_tree_path}')
    Phylo.write(bootstrap_consensus_tree, output_tree_path, 'newick')    

###############################################################################
# TREE INPUT FUNCITONS
###############################################################################

# Function to Get Branch Support For Specific Tree
def get_branch_support(input_tree_path, input_tree_format):
    trees = list(Phylo.parse(input_tree_path, input_tree_format))
    target_tree = trees[0]
    support_tree = get_support(target_tree, trees)
    return(support_tree)


# ###############################################################################
# # DEBUG WORKSPACE
# ###############################################################################
# save_path = 'TEST IMPORT'
# temp_alingment_path = f'{save_path}/combined_MUSCLE_trimAI.fasta'
# tree_from_alignment(temp_alingment_path, 'upgma')
# tree_from_alignment(temp_alingment_path, 'nj')
# gen_boostrap_consensus_tree(temp_alingment_path, 500)
