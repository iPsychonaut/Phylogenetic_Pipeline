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

# Set Executable paths
iqtree_path  = 'E:/iqtree-1.6.12-Windows/iqtree-1.6.12-Windows/bin/iqtree.exe'
mrbayes_path  = 'C:/Users/theda/.spyder-py3/BioPy/mb.3.2.7-win64.exe'

# Open and initiate the Distance Calculator using the Identity model
distance_calculator = DistanceCalculator('identity')

###############################################################################
# ALIGNMENT INPUT FUNCTIONS
###############################################################################
# Function to Generate Trees from an input alignment
def run_iqtree(input_alignment_path, save_path):
    
    # Generate IQ-Tree Files
    os.system(f'{iqtree_path} -s "{input_alignment_path}" -nt AUTO')
        
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
def gen_boostrap_consensus_tree(input_alignment_path, replicate_count):
    print(f'Processing {input_alignment_path} with {replicate_count}x Replicates')
    print('NOTE: THIS CAN TAKE A WHILE')
    
    # Open the alignment file as a MultipleSeqAlignment object 
    with open(input_alignment_path,'r') as aln:
        working_alignment = AlignIO.read(aln, 'fasta')
    global distance_calculator
    boostrap_constructor = DistanceTreeConstructor(distance_calculator)
    output_tree_path = input_alignment_path.replace('.fasta','_bootstrap.tre')
    bootstrap_consensus_tree = bootstrap_consensus(working_alignment, replicate_count, boostrap_constructor, majority_consensus)
    print(f'Saving: {output_tree_path}')
    Phylo.write(bootstrap_consensus_tree, output_tree_path, 'newick')
    
    # Generate and attach Branch Support to the Consensus Tree
    supported_bootstrap_tree = get_branch_support(output_tree_path, 'newick')
    output_supported_tree_path = output_tree_path.replace('.tre','_supported.tre')
    print(f'Saving: {output_supported_tree_path}')
    Phylo.write(supported_bootstrap_tree, output_supported_tree_path, 'newick')    

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
# Function to Generate Trees from an input alignment
# def run_mrbayes(input_alignment_path, save_path):
#   
#     # Execute Mr Bayes command   
#     os.system(f'{mrbayes_path} temp_alingment_path')
#        
#    # Generate and attach Branch Support to the Tree
#    output_tree_path = f'{input_alignment_path}.treefile'
#    supported_tree = get_branch_support(output_tree_path, 'newick')
#    output_supported_tree_path = f'{output_tree_path}_supported.tre'
#    print(f'Saving: {output_supported_tree_path}')
#    Phylo.write(supported_tree, output_supported_tree_path, 'newick')  

# save_path = 'TEST_IMPORT'
# temp_alingment_path = f'{save_path}/combined_MUSCLE_trimAI_relaxed.phylip'
# # gen_nj_tree(temp_alingment_path)
# # gen_boostrap_consensus_tree(temp_alingment_path, 500)
# run_iqtree(temp_alingment_path, save_path)
# save_path = 'C:/Users/theda/.spyder-py3/BioPy/TEST IMPORT'
# temp_alingment_path = f'{save_path}/combined_MUSCLE_trimAI_relaxed.phylip'
