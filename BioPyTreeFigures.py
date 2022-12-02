# -*- coding: utf-8 -*-
'''
Created on Wed Nov 30 11:20:08 2022

@author: theda
'''
# Create Tree Files and Figures

# Import Necessary Libraries
from Bio import AlignIO
from Bio import Phylo
import matplotlib
import matplotlib.pyplot as plt
import os
from Bio import Phylo
import pylab

def get_label(leaf):
    return leaf.name

iqtree_path  = 'E:/iqtree-1.6.12-Windows/iqtree-1.6.12-Windows/bin/iqtree.exe'

def trees_from_alignment(input_alignment, save_path):
    # Generate IqTree2 Files
    os.system(f'{iqtree_path} -s "{input_alignment}" -nt AUTO')
    
    # Open the alignment file as a MultipleSeqAlignment object 
    with open(input_alignment,'r') as aln: 
        alignment = AlignIO.read(aln,'fasta')
    print(type(alignment))
    
    # Open and initiate the Distance Calculator using the Identity model 
    from Bio.Phylo.TreeConstruction import DistanceCalculator 
    calculator = DistanceCalculator('identity')
    
    ## Write the Distance Matrix 
    # distance_matrix = calculator.get_distance(alignment)
    # print(distance_matrix)
    
    # Open and initiate the Tree Constructor 
    from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
    constructor = DistanceTreeConstructor(calculator)
    
    # Build the tree 
    output_tree = constructor.build_tree(alignment)
    output_tree.rooted = True
    #print(output_tree)
    
    # Save the tree to a new file 
    output_xml_path = f'{save_path}/output_tree.xml'
    Phylo.write(output_tree, output_xml_path, 'phyloxml')  
   
    # Convert the tree to a different format (optional)
    output_nex_path = f'{save_path}/output_tree.nex'
    Phylo.convert(output_xml_path, 'phyloxml', output_nex_path, 'nexus')
    
    input_nex = f'{save_path}/output_tree.nex'
    tree = Phylo.read(input_nex, 'nexus')
    tree.ladderize()
    Phylo.draw(tree, label_func=get_label, do_show=False)
    pylab.axis('off')
    pylab.savefig(f'{save_path}/output_tree.svg',
                  format = 'svg', bbox_inches = 'tight', dpi = 300)
    return(output_tree)
    
#def matplotlib_tree_png(output_tree):
    # # Create a basic tree 
    # fig = Phylo.draw(output_tree)
    
    # # Make a better looking tree using the features of matplotlib 
    # fig = plt.figure(figsize=(13, 5), dpi=100) # create figure & set the size 
    # matplotlib.rc('font', size=6)              # fontsize of the leaf and node labels 
    # matplotlib.rc('xtick', labelsize=10)       # fontsize of the tick labels
    # matplotlib.rc('ytick', labelsize=10)       # fontsize of the tick labels
    # #turtle_tree.ladderize()
    # axes = fig.add_subplot(1, 1, 1)
    # Phylo.draw(output_tree, axes=axes)
    # fig.savefig('output_cladogram')