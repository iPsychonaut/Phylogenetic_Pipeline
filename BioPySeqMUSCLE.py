# -*- coding: utf-8 -*-
"""
Created on Wed Nov 30 12:10:04 2022

@author: theda
"""
# Convert a Compiled Fasta into MUSCLE Alignment

# Import Necessary Libraries
import os

def MUSCLE_alignment(combined_path, save_path):   
    # Build out MUSCLE Command Line
    muscle_exe = r"C:/Users/theda/.spyder-py3/BioPy/muscle5.1.win64.exe" #specify the location of your muscle exe file
    output_alignment = f"{save_path}/combined_MUSCLE.fasta"
    muscle_command = f'{muscle_exe} -align "{combined_path}" -output "{output_alignment}"'
    print(muscle_command)
    os.system(muscle_command)
    
    # Display Output Alignment
    # MUSCLE_alignment = AlignIO.read(output_alignment, "fasta") 
    # print(MUSCLE_alignment)
    return(output_alignment)