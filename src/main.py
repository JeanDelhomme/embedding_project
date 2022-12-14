# -*- coding: utf-8 -*-
"""
Created on Thu Sep  8 16:07:07 2022

@author: Jean Delhomme

Main program for the embedding_project.

This program is performing global, local or glocal protein alignment using
fasta and embedding files.

Users should follow the instructions described in the README.md file attached
to the embedding_project repository.

Parameters
----------
argv[1] : str
    The name of the embeding file of the first protein.
    
argv[2] : str
    The name of the fasta file of the first protein.

argv[3] : str
    The name of the embeding file of the second protein.

argv[4] : str
    The name of the fasta file of the second protein.

argv[5] : str
    OPTIONAL, the name of the alignment algorithm. 
    By default, Needleman and Wunsch.
    
    nw : Needleman and Wunsch (global).
    sw : Smith and Waterman (local).
    gl : Glocal.    

Returns
-------
file.txt
    A text file containing the alignment result. The file is created in the
    result repository.

"""

# Importation of common modules.
import sys
import numpy as np

# Importation of the modules used for reading the data files.
import embedding_reader as er
import fasta_reader as fr

# Importation of the modules used for the alignment.
import alignment_matrix as am
import alignment_algorithm as aa



if __name__ == "__main__":

###############################################################################
#                                                                             #
#                            data preparation                                 #
#                                                                             #
###############################################################################

    # Variables for the path of the data files.
    path_embedding = "../data/emb/"
    path_fasta = "../data/fasta/"
    
    # Variables containing the arguments given when lauching the program.
    embedding_file1 = path_embedding + sys.argv[1]
    fasta_file1 = path_fasta + sys.argv[2]
    embedding_file2 = path_embedding + sys.argv[3]
    fasta_file2 = path_fasta + sys.argv[4]
    
    # Creation of embedding array for each protein.
    embedding1 = er.embedding_reader(embedding_file1)
    embedding2 = er.embedding_reader(embedding_file2)
    
    # Creation of fasta sequences for each protein.
    fasta1 = fr.fasta_reader(fasta_file1)
    fasta2 = fr.fasta_reader(fasta_file2)
    
    # Creation of names for each protein.
    prot_name1 = fr.fasta_name(fasta_file1)
    prot_name2 = fr.fasta_name(fasta_file2)
    
    # Calculation of the dot_product between each embedding at each position
    # and construction of an array with those dot_products.
    dot_matrix = np.dot(embedding1, embedding2.T)
  
###############################################################################
#                                                                             #
#                                alignments                                   #
#                                                                             #
###############################################################################
    
    # For Needleman and Wunsch (global alignment).
    if len(sys.argv) < 6 or sys.argv[5] == "nw":
        
        # Produces the alignment matrix needed to find the best path.
        alignment_matrix = am.alignment_matrix_nw(dot_matrix) 
        # Finds the best path and thus, the sequence alignment.
        aa.needleman_wunsch(fasta1, fasta2, prot_name1, prot_name2, \
                            alignment_matrix)

    # For Smith and Waterman (local alignment).
    elif sys.argv[5] == "sw":
        
        # Produces the alignment matrix needed to find the best path.
        alignment_matrix = am.alignment_matrix_sw(dot_matrix)  
        # Finds the best path and thus, the sequence alignment.
        aa.smith_waterman(fasta1, fasta2, prot_name1, prot_name2, \
                            alignment_matrix)
    
    # For a glocal alignment.
    elif sys.argv[5] == "gl":
        
        # Produces the alignment matrix needed to find the best path.
        alignment_matrix = am.alignment_matrix_sw(dot_matrix)  
        # Finds the best path and thus, the sequence alignment.
        aa.glocal(fasta1, fasta2, prot_name1, prot_name2, alignment_matrix)