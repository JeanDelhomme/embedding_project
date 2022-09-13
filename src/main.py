# -*- coding: utf-8 -*-
"""
Created on Thu Sep  8 16:07:07 2022

@author: jeand
"""
# Utiliser import sys et sys.argv pour gérer les arguments que l'on place dans 
# le main.

# try et except pour les messages d'errur.Voir errur complémentaire. 
# cf ch19 et 21.

# ligne de commande 
# python3 main.py 5_3_exonuclease_1bgxt.t5emb 5_3_EXONUCLEASE_1BGXT.fasta 6PF2K_1bif.t5emb 6PF2K_1BIF.fasta nw
# python3 main.py adk_2ak3a.t5emb ADK_2AK3A.fasta 6PF2K_1bif.t5emb 6PF2K_1BIF.fasta nw

import sys
import numpy as np

import embedding_reader as er
import fasta_reader as fr

import alignment_matrix_nw as am_nw
import needleman_wunsch as nw

import alignment_matrix_sw as am_sw
import smith_waterman as sw

import glocal as gl

if __name__ == "__main__":
    
    path_embedding = "../data/emb/"
    path_fasta = "../data/fasta/"
    
    embedding_file1 = path_embedding + sys.argv[1]
    fasta_file1 = path_fasta + sys.argv[2]
    embedding_file2 = path_embedding + sys.argv[3]
    fasta_file2 = path_fasta + sys.argv[4]
    

    embedding1 = er.embedding_reader(embedding_file1)
    embedding2 = er.embedding_reader(embedding_file2)
    
    fasta1 = fr.fasta_reader(fasta_file1)
    fasta2 = fr.fasta_reader(fasta_file2)
    
    prot_name1 = fr.fasta_name(fasta_file1)
    prot_name2 = fr.fasta_name(fasta_file2)
    
    dot_matrix = np.dot(embedding1, embedding2.T)
    
    # For Needleman and Wunsch (global alignment).
    if len(sys.argv) < 6 or sys.argv[5] == "nw":
        
        # Produces the alignment matrix needed to find the best path.
        alignment_matrix = am_nw.alignment_matrix_nw(dot_matrix) 
        # Finds the best path and thus, the sequence alignment.
        nw.needleman_wunsch(fasta1, fasta2, prot_name1, prot_name2, \
                            alignment_matrix)

    # For Smith and Waterman (local alignment).
    elif sys.argv[5] == "sw":
        
        # Produces the alignment matrix needed to find the best path.
        alignment_matrix = am_sw.alignment_matrix_sw(dot_matrix)  
        # Finds the best path and thus, the sequence alignment.
        sw.smith_waterman(fasta1, fasta2, prot_name1, prot_name2, \
                            alignment_matrix)
    
    # For a glocal alignment.
    elif sys.argv[5] == "gl":
        
        # Produces the alignment matrix needed to find the best path.
        alignment_matrix = am_sw.alignment_matrix_sw(dot_matrix)  
        # Finds the best path and thus, the sequence alignment.
        gl.glocal(fasta1, fasta2, prot_name1, prot_name2, alignment_matrix)
        
    
    """
    # ajouter un output sous forme de fichier.
    # ajouter environnement conda.
    # faire les push github.
    """