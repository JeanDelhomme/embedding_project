# -*- coding: utf-8 -*-
"""
Created on Tue Sep 13 15:49:06 2022

@author: Jean Delhomme

This file contains two functions :
    - alignment_matrix_nw
    - alignment_matrix_nw
    
Those functions are able to calculate the alignment matrix used in Needleman
and Wunsch and Smith and Waterman algorithms.
The alignment matrix of Needleman and Wunsch is also used for glocal alignment.

"""

import numpy as np

###############################################################################
#                                                                             #
#                        Needleman and Wunsch + Glocal                        #
#                                                                             #
###############################################################################

def alignment_matrix_nw(array):
    """Creates an alignment matrix from a score matrix.

    Parameters
    ----------
    array : array
        An array containing the score matrix for an alignment.

    Returns
    -------
    array
        An array containing the alignment matrix.
    """
    
    # The gap is fixed at 0.
    gap = 0
    
    # seq1_size is the size of the first sequence. It will be used for the rows
    # of the alignment matix.
    seq1_size = array.shape[0]
    # seq2_size is the size of the second sequence. It will be used for the 
    # columns of the alignment matix.
    seq2_size = array.shape[1]
    
    # The alignment matrix is created and stored with 0s.
    alignment_matrix = np.zeros((seq1_size+1, seq2_size+1))
    
    # Fills the first columns.
    for i in range(seq1_size+1):
        alignment_matrix[i][0] = i * gap
    
    # Fills the first rows.
    for j in range(seq2_size+1):
        alignment_matrix[0][j] = j * gap
    
    # Fills out all other values in the score matrix.
    for i in range(1, seq1_size+1):
        for j in range(1, seq2_size+1):
            # Calculates the score by checking the top, left, and diagonal 
            # cells.
            diagonal = alignment_matrix[i-1][j-1] + array[i-1][j-1]
            top = alignment_matrix[i-1][j] + gap
            left = alignment_matrix[i][j-1] + gap
            # Record the maximum score from the three possible scores 
            # calculated above.
            alignment_matrix[i][j] = max(diagonal, top, left)
    
    return alignment_matrix

###############################################################################
#                                                                             #
#                           Smith and Waterman                                #
#                                                                             #
###############################################################################

def alignment_matrix_sw(array):
    """Creates an alignment matrix from a score matrix according to the Smith
    and Waterman algorithm.

    Parameters
    ----------
    array : array
        An array containing the score matrix for an alignment.

    Returns
    -------
    array
        An array containing the alignment matrix.
    """
    
    # The gap is fixed at 0.
    gap = 0
    
    # seq1_size is the size of the first sequence. It will be used for the rows
    # of the alignment matix.
    seq1_size = array.shape[0]
    # seq2_size is the size of the second sequence. It will be used for the 
    # columns of the alignment matix.
    seq2_size = array.shape[1]
    
    # The alignment matrix is created and stored with 0s.
    alignment_matrix = np.zeros((seq1_size+1, seq2_size+1))
    
    # Fills the first columns.
    for i in range(seq1_size+1):
        alignment_matrix[i][0] = i * gap
    
    # Fills the first rows.
    for j in range(seq2_size+1):
        alignment_matrix[0][j] = j * gap
    
    # Fills out all other values in the score matrix.
    for i in range(1, seq1_size+1):
        for j in range(1, seq2_size+1):
            # Calculates the score by checking the top, left, 
            # and diagonal cells.
            diagonal = alignment_matrix[i-1][j-1] + array[i-1][j-1]
            top = alignment_matrix[i-1][j] + gap
            left = alignment_matrix[i][j-1] + gap
            # Record the maximum score from the three possible scores 
            # calculated above. 
            # In Smith and Waterman, the maximum needs to be > 0.
            alignment_matrix[i][j] = max(diagonal, top, left, 0)
    
    return alignment_matrix