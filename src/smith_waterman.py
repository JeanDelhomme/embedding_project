# -*- coding: utf-8 -*-
"""
Created on Mon Sep 12 11:16:37 2022

@author: jeand
"""

import numpy as np

def smith_waterman(fasta1, fasta2, prot_name1, prot_name2, alignment_matrix):
    """Performs a local alignment following the Smith and Waterman algorithm.

    Parameters
    ----------
    fasta1 : list of sting
        A list of string containing the first fasta sequence.
        
    fasta2 : list of string
        A list of string containing the second fasta sequence.
        
    alignment_matrix : array
        An array containing the alignment matrix.

    Returns
    -------
    str
        A string presenting the alignment.
    """
    
    # strating_position corresponds to the maximum score in the entire table.
    starting_position = np.where(alignment_matrix == np.amax(alignment_matrix))
    max_number = len(list(starting_position[0]))
    
    # We use an if loop in case we have more than one maximum.
    while max_number >= 1:
        
        # result1 and result2 will store the aligned sequences.
        result1 = ""
        result2 = ""
        
        # i and j are set to the maximum score of the table.
        # i for rows and j for columns.
        i = list(starting_position[0])[max_number-1]
        j = list(starting_position[1])[max_number-1]
        
        alignment_score = alignment_matrix[i][j]
    
        # The while loop starts at the bottom-left corner and finishes at the 
        # top-right corner.
        while i > 0 and j > 0:

            score_diagonal = alignment_matrix[i-1][j-1]
            score_left = alignment_matrix[i][j-1]
            score_top = alignment_matrix[i-1][j]
            max_score = max(score_diagonal, score_left, score_top)
        
            # Investigates to know the optimal path.
            # There is 3 solution, the best path is from the diagonal, the top or
            # the left.
            if max_score == score_diagonal:
                result1 += fasta1[i-1]
                result2 += fasta2[j-1]
                i -= 1
                j -= 1
            elif max_score == score_top:
                result1 += fasta1[i-1]
                result2 += '-'
                i -= 1
            elif max_score == score_left:
                result1 += '-'
                result2 += fasta2[j-1]
                j -= 1
        
        # The sequences have been created from the bottom-right to the 
        # top-left. We need to reverse the two sequences.
        result1 = result1[::-1]
        result2 = result2[::-1]
        
        with open (f'../results/Local_{prot_name1}_&_{prot_name2}.txt', "a") as file:
            file.write(
f'Local alignment of {prot_name1} and {prot_name2}\n\n\
Alignment_score = {alignment_score}\n\n\
{prot_name1}\n\
{result1}\n\
{result2}\n\
{prot_name2}\n')
        
        # Reiterates the while loop.
        max_number -= 1