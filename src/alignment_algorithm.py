# -*- coding: utf-8 -*-
"""
Created on Tue Sep 13 15:00:34 2022

@author: Jean Delhomme

This file contains three functions :
    - needleman_wunsch
    - smith_waterman
    - glocal
    
Those functions are able to find an alignment between two proteic sequences
following the corresponding algorithms.

"""

import numpy as np

###############################################################################
#                                                                             #
#                           Needleman and Wunsch                              #
#                                                                             #
###############################################################################

def needleman_wunsch(fasta1, fasta2, prot_name1, prot_name2, alignment_matrix):
    """Performs a global alignment following the Needleman and Wunsch algorithm.

    Parameters
    ----------
    fasta1 : list of sting
        A list of string containing the first fasta sequence.
        
    fasta2 : list of string
        A list of string containing the second fasta sequence.
        
    prot_name1 : string
        The name of the first protein.
        
    prot_name2 : string
        The name of the second protein.
        
    alignment_matrix : array
        An array containing the alignment matrix.

    Returns
    -------
    str
        A string presenting the alignment.
    """
    
    # result1 and result2 will contain the aligned sequences.
    result1 = ""
    result2 = ""
    
    # i and j are set to the size of the two sequences.
    # i for rows and j for columns.
    i = alignment_matrix.shape[0]-1
    j = alignment_matrix.shape[1]-1
    
    # Gives the score of the alignment.
    alignment_score = alignment_matrix[i][j]
    
    # The while loop starts at the bottom-left corner and finishes at the 
    # top-right corner.
    while i > 0 and j > 0:

        score_diagonal = alignment_matrix[i-1][j-1]
        score_left = alignment_matrix[i][j-1]
        score_top = alignment_matrix[i-1][j]
        max_score = max(score_diagonal, score_left, score_top)
        
        # Investigates to know the optimal path.
        # There is 3 solutions, the best path is from the diagonal, the top or
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
    
    # If i = 0 or j = 0, we need to complete the rest of the sequence with gaps.
    while j > 0:
        result1 += fasta1[i-1]
        result2 += '-'
        j -= 1
    while i > 0:
        result1 += '-'
        result2 += fasta2[j-1]
        i -= 1
    
    # The sequences have been created from the bottom-right to the top-left.
    # We need to reverse the two sequences.
    result1 = result1[::-1]
    result2 = result2[::-1]

    # Creates a text file as output and writes the results in it.
    writting_model = (f"Global alignment of {prot_name1} and {prot_name2}\n\n"
                      f"Alignment_score = {alignment_score}\n\n"
                      f"{prot_name1}\n"
                      f"{result1}\n"
                      f"{result2}\n"
                      f"{prot_name2}\n\n")

    file_name = (f"../results/Global_{prot_name1}_&_{prot_name2}.txt")

    with open (file_name, "w") as file:
        file.write(writting_model)

###############################################################################
#                                                                             #
#                           Smith and Waterman                                #
#                                                                             #
###############################################################################

def smith_waterman(fasta1, fasta2, prot_name1, prot_name2, alignment_matrix):
    """Performs a local alignment following the Smith and Waterman algorithm.

    Parameters
    ----------
    fasta1 : list of sting
        A list of string containing the first fasta sequence.
        
    fasta2 : list of string
        A list of string containing the second fasta sequence.
        
    prot_name1 : string
        The name of the first protein.
            
    prot_name2 : string
        The name of the second protein.
        
    alignment_matrix : array
        An array containing the alignment matrix.

    Returns
    -------
    str
        A string presenting the alignment.
    """
    
    # strating_position corresponds to the location of the maximum score in 
    # the entire table.
    starting_position = np.where(alignment_matrix == np.amax(alignment_matrix))
    # max_number corresponds to the number of maximums in the entire table.
    # In case, there is more than one maximum, different alignments are made
    # available.
    max_number = len(list(starting_position[0]))
    
    # We use an if loop in case we have more than one maximum.
    while max_number >= 1:
        
        # result1 and result2 will contain the aligned sequences.
        result1 = ""
        result2 = ""
        
        # i and j are set to the position of the maximum score of the table.
        # i for rows and j for columns.
        i = list(starting_position[0])[max_number-1]
        j = list(starting_position[1])[max_number-1]
        # alignment-score corresponds to the value of the maximum.
        alignment_score = alignment_matrix[i][j]
    
        # The while loop starts at the bottom-left corner and finishes at the 
        # top-right corner.
        while i > 0 and j > 0:

            score_diagonal = alignment_matrix[i-1][j-1]
            score_left = alignment_matrix[i][j-1]
            score_top = alignment_matrix[i-1][j]
            max_score = max(score_diagonal, score_left, score_top)
        
            # Investigates to know the optimal path.
            # There is 3 solutions, the best path is from the diagonal, the top 
            # or the left.
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

        # Creates a text file as output and writes the results in it.
        writting_model = (
                      f"Local alignment of {prot_name1} and {prot_name2}\n\n"
                      f"Alignment_score = {alignment_score}\n\n"
                      f"{prot_name1}\n"
                      f"{result1}\n"
                      f"{result2}\n"
                      f"{prot_name2}\n\n")
        
        file_name = (f"../results/Local_{prot_name1}_&_{prot_name2}.txt")

        with open (file_name, "a") as file:
            file.write(writting_model)
        
        # Reiterates the while loop.
        max_number -= 1

###############################################################################
#                                                                             #
#                                   Glocal                                    #
#                                                                             #
###############################################################################

def glocal(fasta1, fasta2, prot_name1, prot_name2, alignment_matrix):
    """Performs a glocal alignment.
    
    The alignment starts at the maximum score of the last column and ends
    at the first column.

    Parameters
    ----------
    fasta1 : list of sting
        A list of string containing the first fasta sequence.
        
    fasta2 : list of string
        A list of string containing the second fasta sequence.
        
    prot_name1 : string
        The name of the first protein.
            
    prot_name2 : string
        The name of the second protein.
        
    alignment_matrix : array
        An array containing the alignment matrix.

    Returns
    -------
    str
        A string presenting the alignment.
    """
    # strating_position corresponds to the position of the maximum score in 
    # the last column.
    starting_position = np.where(alignment_matrix[:,-1] \
                                 == np.amax(alignment_matrix[:,-1]))
    # max_number is the number of maximums in the last column.
    max_number = len(list(starting_position[0]))
    
    # We use an if loop in case we have more than one maximum.
    while max_number >= 1:
        
        # result1 and result2 will contain the aligned sequences.
        result1 = ""
        result2 = ""
        
        # i and j are set to the maximum score of the table.
        # i for rows and j for columns.
        i = list(starting_position[0])[max_number-1]
        j = int(alignment_matrix.shape[1]-1)
        
        alignment_score = alignment_matrix[i][j]
        
        # The while loop finishes at the first column.
        while j > 0:

            score_diagonal = alignment_matrix[i-1][j-1]
            score_left = alignment_matrix[i][j-1]
            score_top = alignment_matrix[i-1][j]
            max_score = max(score_diagonal, score_left, score_top)
        
            # Investigates to know the optimal path.
            # There is 3 solutions, the best path is from the diagonal, the top 
            # or the left.
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

        # Creates a text file as output and writes the results in it.
        writting_model = (
                      f"Glocal alignment of {prot_name1} and {prot_name2}\n\n"
                      f"Alignment_score = {alignment_score}\n\n"
                      f"{prot_name1}\n"
                      f"{result1}\n"
                      f"{result2}\n"
                      f"{prot_name2}\n\n")
        
        file_name = (f"../results/Glocal_{prot_name1}_&_{prot_name2}.txt")

        with open (file_name, "a") as file:
            file.write(writting_model)
            
        # Reiterates the while loop.
        max_number -= 1
