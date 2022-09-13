# -*- coding: utf-8 -*-
"""
Created on Thu Sep  8 15:49:13 2022

@author: jeand
"""

def dot_product(matrix1, matrix2):
    """Calculates the dot product between two list of vectors.
    
    Each matrix is a list of embedding vectors. 
    The function calculates the dot product for each vector of matrix1 against 
    each vector of matrix2.
    
    The results are organised in a third matrix.

    Parameters
    ----------
    matrix1 : list of list of floats
        A matrix contining x vectors of floats.
        
    matrix2 : list of list of floats
        A matrix contining x vectors of floats.

    Returns
    -------
    list
        A matrix contains the dot product results.
    """
    
    # Sequence is the final list containing the vector. We initiating it.
    sequence = []
    
    # Opens the file in the function argument.
    with open(file, "r") as embedding:
        
        # Reads the file line after line.
        for line in embedding:
            # Creates a list containing each of the 1024 values separted by a 
            # space in the .t5emb file.
            vector = line.split()
            # Converts the strings into floats.
            vector = [float(x) for x in vector]
            # Add vector at the end of the sequence list, creating a list of 
            # lists.
            sequence.append(vector)
            
            import numppy as np
            np.dot(matrix1, matrix2.T)
            
    # The function returns sequence.
    return sequence