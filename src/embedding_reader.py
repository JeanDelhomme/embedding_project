# -*- coding: utf-8 -*-
"""
Created on Thu Sep  8 14:44:59 2022

@author: jeand
"""

import numpy as np

def embedding_reader(file):
    """Transforms a .t5emb file into an array of vectors.
    
    In a .t5emb file, each line corresponds to a residue and contains 1024 
    values. These values represent the physico-chemecal propreties and 
    environnement of the residue.
    
    The function creates a list containing a vector of 1024 values for each 
    residue.

    Parameters
    ----------
    file : str 
        The name of an embedding file.

    Returns
    -------
    array
        An array containing a vector for each residue.
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
       
    # The function returns sequence as an array.
    return np.array(sequence)