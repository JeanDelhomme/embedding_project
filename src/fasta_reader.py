# -*- coding: utf-8 -*-
"""
Created on Thu Sep  8 12:45:51 2022

@author: jeand
"""

def fasta_name(file):
    """Captures the protein name in a fasta file.

    Parameters
    ----------
    file : str 
        The name of a fasta file.

    Returns
    -------
    list
        A list containing the sequence.
    """

    # name will be the name of the protein. We are initiating it as an empty string.
    name = ""

    # Opens the FASTA file.
    with open(file, "r") as fasta:
        
        # Reads each line of the FASTA file.
        for line in fasta:
            # Strips the line of "/n" and adds it to sequence if the line
            # is not beggining with ">".
            if line.startswith(">"):
                name += line.strip()[1:]

    # The function returns sequence.
    return name




def fasta_reader(file):
    """Generates a FASTA sequence in a list from a fasta file.

    Parameters
    ----------
    file : str 
        The name of a fasta file.

    Returns
    -------
    list
        A list containing the sequence.
    """
    
    # sequence will be the fasta sequence in form of a list. We initiating it
    # as an empty string.
    sequence = ""
    
    # Opens the FASTA file.
    with open(file, "r") as fasta:
        
        # Reads each line of the FASTA file.
        for line in fasta:
            # Strips the line of "/n" and adds it to sequence if the line
            # is not beggining with ">".
            if not line.startswith(">"):
                sequence += line.strip()
                
    # The function returns sequence.
    return list(sequence)