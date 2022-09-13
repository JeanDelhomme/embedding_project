# -*- coding: utf-8 -*-
"""
Created on Thu Sep  8 12:45:51 2022

@author: Jean Delhomme

This file contains two functions :
    - fasta_name
    - fasta_reader
    
fasta_name takes the name of a protein from a fasta file.
fasta_reader takes the sequence of a protein from a fasta file.

"""

###############################################################################
#                                                                             #
#                                fasta_name                                   #
#                                                                             #
###############################################################################

def fasta_name(file):
    """Captures the protein name in a fasta file.

    Parameters
    ----------
    file : str 
        The name of a fasta file.

    Returns
    -------
    str
        The name of the protein.
    """

    # name will be the name of the protein. 
    # We are initiating it as an empty string.
    name = ""

    # Opens the FASTA file.
    with open(file, "r") as fasta:
        
        # Reads each line of the FASTA file.
        for line in fasta:
            # Strips the line of "/n" and adds it to name if the line
            # beggins with ">". 
            # ">" is removed.
            if line.startswith(">"):
                name += line.strip()[1:]

    return name

###############################################################################
#                                                                             #
#                              fasta_sequence                                 #
#                                                                             #
###############################################################################

def fasta_reader(file):
    """Generates a FASTA sequence in a list from a fasta file.

    Parameters
    ----------
    file : str 
        The name of a fasta file.

    Returns
    -------
    list
        A list of strings containing the sequence.
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
                
    return list(sequence)