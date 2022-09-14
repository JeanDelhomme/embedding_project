# Embedding_project

## Download this repository

```bash
git clone https://github.com/JeanDelhomme/embedding_project
cd embedding_project/
```

## Install conda environment

Install [conda](https://docs.conda.io/en/latest/miniconda.html).

Install mamba:

```bash
conda install mamba -n base -c conda-forge
```

Create conda environment and install dependendies:

```bash
mamba env create -f binder/environment.yml
```

Load conda environment:

```bash
conda activate 3d-genome-builder
```

## Add data

### fasta files

Two independent fasta files containing the sequence of each protein must be located in `embeding_project/data/fasta/`.

### embedding files

Two independent embedding files containing the encoded sequence of each protein as embeddings must be located in `embeding_project/data/emb/`.

> **Note**
>
> Matching fasta file and embedding files must be from the same proteic sequence.
> Fasta files must use .fasta extension.
> Embedding files must use .t5emb extension and are produced following the [T5 Prot Trans method](https://github.com/agemagician/ProtTrans).

## Run the alignment

Go to the src directory :

```bash
cd embedding_project/src/
```

Run main.py :

```bash
python3 main.py argv[1] argv[2] argv[3] argv[4] argv[5]
```

Parameters :

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

Returns :

file.txt

    A text file containing the alignment result. The file is created in the
    results repository.

## Exemples :

```bash
python3 main.py adk_2ak3a.t5emb ADK_2AK3A.fasta 6PF2K_1bif.t5emb 6PF2K_1BIF.fasta nw
```

Your can try with other sequences and alignment methods :

```bash
python3 main.py 5_3_exonuclease_1bgxt.t5emb 5_3_EXONUCLEASE_1BGXT.fasta 6PF2K_1bif.t5emb 6PF2K_1BIF.fasta nw
python3 main.py 5_3_exonuclease_1bgxt.t5emb 5_3_EXONUCLEASE_1BGXT.fasta 6PF2K_1bif.t5emb 6PF2K_1BIF.fasta sw
python3 main.py 5_3_exonuclease_1bgxt.t5emb 5_3_EXONUCLEASE_1BGXT.fasta 6PF2K_1bif.t5emb 6PF2K_1BIF.fasta gl
```