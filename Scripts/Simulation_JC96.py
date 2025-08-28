#!/usr/bin/env python3

## This script is used to:
## Run 100 simulations of mutational events across the whole genome using the Jukes–Cantor 1996 (JC96) model.


import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import poisson
from Bio import SeqIO

fasta_file = "../../data/tb.ancestor.fasta"  # FASTA file path
# Read the first record in the FASTA file as the reference genome
record = next(SeqIO.parse(fasta_file, "fasta"))
ref_seq = str(record.seq).upper()   # Convert to uppercase
reference_genome = np.array(list(ref_seq))
num_positions = len(reference_genome)  
num_mutations = 2345799    
num_simulations = 100  

lambda_poisson = num_mutations / num_positions  # Expected number of mutations per site

# -------------------------------
# Define bases and JC69 substitution matrix
# -------------------------------
# Assume the reference genome contains only four bases: A, C, G, T
bases = np.array(["A", "C", "G", "T"])
# Construct substitution matrix: for each reference base, possible substitutions (JC69 model, equal probability, cannot mutate to itself)
alt = np.array([
    [1, 2, 3],  # If reference = A (code 0), substitute with C, G, T
    [0, 2, 3],  # If reference = C (code 1), substitute with A, G, T
    [0, 1, 3],  # If reference = G (code 2), substitute with A, C, T
    [0, 1, 2]   # If reference = T (code 3), substitute with A, C, G
])

# Convert reference genome to numeric encoding (A=0, C=1, G=2, T=3)
ref_numeric = np.searchsorted(bases, reference_genome)

# -------------------------------
# Simulate mutation counts at each site per replicate (Poisson distribution)
"simulation_JC69.py" [readonly] 82L, 3730B                           1,1           Top