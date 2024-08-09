import itertools
import time
import random
import numpy as np

# Define the number of positions and possible amino acids
N = 2  # Number of positions 
amino_acids = 'ACDEFGHIKLMNPQRSTVWY'  # 20 standard amino acids, excluding one (e.g., 'C')


"""
import os
import pyrosetta
from pyrosetta import pose_from_pdb, get_fa_scorefxn
from pyrosetta.rosetta.protocols.relax import FastRelax
from pyrosetta.toolbox.mutants import mutate_residue, restrict_non_nbrs_from_repacking

from pyrosetta import *
init()
aa_list = ['G','I','F','A','P','L','V','W','Y','M','S','T','D','E','C','N','Q','R','H','K']

pose = pose_from_pdb("LA.pdb")

## This nested for loop will go though all the variants
for i in range(100):
    pose_idx = i
    for aa in aa_list:
        pose_to_mutate = Pose()
        pose_to_mutate.assign(pose)
        mutate_residue(pose_to_mutate, pose_number, aa, sfxn)
"""
# Define the fitness evaluation function (simplified example)
def evaluate_fitness(variant):
    # Compute fitness based on the variant configuration
    # Placeholder: just sum the indices of the amino acids
    fitness = sum(amino_acids.index(aa) for aa in variant)
    return fitness

# Function to compute the Metropolis acceptance probability
def metropolis_acceptance(delta_fitness, temperature):
    if delta_fitness < 0:
        return True
    else:
        return np.exp(-delta_fitness / temperature) > random.random()

# Start timing
start_time = time.time()

# Initial variant and its fitness
current_variant = tuple(random.choice(amino_acids) for _ in range(N))
current_fitness = evaluate_fitness(current_variant)

# Store the accepted variants
accepted_variants = [(current_variant, current_fitness)]

# Define the temperature for the Metropolis criterion
temperature = 1.0

# Generate and evaluate variants
total_variants = 0
max_variants = 1000000  # Limit the number of variants to explore for feasibility

while total_variants < max_variants:
    new_variant = tuple(random.choice(amino_acids) for _ in range(N))
    new_fitness = evaluate_fitness(new_variant)
    delta_fitness = new_fitness - current_fitness

    if metropolis_acceptance(delta_fitness, temperature):
        current_variant = new_variant
        current_fitness = new_fitness
        accepted_variants.append((current_variant, current_fitness))

    total_variants += 1

# End timing
sampling_time = time.time() - start_time

# Print results
print(f"Sampling time: {sampling_time} seconds")
print(f"Number of variants sampled: {total_variants}")
print(f"Number of variants accepted: {len(accepted_variants)}")

# Print a few accepted variants and their fitness
for i, (variant, fitness) in enumerate(accepted_variants):
    print(f"Variant {i + 1}: {variant}, Fitness: {fitness}")
