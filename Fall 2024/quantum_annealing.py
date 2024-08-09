import time
import numpy as np
from dwave.system import DWaveSampler, EmbeddingComposite
import dimod
from pyrosetta import *
init()
from pyrosetta import pose_from_pdb, get_fa_scorefxn
from pyrosetta.rosetta.protocols.relax import FastRelax
from pyrosetta.toolbox.mutants import mutate_residue, restrict_non_nbrs_from_repacking
# Define the number of positions and possible amino acids
N = 5  # Number of positions 
# num_amino_acids = 19
aa_list = ['G','I','F','A','P','L','V','W','Y','M','S','T','D','E','C','N','Q','R','H','K']

pose = pose_from_pdb("LA.pdb")
# Create a QUBO dictionary
Q = {}

# Define the one-hot constraints and interactions
for pos in range(N):
    for aa in range(num_amino_acids):
        Q[(f'pos{pos}_aa{aa}', f'pos{pos}_aa{aa}')] = -1  # Quadratic term ensuring one-hot encoding
        for aa2 in range(aa + 1, num_amino_acids):
            Q[(f'pos{pos}_aa{aa}', f'pos{pos}_aa{aa2}')] = 2  # Interaction term to penalize multiple amino acids

for pos1 in range(N):
    for pos2 in range(pos1 + 1, N):
        for aa1 in range(num_amino_acids):
            for aa2 in range(num_amino_acids):
                interaction_value = np.random.random()  # Placeholder for actual fitness interaction value - Wait on Ratul's team
                Q[(f'pos{pos1}_aa{aa1}', f'pos{pos2}_aa{aa2}')] = interaction_value

bqm = dimod.BinaryQuadraticModel.from_qubo(Q)
# D-Wave sampler
sampler = EmbeddingComposite(DWaveSampler())

start_time = time.time()

num_reads = num_amino_acids**N  # Number of samples 19^N
sampleset = sampler.sample(bqm, num_reads=num_reads)
dwave_sampling_time = time.time() - start_time

sampled_variants = []
for sample, energy in sampleset.data(['sample', 'energy']):
    sampled_variants.append((sample, energy))

# End timing
sampling_time = time.time() - start_time

# Print results
print(f"Sampling time: {sampling_time} seconds")
print(f"Number of variants sampled: {len(sampled_variants)}")

# Print a few sampled variants and their fitness (energy)
for i, (variant, energy) in enumerate(sampled_variants):
    print(f"Variant {i + 1}: {variant}, Energy: {energy}")
