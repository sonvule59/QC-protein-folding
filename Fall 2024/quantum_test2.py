# File: protein_mutation_comparison.py

import time
import random
from dwave.system import DWaveSampler, EmbeddingComposite
import pyrosetta
from pyrosetta import pose_from_pdb, Pose
from pyrosetta.toolbox.mutants import mutate_residue

# Initialize PyRosetta
pyrosetta.init()

# List of amino acids
aa_list = ['G', 'I', 'F', 'A', 'P', 'L', 'V', 'W', 'Y', 'M', 'S', 'T', 'D', 'E', 'C', 'N', 'Q', 'R', 'H', 'K']
num_aa = len(aa_list)

# Length of the protein's EF-hand 4 and its flanking regions
N = 12

# Load the PDB file
pdb_filename = "LanM_clean.pdb"
pose = pose_from_pdb(pdb_filename)

# Function to create a QUBO matrix
def create_qubo_matrix(N, num_aa):
    size = N * num_aa
    Q = [[0 for _ in range(size)] for _ in range(size)]

    for i in range(size):
        Q[i][i] = -1

    for i in range(size):
        for j in range(i + 1, size):
            if i // num_aa == j // num_aa:
                Q[i][j] = Q[j][i] = 2

    return Q

# QUBO matrix for the problem
Q = create_qubo_matrix(N, num_aa)

# Define the QUBO problem
sampler = EmbeddingComposite(DWaveSampler())

# Quantum annealing: sample solutions
start_time = time.time()
response = sampler.sample_qubo(Q, num_reads=100)
quantum_time = time.time() - start_time

# Extract and count quantum annealing results
explored_variants_quantum = 0
quantum_variants = []

ef4_positions = list(range(108, 120))

for sample in response.samples():
    pose_to_mutate = Pose()
    pose_to_mutate.assign(pose)
    for i in range(N):
        sample_values = [sample[j] for j in range(i * num_aa, (i + 1) * num_aa)]
        aa_index = sample_values.index(1)
        mutate_residue(pose_to_mutate, ef4_positions[i], aa_list[aa_index])

    explored_variants_quantum += 1
    quantum_variants.append(pose_to_mutate.sequence())

# Brute-force approach
start_time = time.time()
explored_variants_bruteforce = 0
bruteforce_variants = []

# Nested loop for brute-force exploration
for aa1 in aa_list:
    for aa2 in aa_list:
        for aa3 in aa_list:
            for aa4 in aa_list:
                for aa5 in aa_list:
                    for aa6 in aa_list:
                        for aa7 in aa_list:
                            for aa8 in aa_list:
                                for aa9 in aa_list:
                                    for aa10 in aa_list:
                                        for aa11 in aa_list:
                                            for aa12 in aa_list:
                                                pose_to_mutate = Pose()
                                                pose_to_mutate.assign(pose)

                                                mutate_residue(pose_to_mutate, ef4_positions[0], aa1)
                                                mutate_residue(pose_to_mutate, ef4_positions[1], aa2)
                                                mutate_residue(pose_to_mutate, ef4_positions[2], aa3)
                                                mutate_residue(pose_to_mutate, ef4_positions[3], aa4)
                                                mutate_residue(pose_to_mutate, ef4_positions[4], aa5)
                                                mutate_residue(pose_to_mutate, ef4_positions[5], aa6)
                                                mutate_residue(pose_to_mutate, ef4_positions[6], aa7)
                                                mutate_residue(pose_to_mutate, ef4_positions[7], aa8)
                                                mutate_residue(pose_to_mutate, ef4_positions[8], aa9)
                                                mutate_residue(pose_to_mutate, ef4_positions[9], aa10)
                                                mutate_residue(pose_to_mutate, ef4_positions[10], aa11)
                                                mutate_residue(pose_to_mutate, ef4_positions[11], aa12)

                                                explored_variants_bruteforce += 1
                                                bruteforce_variants.append(pose_to_mutate.sequence())

bruteforce_time = time.time() - start_time

# Write results to an output file
output_filename = "mutation_comparison_output.txt"
with open(output_filename, "w") as f:
    f.write("Quantum Annealing Results:\n")
    f.write(f"Total variants explored: {explored_variants_quantum}\n")
    f.write(f"Quantum annealing time: {quantum_time} seconds\n")
    for idx, variant in enumerate(quantum_variants):
        f.write(f"Variant {idx+1}: {variant}\n")

    f.write("\nBrute-Force Results:\n")
    f.write(f"Total variants explored: {explored_variants_bruteforce}\n")
    f.write(f"Brute-force time: {bruteforce_time} seconds\n")
    for idx, variant in enumerate(bruteforce_variants):
        f.write(f"Variant {idx+1}: {variant}\n")

print(f"Results written to {output_filename}")
