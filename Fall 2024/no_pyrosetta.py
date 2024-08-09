import time
import numpy as np
from dwave.system import DWaveSampler, EmbeddingComposite

# List of amino acids
aa_list = ['G', 'I', 'F', 'A', 'P', 'L', 'V', 'W', 'Y', 'M', 'S', 'T', 'D', 'E', 'C', 'N', 'Q', 'R', 'H', 'K']
num_aa = len(aa_list)

# Length of the protein's EF-hand 4 and its flanking regions
N = 3

# Initial protein sequence (placeholder)
protein_sequence = "MEXLANM"  # Example sequence; in practice, use the actual sequence

# Positions to mutate (EF-hand 4 region and flanks)
ef4_positions = [3, 4, 5]  # Example positions; adjust based on actual EF-hand 4 positions

# Function to create a QUBO matrix
def create_qubo_matrix(N, num_aa):
    size = N * num_aa
    Q = np.random.rand(size, size) * 2 - 1  # Random QUBO matrix for demonstration
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
variants_quantum = []

for sample in response.samples():
    mutated_sequence = list(protein_sequence)

    for i in range(N):
        aa_index = np.argmax(sample[i*num_aa:(i+1)*num_aa])
        mutated_sequence[ef4_positions[i]] = aa_list[aa_index]

    explored_variants_quantum += 1
    variant = ''.join(mutated_sequence)
    variants_quantum.append(variant)
    print(f"Explored variant {explored_variants_quantum}: {variant}")

print(f"Total variants explored by quantum annealing: {explored_variants_quantum}")
print(f"Quantum annealing time: {quantum_time} seconds")

# Brute-force approach
start_time = time.time()
explored_variants_bruteforce = 0
variants_bruteforce = []

# Nested loop for brute-force exploration
for aa1 in aa_list:
    for aa2 in aa_list:
        for aa3 in aa_list:
            mutated_sequence = list(protein_sequence)

            mutated_sequence[ef4_positions[0]] = aa1
            mutated_sequence[ef4_positions[1]] = aa2
            mutated_sequence[ef4_positions[2]] = aa3

            explored_variants_bruteforce += 1
            variant = ''.join(mutated_sequence)
            variants_bruteforce.append(variant)

print(f"Total variants explored by brute-force: {explored_variants_bruteforce}")
print(f"Brute-force time: {time.time() - start_time} seconds")

# Results comparison
print("\nVariants explored by Quantum Annealing:")
print(variants_quantum)

print("\nVariants explored by Brute-Force:")
print(variants_bruteforce)
