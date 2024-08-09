import time
import numpy as np
from dwave.system import DWaveSampler, EmbeddingComposite

# List of amino acids
aa_list = ['G', 'I', 'F', 'A', 'P', 'L', 'V', 'W', 'Y', 'M', 'S', 'T', 'D', 'E', 'C', 'N', 'Q', 'R', 'H', 'K']
num_aa = len(aa_list)

# Length of the protein's EF-hand 4 and its flanking regions
N = 3

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
for sample in response.samples():
    explored_variants_quantum += 1

# Brute-force approach
start_time = time.time()
explored_variants_bruteforce = 0

# Nested loop for brute-force exploration
for i in range(num_aa):
    for j in range(num_aa):
        for k in range(num_aa):
            explored_variants_bruteforce += 1

bruteforce_time = time.time() - start_time

# Results
print(f"Quantum annealing explored variants: {explored_variants_quantum}")
print(f"Quantum annealing time: {quantum_time} seconds")

print(f"Brute-force explored variants: {explored_variants_bruteforce}")
print(f"Brute-force time: {bruteforce_time} seconds")
