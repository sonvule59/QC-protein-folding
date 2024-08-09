import time
import numpy as np
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
N = 3

# Load the PDB file
pdb_filename = "LanM_clean.pdb"
pose = pose_from_pdb(pdb_filename)

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
    pose_to_mutate = Pose()
    pose_to_mutate.assign(ppose)

    for i in range(N):
        aa_index = np.argmax(sample[i*num_aa:(i+1)*num_aa])
        mutate_residue(pose_to_mutate, ef4_positions[i], aa_list[aa_index])

    explored_variants_quantum += 1
    print(f"Explored variant {explored_variants_quantum}: {pose_to_mutate.sequence()}")

print(f"Total variants explored by quantum annealing: {explored_variants_quantum}")
print(f"Quantum annealing time: {quantum_time} seconds")

# Brute-force approach
start_time = time.time()
explored_variants_bruteforce = 0

# Nested loop for brute-force exploration
for aa1 in aa_list:
    for aa2 in aa_list:
        for aa3 in aa_list:
            pose_to_mutate = Pose()
            pose_to_mutate.assign(pose)

            mutate_residue(pose_to_mutate, ef4_positions[0], aa1)
            mutate_residue(pose_to_mutate, ef4_positions[1], aa2)
            mutate_residue(pose_to_mutate, ef4_positions[2], aa3)

            explored_variants_bruteforce += 1

print(f"Total variants explored by brute-force: {explored_variants_bruteforce}")
print(f"Brute-force time: {time.time() - start_time} seconds")
