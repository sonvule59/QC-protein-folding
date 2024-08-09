import time
from dwave.system import DWaveSampler, EmbeddingComposite

# Function to extract sequence from the PDB file
def extract_sequence_from_pdb(pdb_file_path):
    sequence = ""
    with open(pdb_file_path, 'r') as pdb_file:
        residue_dict = {}
        for line in pdb_file:
            if line.startswith("ATOM"):
                residue_number = int(line[22:26].strip())
                residue_name = line[17:20].strip()
                if residue_number not in residue_dict:
                    residue_dict[residue_number] = residue_name
        # Convert to sequence
        for residue_number in sorted(residue_dict.keys()):
            sequence += three_to_one(residue_dict[residue_number])
    return sequence

# Function to convert 3-letter residue codes to 1-letter codes
def three_to_one(residue_name):
    conversion_dict = {
        "ALA": "A", "CYS": "C", "ASP": "D", "GLU": "E", "PHE": "F", "GLY": "G",
        "HIS": "H", "ILE": "I", "LYS": "K", "LEU": "L", "MET": "M", "ASN": "N",
        "PRO": "P", "GLN": "Q", "ARG": "R", "SER": "S", "THR": "T", "VAL": "V",
        "TRP": "W", "TYR": "Y"
    }
    return conversion_dict.get(residue_name, "X")  # "X" for unknown residues

# Load sequence from the PDB file
pdb_file_path = "EF4_Mex.pdb"
sequence = extract_sequence_from_pdb(pdb_file_path)
N = len(sequence)  # Length of the EF-hand 4 region
num_aa = 19  # Considering 19 possible amino acids

print(f"Extracted Sequence: {sequence}")

# List of amino acids
aa_list = ['G', 'I', 'F', 'A', 'P', 'L', 'V', 'W', 'Y', 'M', 'S', 'T', 'D', 'E', 'C', 'N', 'Q', 'R', 'H', 'K']

# Quantum Computing Approach
def create_qubo_matrix(N, num_aa):
    size = N * num_aa
    Q = [[0 for _ in range(size)] for _ in range(size)]

    for i in range(size):
        Q[i][i] = -1

    for i in range(size):
        for j in range(i + 1, size):
            if i // num_aa == j // num_aa:
                Q[i][j] = 2

    return Q

# Create QUBO matrix
Q = create_qubo_matrix(N, num_aa)

# Measure runtime for Quantum Annealing
start_time = time.time()

# Set up the D-Wave sampler
sampler = EmbeddingComposite(DWaveSampler())

# Perform quantum annealing to explore variants
response = sampler.sample_qubo(Q, num_reads=100)

quantum_time = time.time() - start_time
print(f"Quantum Annealing Time: {quantum_time} seconds")

# Extract explored variants
explored_variants_quantum = set()
for sample in response.samples():
    variant = []
    for i in range(N):
        aa_index = max(range(num_aa), key=lambda x: sample[i * num_aa + x])
        variant.append(aa_list[aa_index])
    explored_variants_quantum.add(tuple(variant))

print(f"Number of unique variants explored by Quantum Annealing: {len(explored_variants_quantum)}")

# Brute-force Approach
# Measure runtime for Brute-force exploration
start_time = time.time()

explored_variants_bruteforce = set()

for aa1 in aa_list:
    for aa2 in aa_list:
        for aa3 in aa_list:
            # Create a sequence variant for the first 3 positions (brute-force for N=3)
            variant = (aa1, aa2, aa3)
            explored_variants_bruteforce.add(variant)

bruteforce_time = time.time() - start_time
print(f"Brute-force Time: {bruteforce_time} seconds")
print(f"Number of unique variants explored by Brute-force: {len(explored_variants_bruteforce)}")

# Write results to an output file
output_filename = "mutation_comparison_output.txt"
with open(output_filename, "w") as f:
    f.write("Quantum Annealing Results:\n")
    f.write(f"Quantum Annealing Time: {quantum_time} seconds\n")
    f.write(f"Number of unique variants explored by Quantum Annealing: {len(explored_variants_quantum)}\n")
    f.write("\nBrute-force Results:\n")
    f.write(f"Brute-force Time: {bruteforce_time} seconds\n")
    f.write(f"Number of unique variants explored by Brute-force: {len(explored_variants_bruteforce)}\n")

print(f"Results written to {output_filename}")

# Comparison of Results (also written to file)
with open(output_filename, "a") as f:
    f.write("\nComparison of Results:\n")
    f.write(f"Quantum Annealing Time: {quantum_time} seconds\n")
    f.write(f"Brute-force Time: {bruteforce_time} seconds\n")
    f.write(f"Number of variants by Quantum: {len(explored_variants_quantum)}\n")
    f.write(f"Number of variants by Brute-force: {len(explored_variants_bruteforce)}\n")
