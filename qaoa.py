import numpy as np
from qiskit.circuit import QuantumCircuit, Parameter
from qiskit.visualization import plot_histogram
from qiskit_algorithms import QAOA
from qiskit.primitives import Sampler, Estimator
from qiskit_optimization import QuadraticProgram
from qiskit_optimization.algorithms import MinimumEigenOptimizer
from qiskit_aer import Aer, AerSimulator
from qiskit_algorithms.optimizers import COBYLA

# # Define the problem
# num_amino_acids = 5
# num_mutations = 19

# # Define the cost function for mutation discrimination capability
# def cost_function(sequence):

#     return np.random.random()  # Replace with actual computation

# qp = QuadraticProgram()
# for i in range(num_amino_acids):
#     for j in range(num_mutations):
#         qp.binary_var(f'x_{i}_{j}')

# #  objective function
# objective = np.zeros((num_amino_acids * num_mutations, num_amino_acids * num_mutations))
# for i in range(num_amino_acids):
#     for j in range(num_mutations):
#         objective[i*num_mutations + j, i*num_mutations + j] = cost_function([i, j])

# qp.minimize(quadratic=objective)

# # Define the QAOA algorithm
# simulator = AerSimulator()
# qaoa = QAOA(optimizer=COBYLA(), quantum_instance=simulator)
# optimizer = MinimumEigenOptimizer(qaoa)

# result = optimizer.solve(qp)

# # Extract the solution
# solution = result.x

# # Print the solution
# print("Optimal mutation sequence:")
# for i in range(num_amino_acids):
#     for j in range(num_mutations):
#         if solution[i*num_mutations + j] == 1:
#             print(f"Amino acid {i}: Mutation {j}")

# #To do
# # Evaluate sequences using ProteinMPNN and PyRosetta for their discrimination capability

sequence = "ZZZZZZZZZZASGADALKALNKDNDDSLEIAEVIHAGATTFTAINPDGDTTLESGETKGRLTEKDWARANKDGDQTLEMDEWLKILRTRFKRADANKDGKLTAAELDSKAGQGVLVMIMKASGADALKALNKDNDDSLEIAEVIHAGATTFTAINPDGDTTLESGETKGRLTEKDWARANKDGDQTLEMDEWLKILRTRFKRADANKDGKLTAAELDSKAGQGVLVMIM"

# Define a segment of interest for mutation
segment_start = 10
segment_length = 5
segment = sequence[segment_start:segment_start + segment_length]

# Number of possible mutations per amino acid (for simplicity, assume 3 possible mutations)
num_mutations = 3

qp = QuadraticProgram()
# num_amino_acids = 5
# num_mutations = 2
# Add binary variables for each amino acid mutation
for i in range(segment_length):
    for j in range(num_mutations):
        qp.binary_var(f'x_{i}_{j}')
# for i in range(num_amino_acids):
#     for j in range(num_mutations):
#         qp.binary_var(f'x_{i}_{j}')

# For simplicity, let's consider binary variables representing mutations
# qp.binary_var('x_0')
# qp.binary_var('x_1')

# Define the cost Hamiltonian
# For simplicity, consider a simple cost function: H_C = Z0 + Z1 + Z0*Z1
# This is a toy example and should be replaced with the actual cost function for protein discrimination

# Define the linear and quadratic terms
# linear_terms = {f'x_{i}_{j}': np.random.rand() for i in range(num_amino_acids) for j in range(num_mutations)}
# quadratic_terms = {(f'x_{i}_{j}', f'x_{k}_{l}'): np.random.rand() 
#                    for i in range(num_amino_acids) 
#                    for j in range(num_mutations) 
#                    for k in range(i+1, num_amino_acids) 
#                    for l in range(num_mutations)}


linear_terms = {f'x_{i}_{j}': float(np.random.rand()) for i in range(segment_length) for j in range(num_mutations)}
quadratic_terms = {(f'x_{i}_{j}', f'x_{k}_{l}'): float(np.random.rand()) 
                   for i in range(segment_length) 
                   for j in range(num_mutations) 
                   for k in range(i+1, segment_length) 
                   for l in range(num_mutations)}



# linear_terms = {'x_0': 1, 'x_1': 1}
# quadratic_terms = {('x_0', 'x_1'): 1}

# Set the objective function
qp.minimize(linear=linear_terms, quadratic=quadratic_terms)
for i in range(segment_length):
    qp.linear_constraint([f'x_{i}_{j}' for j in range(num_mutations)], '=', 1, f'c_{i}')

# Print the objective function to debug
print("Objective function terms:")
print("Linear terms:", linear_terms)
print("Quadratic terms:", quadratic_terms)
print("Constraints:")
for i in range(segment_length):
    print(f'Amino acid {i}:', [f'x_{i}_{j}' for j in range(num_mutations)])

# Define the QAOA algorithm
simulator = AerSimulator()
sampler = Sampler()

qaoa = QAOA(optimizer=COBYLA(maxiter=50), reps=1, sampler=sampler)

# Create the Minimum Eigen Optimizer based on QAOA
optimizer = MinimumEigenOptimizer(qaoa)

# Output

# print("Optimal mutation sequence:")
# print(f"x_0: {solution[0]}")
# print(f"x_1: {solution[1]}")

try:
    result = optimizer.solve(qp)
    # Extract the solution
    solution = result.x

    # Print the solution with debug information
    print("Optimal mutation sequence:")
    for i in range(segment_length):
        mutations = [solution[i*num_mutations + j] for j in range(num_mutations)]
        selected_mutations = [j for j in range(num_mutations) if mutations[j] == 1.0]
        print(f"Amino acid {i}: Selected mutations: {selected_mutations}")

    # Additional debug information
    print("Optimization result:", result)
    print("Objective function value:", result.fval)
except Exception as e:
    print("An error occurred:", e)

# try:
#     result = optimizer.solve(qp)
#     # Extract the solution
#     solution = result.x

#     # Print the solution with debug information
#     print("Optimal mutation sequence:")
#     for i in range(num_amino_acids):
#         mutations = [solution[i*num_mutations + j] for j in range(num_mutations)]
#         selected_mutations = [j for j in range(num_mutations) if mutations[j] == 1.0]
#         print(f"Amino acid {i}: Selected mutations: {selected_mutations}")

#     # Additional debug information
#     print("selected-mutation", selected_mutations)
#     print("Optimization result:", result)
#     print("Objective function value:", result.fval)
# except Exception as e:
#     print("An error occurred:", e)
"""
To do
Include more amino acids and possible mutations from PyRosetta.
"""