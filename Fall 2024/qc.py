from qiskit_aer import Aer, AerSimulator
from qiskit.circuit import QuantumCircuit, Parameter
from qiskit.primitives import Sampler, Estimator
from qiskit.quantum_info import SparsePauliOp, Pauli

import qiskit_nature
from qiskit_algorithms.minimum_eigensolvers import VQE, SamplingVQE, QAOA, NumPyMinimumEigensolver
from qiskit_algorithms.optimizers import COBYLA
from qiskit.circuit.library import RealAmplitudes
# from qiskit.primitives import QuantumInstance

X = Pauli('X')
Z = Pauli('Z')
I = Pauli('I')
# Define Hamiltonian
H = -1 * (Z ^ Z) + -1j * (X ^ I) + 1 * (I ^ X)

# Define the ansatz (variational form)
ansatz = RealAmplitudes(num_qubits=2, reps=2)

# Define the VQE instance
vqe = VQE(ansatz, optimizer=COBYLA(), quantum_instance=AerSimulator())

# Run VQE to find the minimum energy
result = vqe.compute_minimum_eigenvalue(H)

print(f"Minimum energy: {result.eigenvalue}")
print(f"Optimal parameters: {result.optimal_point}")