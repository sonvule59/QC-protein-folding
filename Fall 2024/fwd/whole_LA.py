import os
import pyrosetta
from pyrosetta import pose_from_pdb, get_fa_scorefxn
from pyrosetta.rosetta.protocols.relax import FastRelax
from pyrosetta.toolbox.mutants import mutate_residue, restrict_non_nbrs_from_repacking

# Initialize PyRosetta with the ignore_zero_occupancy flag set to false
pyrosetta.init("-ignore_zero_occupancy false")

# Define the path to the directory containing the PDB file
target_directory = "/work/ratul1/son/Pyrosetta/pyros_whole_LA"

# Change the current working directory to the target directory
os.chdir(target_directory)

# Define the relative path to the PDB file
pdb_file = "LA.pdb"

# Define the output file
output_file = "mutation_analysis_results_dimer_3.txt"
"""
### make an emptr pose and then iterate
## Access file VE2 rosetta helper.py
mutate_residue()

"""



# Check if the file exists
if not os.path.exists(pdb_file):
    print(f"Error: The file '{pdb_file}' does not exist or cannot be accessed.")
else:
    try:
        pose = pose_from_pdb(pdb_file)
    except Exception as e:
        print(f"Error loading PDB file: {e}")
        pose = None

    if pose and pose.total_residue() > 0:
        # Print pose information to debug residue indices
        print("Pose information:")
        print(pose)

        with open(output_file, "w") as f:
            f.write("Pose information:\n")
            f.write(str(pose) + "\n")

            # Print the sequence and residue numbers to verify indices
            f.write("Residue details:\n")
            for i in range(1, pose.total_residue() + 1):
                res = pose.residue(i)
                f.write(f"Pose index {i}: {res.name3()}\n")

            # Define a scoring function
            def get_energy(pose):
                """Calculates and returns the total energy of the pose."""
                sfxn = get_fa_scorefxn()
                return sfxn(pose)

            # Define a relaxation function
            def relax_pose(pose):
                """Applies FastRelax to the pose."""
                fast_relax = FastRelax()
                fast_relax.set_scorefxn(get_fa_scorefxn())
                fast_relax.apply(pose)
                return pose

            # # Define the ranges of residues to mutate
            # peptide_ranges = [
            #     (21, 32),
            #     (45, 56),
            #     (70, 81),
            #     (94, 105),
            #     (130, 141),
            #     (154, 165),
            #     (179, 190),
            #     (203, 214),
            # ]

             # Define the ranges of residues to mutate
            peptide_ranges = [
                (94, 105), ## EF4 
                (108, 119)
            ]

        
            # Store original energy
            original_energy = get_energy(pose)

            # Iterate through each specified range of residues
            for start, end in peptide_ranges:
                for res_num in range(start, end + 1):
                    # Make a copy of the original pose
                    test_pose = pose.clone()
                    
                    # Check if the residue index is valid
                    if res_num > pose.total_residue():
                        f.write(f"Residue {res_num} is out of bounds\n")
                        continue
                    
                    # Mutate the residue to ALA using mutate_residue
                    mutate_residue(test_pose, res_num, "A", pack_radius=5.0)
                    
                    # Relax the structure
                    relax_pose(test_pose)
                    
                    # Calculate the new energy
                    mutated_energy = get_energy(test_pose)
                    
                    # Calculate the change in energy
                    energy_change = mutated_energy - original_energy
                    
                    # Print the results
                    f.write(f"Residue {res_num}: Energy Change = {energy_change:.2f} kcal/mol\n")

            f.write("Mutation analysis complete.\n")
