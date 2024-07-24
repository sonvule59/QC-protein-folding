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