# reverse_mapping [WIP]
Reverse mapping scheme for multi-scale molecular dynamics. Converts coarse-grained (CG) systems to the atomistic level.

## Simplified Algorithm
1. Determines position and orientation of lipid head beads
2. Places corresponding pre-built atomistic structure such that the lipid head positions of the atomistic structure corresponds to that of the CG molecule.
3. If the atomistic molecule is in contact with another molecule (atoms are within a cutoff distance), the molecule is rotated about the z-axis until there are no contacts. If upon rotating by 360 degrees molecules are still in contact, the molecule is shifted in position in a random direction, and the rotation process is conducted again.
4. If a final structure with no contacts is reached, the configuration is saved to a .gro file.

## Rigorous Algorithm
1. Determines which atomistic moieties map to a CG bead.
2. Looping through each CG bead in each residue, places atomistic moiety such that the center of mass (COM) of the atomistic moiety matches the position of the CG bead.
3. Energy minimization of the atomistic residue.
4. Saves final configuration to .gro file format.
