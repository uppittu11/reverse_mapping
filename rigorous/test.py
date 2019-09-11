import reverse_mapping as revm
import mbuild as mb
import mdtraj as md
import time
import sys

sys.setrecursionlimit(10000)

def from_traj(compound, traj):
    atom_mapping = dict()
    for residue in traj.top.residues:
        res_compound = mb.compound.Compound()
        for atom in residue.atoms:
            new_atom = mb.Particle(name=str(atom.name),
                    pos=traj.xyz[-1, atom.index])
            res_compound.add(new_atom)
            atom_mapping[atom] = new_atom
        res_compound.name = '{0}'.format(residue.name)
        compound.add(res_compound)

    for mdtraj_atom1, mdtraj_atom2 in traj.topology.bonds:
        atom1 = atom_mapping[mdtraj_atom1]
        atom2 = atom_mapping[mdtraj_atom2]
        compound.add_bond((atom1, atom2))

    compound.periodicity = traj.unitcell_lengths[0]
    return compound


# Load in your CG system
traj = md.load('cg-traj.xtc', top='cg-traj.gro')[-1]
print('Loaded CG frame')
"""
# CG length conversion
traj.xyz *= .6
traj.unitcell_lengths *= .6
"""
# get rid of waters
#traj = traj.atom_slice(traj.top.select('name water'))

# Now select the residues/atoms i want to keep
atoms_i_want = []

# Put everything in the box
#anchor = traj.top.find_molecules()
traj = traj.center_coordinates()
#traj.image_molecules(inplace=True, anchor_molecules=anchor)
"""
# get the mhead beads
heads = traj.top.select('name " mhead2"')
waters = traj.top.select('name " water"')

# if the mhead is within the first quadrant we add it to the mix
for res, head in enumerate(heads):
    pos = traj.xyz[0,head,:]
    if (pos[0] < traj.unitcell_lengths[0,0]/16) and (pos[1] < traj.unitcell_lengths[0,1]/16):
        atoms_i_want += list(traj.top.select('residue {}'.format(res)))

for res, water in enumerate(waters):
    pos = traj.xyz[0,water,:]
    if (pos[0] < traj.unitcell_lengths[0,0]/16) and (pos[1] < traj.unitcell_lengths[0,1]/16):
        atoms_i_want += list(traj.top.select('index {}'.format(water)))
"""
print('Collected atoms')
"""
# cut out only the atoms i want to keep
traj = traj.atom_slice(atoms_i_want)

for ind in range(traj.top.n_atoms):
    traj.top.atom(ind).name = traj.top.atom(ind).name[0]

# re-center the molecules
traj = traj.center_coordinates()
#traj = traj.image_molecules(traj.top.find_molecules())
"""
# make our cg system into a mBuild compound
cg = mb.compound.Compound()
cg = from_traj(cg, traj)
cg.translate_to(pos=[0,0,0])

# Give names to the subcomponents
for index, subcompound in enumerate(cg.children):
    if subcompound.n_particles > 2:
        cg.children[index].name = 'ucer3'
    else:
        cg.children[index].name = 'water'

# Save original cg structure
resnames = []
for child in cg.children:
    if child.name == 'ucer3':
        resnames += ['c']
    if child.name == 'water':
        resnames += ['w']
cg.save('my_cg.gro', box=cg.boundingbox, residues=resnames)

# Load in atomistic target structures:
# 1. UCER3
ucer3_md = md.load('ucer3.gro')
ucer3 = mb.Compound()
ucer3.from_trajectory(ucer3_md)
for i in ucer3:
    if i.name[0] in {'C', 'N', 'O', 'P', 'S'}:
        pir = ucer3.particles_in_range(i, .16)
        for j in pir[1:]:
            ucer3.add_bond((i, j))

# 2. Water
water_md = md.load('water.gro')
water = mb.Compound()
water.from_trajectory(water_md)
for i in water:
    if i.name[0] in {'C', 'N', 'O', 'P', 'S'}:
        pir = water.particles_in_range(i, .16)
        for j in pir[1:]:
            water.add_bond((i, j))

# put the atomistic target structures into the dictionary
target = dict()
target['ucer3'] = ucer3
target['water'] = water

print('Loaded target structures')

# get the mapping moieties for these molecule types:
mapping_moieties = dict()

# 1. UCER3
mapping_moieties['ucer3'] =  [[67, 68, 69, 70, 71, 72, 73],
                              [58, 59, 60, 61, 62, 63, 64, 65, 66],
                              [49, 50, 51, 52, 53, 54, 55, 56, 57],
                              [40, 41, 42, 43, 44, 45, 46, 47, 48],
                              [31, 32, 33, 34, 35, 36, 37, 38, 39],
                              [22, 23, 24, 25, 26, 27, 28, 29, 30],
                              [13, 14, 15, 16, 17, 18, 19, 20, 21],
                              [4, 5, 6, 7, 8, 9, 10, 11, 12],
                              [0, 1, 2, 3],
                              [74, 75, 76, 77, 78, 81, 82],
                              [85, 86, 89, 90, 91, 92, 93, 94],
                              [95, 96, 97, 98, 99, 100, 101, 102, 103],
                              [104, 105, 106, 107, 108, 109, 110, 111, 112],
                              [113, 114, 115, 116, 117, 118, 119, 120, 121],
                              [122, 123, 124, 125, 126, 127, 128, 129, 130, 131],
                              [79, 80],
                              [83, 84],
                              [87, 88]]

# 2. Water
mapping_moieties['water'] = [[0, 1, 2]]

print('Starting reverse mapping on {} residues'.format(len(cg.children)))
# run reverse mapping and time it

start = time.time()
reverse_mapped = revm.reverse_map(coarse_grained=cg, mapping_moieties=mapping_moieties, target=target, solvent_name='water',
        sol_per_bead=4, sol_cutoff=2.0, parallel=True)
end = time.time()

# print and save.
print("reverse mapping took {} min or {} per residue.".format((end-start)/60, (end-start)/len(cg.children)))

reverse_mapped.translate_to(reverse_mapped.boundingbox.lengths / 2)
resnames = [child.name for child in reverse_mapped.children]
reverse_mapped.save('my_reverse_mapped.gro', box=reverse_mapped.boundingbox, residues=resnames)
