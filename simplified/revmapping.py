import mdtraj as md
import mbuild as mb
import numpy as np
import matplotlib.pyplot as plt

# Loads trajectory and scales it by CG scaling factor and relaxation scaling factor
ref_orig = md.load('final_wrapped.hoomdxml')

ref = ref_orig
ref = ref.atom_slice(ref.top.select('not name water'))
ref.xyz *= .6 
ref.unitcell_lengths *= .6

ref.xyz[:,:,:2] *= 1.0
ref.unitcell_lengths[:,:2] *= 1.0

ref.center_coordinates()
keep_residues = []
for lipid_head in ref.top.select('name mhead2 head chead'):
    if ref.xyz[-1, lipid_head,0] > 1 and ref.xyz[-1, lipid_head,1] > 1:
        keep_residues += [ref.top.atom(lipid_head).residue.index]
keep_residues = 'residue '+ ' '.join([str(i) for i in keep_residues])
ref = ref.atom_slice(ref.top.select(keep_residues))

# get residue names for each residue
# also get coordinates of head for each residue
lipid_coordinates = []

for residue in ref.top.residues:
    atom_names = [atom.name for atom in residue.atoms]
    if 'mhead2' in atom_names:
        if len(atom_names) == 17:
            angle = md.compute_angles(traj=ref, 
                                      angle_indices=[(residue.atom(0).index, 
                                                      residue.atom(9).index,
                                                      residue.atom(14).index)])
            if angle[0] > np.pi/2:
                residue.name = 'ucer2_extended'
            else:
                residue.name = 'ucer2_hairpin'
        elif len(atom_names) == 14:
            angle = md.compute_angles(traj=ref, 
                                      angle_indices=[(residue.atom(0).index, 
                                                      residue.atom(6).index,
                                                      residue.atom(11).index)])
            if angle[0] > np.pi/2:
                residue.name = 'ecer2_extended'
            else:
                residue.name = 'ecer2_hairpin'
    elif 'chead' in atom_names:
        residue.name = 'chol'
        
    elif 'head' in atom_names:
        if len(atom_names) == 9:
            residue.name = 'ffa24'
        elif len(atom_names) == 6:
            residue.name = 'ffa16'
    else:
        print("Unsupported molecule type")
        assert 1==0
        
for residue in ref.top.residues:
    # determine if lipid tails point up or down
    if 'cer' in residue.name:
        if ref.xyz[-1,residue.atom(0).index,2] > ref.xyz[-1,residue.atom(6).index,2]:
            orientation = 'up'
        else:
            orientation = 'down'
    if 'ffa' in residue.name:
        if ref.xyz[-1, residue.atom(0).index,2] > ref.xyz[-1,residue.atom(-1).index,2]:
            orientation = 'up'
        else:
            orientation = 'down'
    if 'chol' in residue.name:
        if ref.xyz[-1, residue.atom(0).index,2] < ref.xyz[-1,residue.atom(-3).index,2]:
            orientation = 'up'
        else:
            orientation = 'down'
    residue.name = residue.name + "-" + orientation
    
# get the coordinates and lipid type for each lipid
for residue in ref.top.residues:
    if 'ucer' in residue.name:
        lipid_coordinates += [(ref.xyz[-1,residue.atom(9).index], residue.name)]
    elif 'ecer' in residue.name:
        lipid_coordinates += [(ref.xyz[-1,residue.atom(6).index], residue.name)]
    elif 'ffa' in residue.name:
        lipid_coordinates += [(ref.xyz[-1,residue.atom(-1).index], residue.name)]
    elif 'chol' in residue.name:
        lipid_coordinates += [(ref.xyz[-1,residue.atom(0).index], residue.name)]
    else:
        print("Unsupported molecule type")
        assert 1==0

# define atomistic molecules
lipid_types = dict()
# hairpin ucer2
filename = '/Users/parashara/Documents/devel/forcefield/charmm-ff/ucer2_hairpin_charmm/ucer2.mol2'
structure = mb.load(filename)
structure.translate_to([0,0,0])
lipid_types.update({"ucer2_hairpin":[74, structure]})

# extended ucer2
filename = '/Users/parashara/Documents/devel/forcefield/charmm-ff/ucer2_extended/ucer2.mol2'
structure = mb.load(filename)
structure.translate_to([0,0,0])
lipid_types.update({"ucer2_extended":[74, structure]})

# ffa24
filename = '/Users/parashara/Documents/devel/forcefield/charmm-ff/c24ffa/ffa24.mol2'
structure = mb.load(filename)
structure.translate_to([0,0,0])
lipid_types.update({"ffa24":[23, structure]})

# chol
filename = '/Users/parashara/Documents/devel/forcefield/charmm-ff/cholesterol/chol.mol2'
structure = mb.load(filename)
structure.translate_to([0,0,0])
lipid_types.update({"chol":[0, structure]})

# water
filename = '/Users/parashara/Documents/devel/forcefield/charmm-ff/tip3p_pppm/SOL_new.mol2'
structure = mb.load(filename)
structure.translate_to([0,0,0])
solvent = mb.clone(structure)

bilayer = mb.Compound()

tilt = np.deg2rad(11)
tilt_direction = np.random.rand(3)
tilt_direction[2] = 0

wpl = 20
n_waters = ref.top.n_residues * wpl

solvent_density = 1000
solvent_mass = 18.01

for molecule_index, molecule in enumerate(lipid_coordinates):
    if molecule_index % 1 == 0:
        print('{} of {}'.format(molecule_index, len(lipid_coordinates)) + ' '*40, end='\r')
    lipid_name, orientation = molecule[1].split('-')
    ref_atom = lipid_types[lipid_name][0]
    lipid = mb.clone(lipid_types[lipid_name][1])
    lipid.translate_to([0,0,0])
    theta = np.random.rand()*2.0*np.pi
    lipid_final_pos = np.array(molecule[0])
    lipid.rotate(around=[0,0,1], theta=theta)
    if orientation == 'up':
        lipid.rotate(around=[1, 0, 0], theta=np.pi)
    lipid.rotate(around=tilt_direction, theta=tilt)
    lipid.translate(-lipid[ref_atom].pos)
    lipid.translate(lipid_final_pos)
    
    
    # re-orients if two atoms are in contact
    contact = False
    neighbors_pos = np.array([[0,0,0]])
    lipid_pos = lipid.xyz
    for child in bilayer.children:
        if np.linalg.norm(child.pos[:2]-lipid.pos[:2]) < 1.5:
            neighbors_pos = np.append(neighbors_pos, child.xyz, axis=0)
    if len(neighbors_pos) > 0:
        for i in lipid_pos:
            if np.sum(np.linalg.norm(neighbors_pos-i, axis=1) < 0.05) > 0:
                contact = True
        n_rot = 1
        n_wiggle = 0
        while contact:
            print('{} of {} -- wiggling'.format(molecule_index, len(lipid_coordinates)), end='\r')
            lipid.rotate(around=tilt_direction, theta=-tilt)
            lipid.translate_to([0,0,0])
            lipid.rotate(around=[0,0,1], theta=-theta)
            theta += np.deg2rad(5)
            lipid.rotate(around=[0,0,1], theta=theta)
            lipid.rotate(around=tilt_direction, theta=tilt)
            lipid.translate_to([0,0,0])
            lipid.translate(-lipid[ref_atom].pos)
            lipid.translate(lipid_final_pos)
            lipid_pos = lipid.xyz
            contact = False
            for i in lipid_pos:
                if np.sum(np.linalg.norm(neighbors_pos-i, axis=1) < 0.05) > 0:
                    contact = True
            n_rot += 1
            if n_rot > 360/5+1:
                lipid_final_pos += np.random.randint(2,size=3)*0.02
                n_rot = 0
                n_wiggle += 1
                if n_wiggle == 10:
                    print("Molecule {} doesn't fit!")
                    assert 1 == 2
    bilayer.add(lipid)
bilayer.translate_to([0, 0, 0])

#bilayer = mb.clone(quarter_bilayer)
bilayer.translate_to([0, 0, 0])
lipid_box = bilayer.boundingbox
# make solvent boxes
n_waters = 363 * 20
n_solvent = int(n_waters/4)
solvent_volume = np.cbrt(n_solvent*solvent_mass/(6.02e23)/solvent_density/1000) *10e9
solvent_z = solvent_volume / (lipid_box.lengths[0] * lipid_box.lengths[1]) * 1.2
solvent_box1 = mb.Box(mins=[lipid_box.mins[0], lipid_box.mins[1], lipid_box.maxs[2]],
                      maxs=[lipid_box.maxs[0], lipid_box.maxs[1], lipid_box.maxs[2]+solvent_z])
solvent_box2 = mb.Box(mins=[lipid_box.mins[0], lipid_box.mins[1], lipid_box.mins[2]-solvent_z],
                      maxs=[lipid_box.maxs[0], lipid_box.maxs[1], lipid_box.mins[2]])
print(solvent_volume)
print(solvent_z)
print(lipid_box.maxs[2]+solvent_z)
print(solvent_box1)
print(solvent_box2)

solvent1 = mb.fill_box(compound=solvent, n_compounds=n_solvent, box=solvent_box1, overlap=0.5)
print(2)
solvent2 = mb.fill_box(compound=solvent, n_compounds=n_solvent, box=solvent_box2, overlap=0.5)

# add all components
system = mb.Compound()
system.add(mb.clone(bilayer))
system.add(solvent1)
system.add(solvent2)

system.translate_to(system.boundingbox.lengths/2.0)

system.save('rev_mapped.hoomdxml', overwrite=True, box=mb.Box(mins=[0, 0, 0], maxs=bilayer.boundingbox.maxs+0.1))

#bilayer.save('rev_mapped_lipids.hoomdxml', overwrite=True, box=bilayer.boundingbox)
system_md = md.load('rev_mapped.hoomdxml')
print('loaded')
system_md.unitcell_lengths = system.boundingbox.lengths
system_md.xyz /= 10
system_md.unitcell_angles = [90, 90, 90]

cer_count = 0
chol_count = 0
ffa_count = 0
water_count = 0

for residue in system_md.top.residues:
    if 'CLN' in [atom.name for atom in residue.atoms]:
        residue.name = 'ucer2'
        cer_count += 1
    elif 'OCL' in [atom.name for atom in residue.atoms]:
        residue.name = 'ffa24'
        ffa_count += 1
    elif 'O1' in [atom.name for atom in residue.atoms]:
        residue.name = 'tip3p'
        water_count += 1
    else:
        residue.name = 'chol'
        chol_count += 1

print(cer_count, chol_count, ffa_count, water_count)
system_md.save('rev_mapped.gro')