from collections import OrderedDict

from warnings import warn

import mbuild as mb
from mbuild.compound import Compound
from mbuild.compound import clone
from mbuild.coordinate_transform import force_overlap
from copy import deepcopy
import multiprocessing as mp
import numpy as np

import pdb

__all__ = ['reverse_map']

def reverse_map(coarse_grained, mapping_moieties, target=None, solvent_name=None, sol_per_bead=4, sol_cutoff=2, scaling_factor=5, parallel=True):
    """ Reverse map an mb.Compound

    Parameters
    ---------
    coarse_grained : mb.Compound
        original structure. Generated from a MDTraj trajectory
    mapping_moieties : dictionary
        Relate CG molecule names to a list of finer-detailed mbuild
        Compounds. Care must be taken that bead indices match with
        list indices.
    target_structure : dictionary
        mb.Compound, optional, default=False
        A target atomistic structure which can be used to reconstruct
        bonding.
        Bond network in the reverse-mapped structure will be completely
        overridden by the bond network from the target atomistic structure
        Care must be taken that atom indices match perfectly

    """

    aa_system = Compound()

    not_solvent = [mol for mol in coarse_grained.children if mol.name != solvent_name]
    is_solvent  = [mol for mol in coarse_grained.children if mol.name == solvent_name]

    print("There are {} non-solvent molecules and {} solvent molecules.".format(len(not_solvent), len(is_solvent)))

    # For each bead, replace it with the appropriate mb compound
    # Iterate through each molecule (set of particles that are bonded together)
    if parallel:
        pool = mp.Pool(processes=mp.cpu_count())

        # get the solvent molecules mapped in parallel
        inp = zip(is_solvent,
                [target[solvent_name]]*len(is_solvent),
                [sol_per_bead]*len(is_solvent),
                [sol_cutoff]*len(is_solvent))
        chunksize = int(len(is_solvent) / mp.cpu_count()) + 1
        solvent_list = pool.starmap(reverse_map_solvent, inp, chunksize)
        # name the solvents

        # get the non_solvent molecules mapped in parallel
        inp = zip(not_solvent,
                [target]*len(not_solvent),
                [mapping_moieties]*len(not_solvent))
        chunksize = int(len(not_solvent) / mp.cpu_count()) + 1
        molecule_list = pool.starmap(reverse_map_molecule, inp, chunksize)


        # put put solvents in one list
        solvent_molecule_list = []
        for i in solvent_list:
            solvent_molecule_list += i

        # put lipids in a box and get the box size
        for molecule in molecule_list:
            aa_system.add(molecule)

        print(aa_system.boundingbox)

        # put everything in a box
        for molecule in solvent_molecule_list:
            aa_system.add(molecule)

    else:
        [aa_system.add(reverse_map_molecule(molecule, target, mapping_moieties)) for molecule in not_solvent]
        solvent_compound = reverse_map_solvent(is_solvent, target[solvent_name], sol_per_bead, sol_cutoff)
        [aa_system.add(molecule) for molecule in solvent_compound.children]


    return aa_system

def reverse_map_molecule(molecule, target, mapping_moieties):
    print('unique')
    cg_molecule = clone(molecule) # CG molecule
    aa_template = target[molecule.name] # full aa Compound for molecule
    aa_moieties = mapping_moieties[molecule.name] # list of lists of indices for each bead
    aa_molecule = Compound() # this will have the final aa molecule
    cg_to_aa = [] # list of tuples containing (real index, aa atom)

    # now cycle through beads
    for index, bead in enumerate(cg_molecule.particles()):
        aa_atoms = Compound() # placeholder for molecule atoms
        [aa_atoms.add(clone(aa_template.children[i])) for
            i in aa_moieties[index]]
        aa_atoms.translate_to(bead.pos) # shift to cg_bead position
        cg_to_aa += list(zip(aa_moieties[index], aa_atoms.children))

    # sort atoms in cg_to_aa and add them to the aa_molecule
    cg_to_aa = sorted(cg_to_aa)
    for atom in cg_to_aa:
        aa_molecule.add(clone(atom[1]))

    # add bonds from the template
    aa_template = aa_template.to_trajectory()
    for i,j in aa_template.top.bonds:
        aa_molecule.add_bond([aa_molecule[i.index], aa_molecule[j.index]])

    # equilibrate molecule and shift back to center
    # if the atom names match OpenBabel forcefield naming convention:
    try:
        aa_molecule.energy_minimization(steps=2500)

    # otherwise rename with just element names:
    except:
        atomnames = [i.name for i in aa_molecule] # get the atomnames
        for atom in aa_molecule: # make the atomnames elements
            atom.name=atom.name[0]

        aa_molecule.energy_minimization(steps=2500)

        for i, atom in enumerate(atomnames):
            aa_molecule[i].name = atomnames[i]

    #aa_molecule.translate(center)
    aa_molecule.name = molecule.name
    return aa_molecule

def reverse_map_solvent(cg_molecule, target, sol_per_bead=4, cutoff=2, scaling_factor=5):
    """
    molecules: list of water beads
    target: single atomistic solvent molecule
    sol_per_bead: number of atomistic solvent molecules per CG bead
    cutoff: max distance an atomistic molecule can be placed from the
        center of the CG bead
    """

    solvent = Compound() # will contain all the solvent molecules in a bead
    solvent_molecule = Compound() # placeholder for each single molecule

    # for each atomistic molcule
    print("unique")
    for i in range(sol_per_bead):
        # get a random vector by which to shift the atomistic molecule
        """
        randx = cutoff * (1 - 2 * np.random.rand())
        randy = np.sqrt(cutoff**2 - randx**2) * (1 - 2 * np.random.rand()) # randy bobandy
        randz = np.sqrt(cutoff**2 - randx**2 - randy**2) * (1 - 2 * np.random.rand())
        """
        randx = 0.0
        randy = 0.0
        randz = 0.0
        if i == 0:
            randx += 1.0
        elif i == 1:
            randx -= 1.0
        elif i == 2:
            randy += 1.0
        elif i == 3:
            randy -= 1.0

        shift_vec = np.array([randx, randy, randz])
        shift_vec *= 0.2
        #np.random.shuffle(shift_vec)

        # get random angles to spin the solvent molecule by
        theta = 2 * np.pi * np.random.rand()
        phi = np.pi * np.random.rand()

        # make a solvent molecule compound and shift it to the correct position
        solvent_molecule = clone(target)
        solvent_molecule.translate_to(cg_molecule.pos + shift_vec)
        solvent_molecule.spin(theta, [0,0,1])
        solvent_molecule.spin(phi, [1,0,0])

        # add the molecule to the bead compound
        solvent.add(solvent_molecule)

    # time to minimize the energy!
    # try with current atom names
    """
    try:
        solvent.energy_minimization(steps=500)

    # otherwise rename with just element names:
    except:
        atomnames = [i.name for i in solvent] # get the atomnames
        for atom in solvent: # make the atomnames elements
            atom.name=atom.name[0]

        solvent.energy_minimization(steps=500)

        for i, atom in enumerate(atomnames):
            solvent[i].name = atomnames[i]
    """
    # scale the solvent by 5


    # get a list of individual molecules (so that it separates them into residues)
    solvent = [clone(child) for child in solvent.children]

    for i, solvent_compound in enumerate(solvent):
        solvent[i].name = cg_molecule.name
        solvent[i].translate_to(np.array(solvent[i].pos)*scaling_factor)
    return solvent




"""

        new_molecule = Compound()
        # Rather than sort through the molecule, which may be unsorted
        # Look at the parent's particles, which will be sorted
        for bead in molecule[0].parent.particles():
            new_atom = clone(mapping_moieties[bead.name])
            cg_to_aa[bead] = new_atom
            new_atom.translate(bead.pos)
            new_molecule.add(new_atom)
        aa_system.add(new_molecule)

    # Go back and include bonds
    if target_structure:
        # If a target atomistic structure is provided, just its bond graph
        # to the reverse-mapped structure
        aa_system.root.bond_graph = None
        target_traj = target_structure.to_trajectory()

        for (i,j) in target_traj.topology.bonds:
            aa_system.add_bond([aa_system[i.index], aa_system[j.index]])

    else:
        # If no target atomistic structure is provided, look at each molecule,
        # working inwards from the ends of the molecule

        cg_bonds = list(coarse_grained.bonds())
        # Repeatedly iterate through the coarse grained bonds, but only bond
        # particles that have a certain number of available ports
        while len(cg_bonds) > 0:
            for p_i, p_j in cg_bonds:
                if 0 < len(cg_to_aa[p_i].available_ports()) <= 1 or \
                    0 < len(cg_to_aa[p_j].available_ports()) <= 1:
                            p_i_port, p_j_port = _find_matching_ports(cg_to_aa[p_i],
                                cg_to_aa[p_j])
                            force_overlap(cg_to_aa[p_i], from_positions=p_i_port,
                                to_positions=p_j_port)
                            cg_bonds.remove((p_i, p_j))


    # Put molecules back after energy minimization
    for cg_particle, aa_particles in cg_to_aa.items():
        aa_particles.translate_to(cg_particle.pos)
"""




def _find_matching_ports(i, j):
    """ Find corresponding ports on two mBuild compounds"""

    i_ports = i.available_ports()
    j_ports = j.available_ports()
    i_port_names = [p.name for p in i.available_ports()]
    j_port_names = [p.name for p in j.available_ports()]
    common_name = list(set(i_port_names).intersection(j_port_names))
    if len(common_name) != 1:
        warn("{} ports were found with corresponding names for"
                " particles {} and {}".format(len(common_name), i,j))
    i_port = [p for p in i.available_ports() if p.name == common_name[0]]
    j_port = [p for p in j.available_ports() if p.name == common_name[0]]
    #for j_port in j_ports:
        #if j_port.name == i_port.name:
            #return i_port, j_port
    return i_port[0], j_port[0]


