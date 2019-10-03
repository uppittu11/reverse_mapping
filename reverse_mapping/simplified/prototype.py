import os
import glob
import numpy as np
import mbuild as mb
import mdtraj as md

class Prototype(object):
    """
    Class to store information about each molecule type.

    Stores information about the indices for the heads/tails to be used
    to classify lipid moleules. Also stores the atomistic
    prototypes corresponding to the compound

    Args:
    -----------
    name : string
        The name of the molecule
    cg_head : int
        The index for the head bead of the CG lipid
    aa_head : int
        The index for the head bead of the atomistic lipid
    tail : int
        The index for the terminal tail bead of the CG lipid
    atomlist : list
        A list of atoms in the CG molecule. Can be used to identify
        CG atoms
    compound : mb.Compound
        The location of the atomistic prototype structure
    
    """

    def __init__(self, name, cg_head, aa_head, tail, atomlist, compound):
        self._name = name
        self._cg_head = cg_head
        self._aa_head = aa_head
        self._tail = tail
        self._atomlist = atomlist
        self._compound = compound

    @property
    def name(self):
        return self._name
    
    @name.setter
    def name(self, name):
        assert type(name) == str
        self._name = name
    
    @property
    def cg_head(self):
        return self._cg_head
    
    @cg_head.setter
    def cg_head(self, index):
        index = int(index)
        self._cg_head = index

    @property
    def aa_head(self):
        return self._aa_head
    
    @aa_head.setter
    def aa_head(self, index):
        index = int(index)
        self._aa_head = index

    @property
    def tail(self):
        return self._tail
    
    @tail.setter
    def tail(self, index):
        index = list(index)
        self._head = index
    
    @property
    def atomlist(self):
        return self._atomlist
    
    @atomlist.setter
    def atomlist(self, atomlist):
        assert type(atomlist) == list
        self._atomlist = atomlist
    
    @property
    def compound(self):
        return self._compound
    
    @compound.setter
    def compound(self, compound):
        assert isinstance(compound, mb.Compound)
        self._compound = compound

class CeramidePrototype(Prototype):
    """
    Class to store information about each ceramide molecule type.

    Stores information about the indices for the heads/tails to be used
    to classify lipid moleules. Also stores the atomistic
    prototypes corresponding to the compound. Allows 

    Args:
    -----------
    name : string
        The name of the molecule
    cg_head : int
        The index for the head bead of the CG lipid
    aa_head : int
        The index for the head bead of the atomistic lipid
    tail : int
        The index for the terminal tail bead of the CG lipid
    atomlist : list
        A list of atoms in the CG molecule. Can be used to identify
        CG atoms
    compound : mb.Compound
        The location of the atomistic prototype structure
    
    """
    @property
    def tail(self):
        return self._tail
    
    @tail.setter
    def tail(self, index):
        index = list(index)
        self._head = index