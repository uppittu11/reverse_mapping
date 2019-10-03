import os
import glob
import numpy as np
import mbuild as mb
import mdtraj as md

class Molecule(object):
    """
    Class to store information about each molecule type.

    Stores information about the indices for the heads/tails to be used
    to classify lipid moleules. Also stores the location of atomistic
    prototypes to be used later

    Arg:
    -----------
    name : string, optional, default=""
        The name of the molecule
    cg_head : int, optional, default=0
        The index for the head bead of the CG lipid
    aa_head : int, optional, default=0
        The index for the head bead of the atomistic lipid
    tail : list, optional, default=[0]
        The indices for the terminal tail beads of the CG lipid
    atomlist : list, optional, default=[]
        A list of atoms in the CG molecule. Can be used to identify
        CG atoms
    prototype : str, optional, default=""
        The location of the atomistic prototype structure
    
    """

    def __init__(self, name="", head=0, tail=0, atomlist=[],
                    prototype=""):
        self._name = name
        self._head = head
        self._tail = tail
        self._atomlist = atomlist
        self._prototype = prototype

    @property
    def name(self):
        return self._name
    
    @name.setter
    def name(self, name):
        assert type(name) == str
        self._name = name
    
    @property
    def head(self):
        return self._head
    
    @head.setter
    def head(self, index):
        index = int(index)
        self._head = index
    
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
    def prototype(self):
        return self._prototype
    
    @prototype.setter
    def prototype(self, compound):
        assert type(compound) == mb.Compound
        self.prototype = mb.clone(compound)