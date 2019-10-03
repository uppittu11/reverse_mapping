import numpy as np
import mbuild as mb

class Molecule(object):
    """
    Class to store information about each molecule type.

    Arg:
    -----------
    name : string
        The name of the molecule
    pos : arraylike, shape=(3)
        The position of the molecule. This is typically where the head
        of the lipid is placed.
    prototype : mb.Compound
        The atomistic compound representing the molecule
    
    """

    def __init__(self, name, pos, prototype):
        def _validate_inputs():
            assert isinstance(name, str)
            pos = np.array(pos)
            assert isinstance(pos, np.array())
            assert pos.shape[-1] == 3
            assert isinstance(prototype, mb.Compound)

        _validate_inputs()
        self._name = name
        self._pos = np.array(pos)
        self._prototype = mb.clone(prototype)

    @property
    def name(self):
        return self._name
    
    @name.setter
    def name(self, name):
        assert isinstance(name, str)
        self._name = name
    
    @property
    def pos(self):
        return self._pos
    
    @pos.setter
    def pos(self, pos):
        assert isinstance(pos, np.array())
        assert pos.shape[-1] == 3
        self._pos = np.array(pos)
    
    @property
    def prototype(self):
        return self._prototype
    
    @prototype.setter
    def prototype(self, compound):
        assert type(compound) == mb.Compound
        self.prototype = mb.clone(compound)  