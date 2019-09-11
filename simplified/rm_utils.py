import mdtraj as md
import mbuild as mb
import numpy as np

def load_original(filename, scaling_factor=1):
    configuration = md.load(filename)
    configuration = configuration.atom_slice(
                        ref.top.select('not name water'))
    configuration.xyz *= .6 
    configuration.unitcell_lengths *= .6

    configuration.xyz[:,:,:2] *= 1.0
    configuration.unitcell_lengths[:,:2] *= 1.0

    configuration.center_coordinates()

    return configuration
