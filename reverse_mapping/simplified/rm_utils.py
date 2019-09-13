import mdtraj as md
import mbuild as mb
import numpy as np

def load_original(filename, scaling_factor=1, cg_unit_conversion=0.6):
    conf = md.load(filename)
    conf = conf.atom_slice(conf.top.select('not name water'))

    conf.xyz *= cg_unit_conversion
    conf.unitcell_lengths *= cg_unit_conversion

    if np.prod(np.array(scaling_factor).shape) == 1:
        conf.xyz *= scaling_factor
        conf.unitcell_lengths *= 1.0
    elif np.prod(np.array(scaling_factor).shape) == 3:
        conf.xyz = conf.xyz * scaling_factor[np.newaxis, np.newaxis, :]
        conf.unitcell_lengths = (conf.unitcell_lengths * 
                                    scaling_factor[np.newaxis, :])

    conf = conf.center_coordinates()

    return conf


