import os.path
import glob
import reverse_mapping

def load_defaults():
    for filename in glob.glob('*.mol2'):
        name = filename.split(".")[0]
        file = os.path.join(os.path.dirname(__file__), filename)
        molecule = 

