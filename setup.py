from setuptools import setup, find_packages
import sys

try:
    import mdtraj
except ImportError:
    print('Building and running analysis requires mdtraj. See '
          'http://mdtraj.org/latest/installation.html for help!')
    sys.exit(1)

try:
    import mbuild
except ImportError:
    print('Building and running analysis requires mbuild. See '
          'https://mosdef.org/mbuild/ for help!')
    sys.exit(1)

setup(name='reverse_mapping',
      version='0.1',
      description=('Reverse mapping scripts to convert CG multilayer' +
                  ' lipid structures to atomistic configurations'),
      url='https://github.com/uppittu11/reverse_mapping',
      author='Parashara Shamaprasad',
      author_email='p.shama@vanderbilt.edu',
      license='MIT',
      packages=['reverse_mapping.simplified'],
      package_dir={'reverse_mapping': 'reverse_mapping/'},
      include_package_data=True,
      install_requires=["mdtraj"],
)