# Import pyNeuroChem
import os
import sys
#Set this path to fit your personal installation
sys.path.append('/home/USER/qr_ani/old/ASE_ANI/lib')
from ase_interface import ANI
from libtbx import easy_run
import numpy as np
import  ase
import time
import subprocess
from ase import units
from ase.io import read, write
from ase.optimize import BFGS, LBFGS
from ase.calculators.general import Calculator
from libtbx import easy_run

# Set required files for pyNeuroChem
#anipath  = '/home/jujuman/Dropbox/ChemSciencePaper.AER/ANI-c08e-ntwk'
#cnstfile = anipath + '/rHCNO-4.6A_16-3.1A_a4-8.params'
#saefile  = anipath + '/sae_6-31gd.dat'
#nnfdir   = anipath + '/networks/'

# Construct pyNeuroChem class
#nc = pync.molecule(cnstfile, saefile, nnfdir, 0)

class Ani(Calculator):
    def __init__(self,label="ase",atoms=None,coordinates='tmp_ase.pdb',**kwargs):
        self.label=label
        coordinates=os.path.dirname(label)+"/"+ coordinates
        self.coordinates=coordinates
        self.atoms = atoms
        self.energy_free = None
        self.forces = []

    def run_qr(self,atoms,coordinates,charge,pointcharges,command=None,define_str=None):
        print " RUNNING ANI"
 	self.atoms=atoms
	self.coordinates=coordinates
        self.charge=charge
        self.pointcharges=pointcharges
 	self.command=command
	self.define_str=define_str
	mol = atoms
	mol.set_calculator(ANI())
	print(mol)
	energy = mol.get_potential_energy()
	force = mol.get_forces()
	force = force.astype(np.float64)
	self.energy_free = energy
	self.forces = force
	print('Final Energy: ', energy)
	print('forces are : ',force)
	
  
    def get_command(self):
        return "ANI" 
   
    def set_label(self,label):
        self.label=label
