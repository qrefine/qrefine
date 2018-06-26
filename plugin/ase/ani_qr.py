import os
import sys
from ani.ase_interface import ANI
import numpy as np
from ase.calculators.general import Calculator



class Ani(Calculator):
    def __init__(self,label="ase",atoms=None,coordinates='tmp_ase.pdb',**kwargs):
        self.label=label
        coordinates=os.path.dirname(label)+"/"+ coordinates
        self.coordinates=coordinates
        self.atoms = atoms
        self.energy_free = None
        self.forces = []
	self.ani = ANI()

    def run_qr(self,atoms,coordinates,charge,pointcharges,command=None,define_str=None):
        print " RUNNING ANI"
 	self.atoms=atoms
	self.coordinates=coordinates
        self.charge=charge
        self.pointcharges=pointcharges
 	self.command=command
	self.define_str=define_str
	mol = atoms
	mol.set_calculator(self.ani)
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
