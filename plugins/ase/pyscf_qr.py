import os
import string
import numpy as np
from ase.io import write
from ase.units import kcal
from ase.units import  mol as unit_mol
from ase.units import Hartree, Bohr
from ase.calculators.general import Calculator
import copy
import pyscf
from pyscf import gto,scf,grad,dft

class Pyscf(Calculator):
    def __init__(self, label='ase',basis='321g',charge=0,spin=0,method="hf", **kwargs):
        mol = pyscf.gto.Mole()
        mol.charge = charge
        mol.spin = spin
        mol.output = None#label+'.out'
        mol.chkfile = label+'.chk'
        mol.verbose = 1
        self.mol = mol
        self.method = method
        # set user values
        self.set(**kwargs)
         
    def set(self, **kwargs):
        """
        Sets the parameters on the according keywords
        Raises RuntimeError when wrong keyword is provided
        """
        for key in kwargs:
            if key is "charge":
            	self.mol.charge = kwargs[key]
            if key is "spin":
                self.mol.spin = kwargs[key]
            if key is "method":
                self.method = kwargs[key]
            if key is "basis":
                self.mol.basis = kwargs[key]
	    	
    def get_version(self):
        return self.version

    def initialize(self, atoms):
        pass

    def run(self):
            mol = self.mol
            mol.atom = [[atom.symbol, atom.position] for atom in self.atoms]
            mol.build()
            if self.method == 'hf':
                mscf = pyscf.scf.RHF(mol)
                mscf.max_cycle=200
                mscf.max_memory=1000
                self.command = mscf.run(conv_tol=1e-8,conv_tol_grad=1e-12).apply(pyscf.grad.RHF)
            if self.method == 'dft':
                mdft =  pyscf.dft.RKS(mol)
                mdft.max_cycle=200
                mdft.max_memory=1000
                self.command = mdft.run(conv_tol=1e-8,conv_tol_grad=1e-12,xc='bp86').apply(pyscf.grad.RKS)
            result = self.command.run()
            self.energy = result._scf.e_tot*((Hartree)/(kcal / unit_mol))
            self.energy_zero = self.energy
            self.energy_free = self.energy
            grad = np.array(result.grad())
            self.forces = grad*(-(Hartree/Bohr)/(kcal / unit_mol))         

    def read_energy(self, fname):
        return self.energy

    def read_forces(self, fname):
        return self.forces
        
    def atoms_are_equal(self, atoms_new):
        ''' (adopted from jacapo.py)
        comparison of atoms to self.atoms using tolerances to account
        for float/double differences and float math.
        '''
    
        TOL = 1.0e-6  # angstroms

        # check for change in cell parameters
        test = len(atoms_new) == len(self.atoms)
        if test is not True:
            return False
        
        # check for change in cell parameters
        test = (abs(self.atoms.get_cell() - atoms_new.get_cell()) <= TOL).all()
        if test is not True:
            return False
        
        old = self.atoms.arrays
        new = atoms_new.arrays
        
        # check for change in atom position
        test = (abs(new['positions'] - old['positions']) <= TOL).all()
        if test is not True:
            return False
        
        # passed all tests
        return True

    def update(self, atoms_new):
        if not self.atoms_are_equal(atoms_new):
            self.atoms = atoms_new.copy()
            self.run()
