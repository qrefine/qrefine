"""
AIMNet2 plugin for Q|R.
Based on the AIMNet2 ASE calculator by the authors of the AIMNet2 model.
"""

from __future__ import print_function
from __future__ import absolute_import
import numpy as np
import ase.units as ase_units
from aimnet2calc import AIMNet2ASE



class AIMNet2Calculator(AIMNet2ASE):
    """ Modification of the AIMNet2ASE class to work with Q|R.
    """
    def run_qr(self, atoms, coordinates, charge, pointcharges, define_str=None):
        """
        This method is called every time an energy and forces are needed.
        The Q|R code calls this method at each step of LBFGS.
        Args:
        atoms (ase.atoms.Atoms) an updated set of atoms.
        TODO: coordinates (numpy.ndarray) an updated set of coordinates. are these even being used?
        charge (int) the charge on the molecular system.
        TODO: pointcharges (?) are these even being used?
        command (str) not used in this calculator
        define_str() not used in this calculator, legacy from Turbomole?
        """
        atoms.calc = self
        self.atoms = atoms
        self.set_charge(charge)
        unit_convert = ase_units.kcal / ase_units.mol
        self.calculate(atoms, properties=['energy', 'forces'])
        self.energy_free = self.results['energy'] * unit_convert
        self.forces = self.results['forces'].astype(np.float64) * unit_convert

        if 0: # we need debugging flag here to switch on and off.
            print(("AIMNet2: ",self.method))
            print(('Energy:',self.energy_free))
            print(('Force:', self.forces))
