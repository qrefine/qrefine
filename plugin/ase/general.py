"""
Minimal base class for the legacy Q|R ASE calculator plugins.

This is a verbatim copy of ``ase.calculators.general.Calculator`` from
ASE 3.22.1 (https://gitlab.com/ase/ase, LGPL-2.1-or-later), with the
unused ``read`` helper dropped.  ASE deprecated that module for years and
removed it in ASE 3.23, so the Q|R plugins that subclass it (MOPAC, ORCA,
TeraChem, Turbomole, PySCF, TorchANI, ANI, TensorMol, FAIRChem) now import
it from here instead of from ASE.

Q|R itself never calls the ASE calculator interface on these objects: it
calls ``run_qr`` and the ``set_*`` helpers and then reads the ``energy_free``
and ``forces`` attributes the plugins populate.  The methods below are kept
only so the plugins keep behaving exactly as they did under ASE 3.22.
"""
from ase.calculators.calculator import PropertyNotImplementedError


class Calculator:
    def __init__(self):
        return

    def set_atoms(self, atoms):
        self.atoms = atoms.copy()

    def get_atoms(self):
        atoms = self.atoms.copy()
        atoms.calc = self
        return atoms

    def get_name(self):
        """Return the name of the calculator (string)."""
        return self.name

    def get_version(self):
        """Return the version of the calculator (string)."""
        raise NotImplementedError

    def get_potential_energy(self, atoms, force_consistent=False):
        self.update(atoms)
        if force_consistent:
            return self.energy_free
        else:
            return self.energy_zero

    def get_forces(self, atoms):
        self.update(atoms)
        return self.forces

    def get_stress(self, atoms):
        self.update(atoms)
        if self.stress is not None:
            return self.stress
        else:
            raise PropertyNotImplementedError

    def initialize(self, atoms):
        """Prepare the input files required to start the program."""
        raise NotImplementedError
