import numpy as np
import ase
import time
from ase import units
from ase.calculators.calculator import Calculator, all_changes
import requests


class RestAPICalculator(Calculator):
    """ASE calculator for a RestAPI server.
    Arguments:
    url: <hostname:port> of the server
    Server API requirement:
    /calc: accepts atomic numbers and atomic coordinate lists as json and return energy and forces (POST)
    The user is responsible to load the ML model at the server before running QR!
    ToDo:
     - charges
     - pointcharges
     - energy-only (current server does not support it)
    """

    implemented_properties = ["energy", "forces"]

    def __init__(self, url=None):
        super().__init__()
        self.url = url
        self.charge = 0
        self.spin = 1
        self.nproc = 1

    def calculate(self, atoms=None, properties=["energy", "forces"], system_changes=all_changes):
        atoms_json = {
            "species": atoms.get_atomic_numbers().tolist(),
            "coordinates": atoms.get_positions().tolist(),
            "spin": self.spin,
            "charge": self.charge,
            "ncpu": self.nproc
        }
        try:
            response = requests.post(url=f"{self.url}/calc", json=atoms_json)
            if response.status_code != int(200):
                raise Exception(f"Server-side error in calculation: \n {response.status_code}: {response.text}\n")
        except requests.exceptions.HTTPError as error:
            print(error)
        self.results["energy"] = np.array(response.json()["energy"])
        if "forces" in properties:
            self.results["forces"] = np.array(response.json()["forces"])

    def run_qr(self, atoms,coordinates=None,charge=None,pointcharges=None,define_str=None):
        """ Expects energy and forces in kcal/mol and kcal/mol/A """
        self.atoms = atoms
        self.charge = int(charge)
        # unit_convert = ase.units.kcal / ase.units.mol
        unit_convert = 1
        self.calculate(atoms, properties=["energy", "forces"])
        self.energy_free = self.results["energy"] * unit_convert
        self.forces = self.results["forces"].astype(np.float64) * unit_convert

    def set_nproc(self, nproc):
      self.nproc = int(nproc)