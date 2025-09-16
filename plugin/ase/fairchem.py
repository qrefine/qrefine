from fairchem.core import pretrained_mlip
from fairchem.core import FAIRChemCalculator as FAIRCalc
from ase.calculators.general import Calculator
import numpy as np
import ase.units as ase_units
import torch
from fairchem.core.units.mlip_unit.api.inference import InferenceSettings
settings = InferenceSettings(
    tf32=True,
    activation_checkpointing=False,
    merge_mole=False,
    compile=False,
    wigner_cuda=False,
    external_graph_gen=False,
    internal_graph_gen_version=2,
)


device = 'cuda' if torch.cuda.is_available() else 'cpu'
predictor = pretrained_mlip.get_predict_unit(
    "uma-sm",
     device=device,
     inference_settings=settings
     )
# calc = FAIRCalc(predictor,task_name='omol')

class FAIRChemCalculator2(Calculator):
    def __init__(self):
        # device = 'cuda' if torch.cuda.is_available() else 'cpu'
        # predictor = pretrained_mlip.get_predict_unit("uma-sm", device=device)
        self.label = None
        self.calc = FAIRCalc(predictor,task_name='omol')

    def run_qr(self, atoms, coordinates, charge, pointcharges, define_str=None):

        self.atoms = atoms
        atoms.info.update({"spin": 1, "charge": int(charge)})
        atoms.calc = self.calc
        unit_convert = ase_units.kcal / ase_units.mol
        self.energy_free = atoms.get_potential_energy()*unit_convert
        self.forces = atoms.get_forces().astype(np.float64)*unit_convert

        if 1: # we need debugging flag here to switch on and off.
            print(('Energy:',self.energy_free))
            print(('Force:', self.forces))

    def set_label(self, label):
        self.label = label

    # def set_charge(self,charge):
    #     self.atoms.info.update({"spin": 1, "charge": int(charge)})

    # def set_method(self, method):
    #   self.method = method

# class FAIRChemCalculator(FAIRCalc(predictor,task_name='omol')):
#     """
#      Wrapper around FAIRChemCalculator
#     """
#     def run_qr(self, atoms, coordinates, charge, pointcharges, define_str=None):
#         """
#         This method is called every time an energy and forces are needed.
#         The Q|R code calls this method at each step of LBFGS.
#         Args:
#         atoms (ase.atoms.Atoms) an updated set of atoms.
#         TODO: coordinates (numpy.ndarray) an updated set of coordinates. are these even being used?
#         charge (int) the charge on the molecular system.
#         TODO: pointcharges (?) are these even being used?
#         command (str) not used in this calculator
#         define_str() not used in this calculator, legacy from Turbomole?
#         """

#         self.atoms = atoms
#         self.atoms.calc = self
#         self.atoms.set_charge(charge)
#         unit_convert = ase_units.kcal / ase_units.mol
#         self.energy_free = self.atoms.get_potential_energy()*unit_convert
#         self.forces = self.atoms.get_forces().astype(np.float64)*unit_convert

#         if 0: # we need debugging flag here to switch on and off.
#             print(('Energy:',self.energy_free))
#             print(('Force:', self.forces))