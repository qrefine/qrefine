"""
TensorMol is a Neural Network Augmented with Long-Range Physics.

*** Warning: it appears that this project is abandoned:
             https://github.com/jparkhill/TensorMol/issues/39

"""

from __future__ import print_function
from __future__ import absolute_import

import os, sys
import numpy as np
from TensorMol import *
from ase.calculators.general import Calculator

os.environ["CUDA_VISIBLE_DEVICES"] = ""  # set to use CPU

class TensormolCalculator(Calculator):

    def __init__(self):
        self.label =  None
        self.atoms = None
        self.coordinates = None
        self.atoms = MSet("morphine", center_=False)
        self.atoms.ReadXYZ("morphine")
        self.manager = self.GetChemSpiderNetwork(self.atoms, False)  # load chemspider network
        self.m = self.atoms.mols[0]
        if(0):
            print("Atoms {}  ".format(m.atoms))
            print()
            print("Coords {} ".format( m.coords))

    def run_qr(self,m):
        e, g = self.EnAndForce(m.coords, DoForce=True)
        if (0):
            print("Energy: {} ".format(e))
            print()
            print("Gradient: {}".format(g.tolist()))
        return e,g

    def get_command(self):
        return "TensorMol"

    def set_label(self,label):
        self.label = label

    def set_method(self,method):
        self.method = method

    # Make wrapper functions for energy, force and dipole
    def EnAndForce(self, x_, DoForce=True):
        mtmp = Mol(self.m.atoms, x_)
        Etotal, Ebp, Ebp_atom, Ecc, Evdw, mol_dipole, atom_charge, gradient = self.manager.EvalBPDirectEEUpdateSingle(mtmp,
                                                                                                                 PARAMS[
                                                                                                                     "AN1_r_Rc"],
                                                                                                                 PARAMS[
                                                                                                                     "AN1_a_Rc"],
                                                                                                                 PARAMS[
                                                                                                                     "EECutoffOff"],
                                                                                                                 True)
        energy = Etotal[0]
        force = gradient[0]
        if DoForce:
            return energy, force
        else:
            return energy

    EnergyForceField = lambda x: EnAndForce(x)

    def GetChemSpiderNetwork(self, a, Solvation_=False):
        TreatedAtoms = np.array([1, 6, 7, 8], dtype=np.uint8)
        PARAMS["tf_prec"] = "tf.float64"
        PARAMS["NeuronType"] = "sigmoid_with_param"
        PARAMS["sigmoid_alpha"] = 100.0
        PARAMS["HiddenLayers"] = [2000, 2000, 2000]
        PARAMS["EECutoff"] = 15.0
        PARAMS["EECutoffOn"] = 0
        PARAMS["Elu_Width"] = 4.6  # when elu is used EECutoffOn should always equal to 0
        PARAMS["EECutoffOff"] = 15.0
        PARAMS["AddEcc"] = True
        PARAMS["KeepProb"] = [1.0, 1.0, 1.0, 0.7]
        PARAMS["KeepProb"] = [1.0, 1.0, 1.0, 0.7]
        d = MolDigester(TreatedAtoms, name_="ANI1_Sym_Direct",
                        OType_="EnergyAndDipole")  # Initialize a digester that apply descriptor for the fragme
        tset = TensorMolData_BP_Direct_EE_WithEle(a, d, order_=1, num_indis_=1, type_="mol", WithGrad_=True)
        if Solvation_:
            PARAMS["DSFAlpha"] = 0.18
            manager = TFMolManage("chemspider12_solvation", tset, False,
                                  "fc_sqdiff_BP_Direct_EE_ChargeEncode_Update_vdw_DSF_elu_Normalize_Dropout", False,
                                  False)
        else:
            PARAMS["DSFAlpha"] = 0.18 * BOHRPERA
            manager = TFMolManage("chemspider12_nosolvation", tset, False,
                                  "fc_sqdiff_BP_Direct_EE_ChargeEncode_Update_vdw_DSF_elu_Normalize_Dropout", False,
                                  False)
        return manager
