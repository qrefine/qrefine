import sys
import pickle
from libtbx import easy_pickle
from fragment_utils import energy_gradients_from_qm

def qsub(index,qm_engine_name, charge, pdb_file, system_size, basis,
          charge_file=None,clustering_info_file="./ase/tmp/clustering_info_pickle"):
    unserialized_data = easy_pickle.load(clustering_info_file)
    fragment_atoms = unserialized_data["fragment_atoms"]
    cluster_atoms = unserialized_data["cluster_atoms"]
    e_and_g = energy_gradients_from_qm(i=int(index),qm_engine_name=qm_engine_name, qm_charge=charge,
                                       pdb_file=pdb_file, system_size=int(system_size), charge_file=charge_file,
                                       cluster_atoms=cluster_atoms, fragment_atoms=fragment_atoms, basis=basis)
    return_data = {}
    return_data["energy"] = e_and_g[0]
    return_data["gradients"] = e_and_g[1]
    return_file = './ase/tmp/'+ str(index)+".pickle"
    # Store data (serialize)
    easy_pickle.dump(file_name=return_file, obj=return_data)

if __name__=="__main__":
  args = sys.argv[1:]
  del sys.argv[1:]
  qsub(*tuple(args))
