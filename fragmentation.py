from __future__ import division
from __future__ import print_function
import sys
import time
import os.path
import libtbx
import iotbx.pdb
from libtbx.program_template import ProgramTemplate
from qrefine import qr


from qrefine.fragment import fragments
from qrefine.fragment import get_qm_file_name_and_pdb_hierarchy
from qrefine.fragment import charge
from qrefine.fragment import write_mm_charge_file


class Program(ProgramTemplate):

  description = """
  qr.fragment break up a system into many small pieces


  Example:
  qr.fragment model.pdb  [<param_name>=<param_value>] ...
  """


  datatypes = ['model', 'phil', ]

  master_phil_str = qr.master_phil_str

  def validate(self):
    print('Validate inputs:', file=self.logger)
    self.data_manager.has_models(
      expected_n=1,
      exact_count=True,
      raise_sorry=True)


  def run(self):
    self.header("Fragment start")
    ph = self.data_manager.get_model().get_hierarchy()
    cs = self.data_manager.get_model().crystal_symmetry()

    fq = fragments(
      pdb_hierarchy=ph,
      crystal_symmetry=cs,
      charge_embedding=True,
      debug=True,
      qm_engine_name="terachem")
    print("Residue indices for each cluster:\n", fq.clusters, file=log)
    fq_ext = fq.get_fragment_extracts()
    for i in range(len(fq.clusters)):
        # add capping for the cluster and buffer
        print("capping frag:", i, file=log)
        get_qm_file_name_and_pdb_hierarchy(
                            fragment_extracts=fq_ext,
                            index=i)
        print("point charge file:", i, file=log)
        #write mm point charge file
        write_mm_charge_file(fragment_extracts=fq_ext,
                                        index=i)
