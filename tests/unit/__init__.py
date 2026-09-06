from __future__ import absolute_import, division, print_function
from libtbx import easy_run
import os, pathlib
from iotbx import reflection_file_reader
import mmtbx.model
import iotbx.pdb
from libtbx.utils import null_out
import qrefine
from qrefine import qr

class run_fmodel(object):
  def __init__(self, args, prefix):
    self.mtz     = "%s.mtz"%prefix
    self.std_out = "%s.log"%prefix
    args += ["output.file_name=%s.mtz"%prefix, "type=real", "label=F-obs",
             "r_free=0.01", "random_seed=2679941"]
    self.cmd = " ".join(["phenix.fmodel"] + args + [">%s"%self.std_out])
    assert easy_run.call(self.cmd)==0
    assert os.path.isfile(self.mtz)
    #
    miller_arrays = reflection_file_reader.any_reflection_file(file_name =
      self.mtz).as_miller_arrays()
    self.f_obs, self.flags = None, None
    for ma in miller_arrays:
      if("F-obs" in ma.info().label_string()):        self.f_obs = ma
      if("R-free-flags" in ma.info().label_string()):
        self.flags = ma.array(data = ma.data()==1)
    assert [self.f_obs, self.flags].count(None)==0

class run_qrefine(object):
  def __init__(self, args, prefix):
    self.log     = "%s_real_space_refined_000.log"%prefix
    self.std_out = "%s.log"%prefix
    #
    self.cmd = " ".join(
      ["qr.refine"] + args + [">%s"%self.std_out])
    print(self.cmd)
    assert easy_run.call(self.cmd)==0
    self.pdb = str(next(pathlib.Path('pdb').glob('*_refined.pdb'), None))
    pdb_inp = iotbx.pdb.input(file_name = self.pdb)
    self.model = mmtbx.model.manager(model_input = pdb_inp, log = null_out())
    self.model.setup_scattering_dictionaries(
      scattering_table = qr.get_default_params().scattering_table)
    #
    #log_lines = self.get_lines(fn = self.log)
    #buf_lines = self.get_lines(fn = self.std_out)
    #if([log_lines,buf_lines].count(None)==0 and
    #   len(log_lines) != len(buf_lines)):
    #  if 0:
    #    for line1, line2 in zip(log_lines, buf_lines):
    #      if line1!=line2:
    #        print(line1)
    #        print(line2)
    #        assert 0
    #  raise RuntimeError("Log and stdout are different!")

  def get_lines(self, fn):
    isf = os.path.isfile(fn)
    if(self.sorry_expected and not isf): return
    assert isf, fn
    result = []
    with open(fn,"r") as fo:
      for l in fo.readlines():
        l = l.strip("\n")
        result.append(l)
    return result
