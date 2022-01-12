from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
import time
from .test_reg_01_finalise  import test_finalise
from .test_reg_02_chunk     import test_chunk
from .test_reg_03_restraint import test_restraint
from .test_reg_04_refine    import test_refine

def regression_tests():
  return  {
    " finalise"      : test_finalise() ,
    " chunk"         : test_chunk(),
    " restraints"    : test_restraint(),
    " refinement"    : test_refine(),
    }

def run():
  for name, regression_test in regression_tests().items():
     t0 = time.time()
     print("Testing: {} ".format(name))
     regression_test.run()
     print("Total time: %6.2f"%(time.time()-t0))

if(__name__ == "__main__"):
  run()
