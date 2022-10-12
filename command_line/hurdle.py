from __future__ import division, print_function
# LIBTBX_SET_DISPATCHER_NAME phenix.development.hurdle

import os, sys

from iotbx.cli_parser import run_program
from qrefine.hurdle import Program

# =============================================================================

if (__name__ == '__main__'):
  results = run_program(program_class=Program)
