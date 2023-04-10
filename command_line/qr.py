# LIBTBX_SET_DISPATCHER_NAME qr.refine
from __future__ import absolute_import, division, print_function

from iotbx.cli_parser import run_program
from qrefine import qr

if __name__ == '__main__':
  run_program(program_class=qr.Program)
