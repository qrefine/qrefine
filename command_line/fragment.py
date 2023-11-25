from __future__ import absolute_import, division, print_function
# LIBTBX_SET_DISPATCHER_NAME qr.fragment

from iotbx.cli_parser import run_program
from qrefine import fragmentation

if __name__ == '__main__':
  run_program(program_class=fragmentation.Program)
