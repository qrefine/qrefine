# LIBTBX_SET_DISPATCHER_NAME qr.fetch
from __future__ import absolute_import, division, print_function

from iotbx.cli_parser import run_program
from qrefine import fetch

if __name__ == '__main__':
  run_program(program_class=fetch.Program)
