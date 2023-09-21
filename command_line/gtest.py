# LIBTBX_SET_DISPATCHER_NAME qr.gtest
from iotbx.cli_parser import run_program
from qrefine import gtest

if __name__ == '__main__':
  run_program(program_class=gtest.Program)
