# LIBTBX_SET_DISPATCHER_NAME qr.restraint
from __future__ import absolute_import, division, print_function
from iotbx.cli_parser import run_program
from qrefine import qr

# qr.restraint: one gradient evaluation of the restraints (QM by default)
# and one minimizer step, as a quick check that the restraints work on a
# model. Same arguments as qr.refine.
default_args = ['mode=opt', 'number_of_micro_cycles=1',
                'max_iterations_refine=1', 'clustering=False']

if __name__ == '__main__':
  run_program(program_class=qr.Program, args=default_args,
              allow_default_args=True)
