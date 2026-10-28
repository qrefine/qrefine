from __future__ import absolute_import, division, print_function
# LIBTBX_SET_DISPATCHER_NAME qr.finalise
import sys
from iotbx.cli_parser import run_program
from mmtbx.programs import reduce2

note = """\
qr.finalise now runs mmtbx.reduce2 qmready=True. The former qr.finalise
options (action, calculate_charge, reduce, options.neutron,
hydrogen_atom_occupancies, append_to_end_of_model, skip_validation,
stop_if_poor_model) are gone: the total charge is always calculated and put
into the output file name (<prefix>H_charge<q>.pdb). All other arguments are
passed to mmtbx.reduce2 (for example capping=True); run mmtbx.reduce2
directly for its full option list.
"""

def run(args):
  print(note)
  return run_program(program_class=reduce2.Program,
                     args=list(args) + ['qmready=True'])

if __name__ == '__main__':
  run(sys.argv[1:])
