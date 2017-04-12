from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME qr.help
import argparse

def keywords():
  """ Commands in qrefine:
      - start (including restart, check{dry run, and 1SCF, qmready} )
      - pause
      - stop
      - show
      - help
      - test
      - example
  """

def commands():
  """ Options and keywords in qrefine:
      - qm_calculator
      - macro_cycles
      - micro_cycles
      - max_bond_rmsd
      - refine_sites
      - refine_adp
      - cluster_qm
      - charge_embedding
      - cluster
  """

def settings():
  """ Settings in qrefine:
      - paths
      - ??
  """

def run():
  """ Help for qrefine:
       please use either --commands, --keywords, or --settings for more info
  """

if (__name__ == "__main__"):
  parser = argparse.ArgumentParser(description='qr.help')
  parser.add_argument('--commands', action='store_true',
                                    default=False,
                                    help='display the set of available commands ')
  parser.add_argument('--keywords', action='store_true',
                                    default=False,
                                    help='display the set of available keywords ')
  parser.add_argument('--settings', action='store_true',
                                    default=False,
                                    help='display the set of default settings   ')
  args = parser.parse_args()

  if(args.commands):   print commands.__doc__
  if(args.keywords):   print keywords.__doc__
  if(args.settings):   print settings.__doc__
  if(not args.commands and not args.keywords and not args.settings):
    print run.__doc__
