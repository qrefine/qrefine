from __future__ import print_function
# LIBTBX_SET_DISPATCHER_NAME qr.development.qm_to_phenix_pdb
import os, sys

from iotbx import pdb
from qrefine.plugin.ase import gaussian_qr

'''
Optimization completed.
    -- Stationary point found.
                           ----------------------------
                           !   Optimized Parameters   !
                           ! (Angstroms and Degrees)  !
                           GradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGrad

                          Input orientation:
 ---------------------------------------------------------------------
 Center     Atomic      Atomic             Coordinates (Angstroms)
 Number     Number       Type             X           Y           Z
 ---------------------------------------------------------------------
      1          7           0       -3.647003    3.153409    0.276666
      2          1           0       -4.511792    2.982202    0.802624
'''

def get_input_orientation(lines, start):
  rc = []
  start = lines.find('Input orientation:', start)
  for i, line in enumerate(lines[start:].splitlines()):
    if i<5: continue
    if line.find('-----')>-1: break
    tmp = line.split()
    rc.append((float(tmp[-3]), float(tmp[-2]), float(tmp[-1])))
  return rc

def read_cartesian_coordinates(filename):
  f=open(filename, 'r')
  lines=f.read()
  f.close()
  #read_input_orientation=False
  coords = []
  start=0
  for i, line in enumerate(lines.splitlines()):
    start+=len(line)
    if line.find('-- Stationary point found.')>-1:
      rc = get_input_orientation(lines, start)
      coords.append(rc)
  return coords

def main(input_filename, master):
  print(input_filename, master)
  coords = read_cartesian_coordinates(input_filename)
  pdb_inp = pdb.input(master)
  hierarchy = pdb_inp.construct_hierarchy()
  hierarchy.show()
  for i in range(len(coords)):
    assert len(hierarchy.atoms())==len(coords[i])
    for atom, new in zip(hierarchy.atoms(), coords[i]):
      print(atom.quote(), atom.xyz, new)
      atom.xyz = new
    output_filename='%s_%05d.pdb' % (os.path.splitext(input_filename)[0], i+1)
    print(output_filename)
    hierarchy.write_pdb_file(output_filename)

if __name__ == '__main__':
  main(*tuple(sys.argv[1:]))
