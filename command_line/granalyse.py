from __future__ import division
from __future__ import print_function
# LIBTBX_SET_DISPATCHER_NAME qr.granalyse
import numpy as np
import os
import sys
from qrefine.utils.mathbox import get_grad_mad, get_grad_angle
from libtbx.utils import Usage
import iotbx.pdb




def get_help():
  raise Usage("""
    qr.granalyse analyses gradients obtained from 'qr.refine mode=gtest' runs.
    Will write the gradient or difference gradient into the PDB.
    The reference gradient is found automatically if not explicitly specified.

    Examples:
    i)  qr.granalyse model.pdb  
    ii) qr.granalyse model.pdb --ref 3-15.npy

    Options:
      --ref <npy files> (set reference gradient)
      --occ write to occupancy field (instead of beta) 
      --grad (write gradient instead of difference gradient into pdb)
      --help  (print this help)
    """)
  sys.exit(0)
  return  


def id_file(filename):
  g_mode=int(filename[0])
  clustersize=int(filename[2:].split('.')[0])
  return g_mode, clustersize

def get_deviations(log,ref_grad,grad):
  import ase.units as ase_units
  '''
  calculates several deviation metrics for two gradient vectors
  '''
  # convert from unit less to kcal/mol
  unit_convert = ase_units.mol/ase_units.kcal # ~ 23.06
  ref_grad=ref_grad/unit_convert
  grad=grad/unit_convert
  ref_max=max(abs(i) for i in ref_grad)
  ref_min=min(abs(i) for i in ref_grad)
  ref_gnorm=np.linalg.norm(ref_grad)
  # grad=np.array(grad)
  
  print(' d(angle)  %f'  %(get_grad_angle(grad,ref_grad) ),file=log)
  print(' d(gnorm)  %f'  %(abs(np.linalg.norm(grad)-ref_gnorm) ),file=log)
  print(' d(max_g)  %f'  %(abs(max(abs(i) for i in grad)-ref_max) ),file=log)
  print(' d(min_g)  %f'  %(abs(min(abs(i) for i in grad)-ref_min) ),file=log)
  print(' RMSD      %f'  %(rmsd(grad,ref_grad)),file=log)
  print(' MAD       %f'  %(get_grad_mad(grad,ref_grad) ),file=log) 
  print(' ',file=log)

def get_grad_delta(ref_grad,grad):
  return np.abs(ref_grad-grad)

def get_grad_wdelta(ref_grad,grad):
  # weighted delta: delta_i= (g_i - g_i^ref)|g_i^ref|
  # *100 would be a 'percentage'
  dim3=int(ref_grad.shape[0])
  dim=int(dim3/3)
  d=np.zeros(dim)
  ref=np.reshape(ref_grad,(3,dim))
  g=np.reshape(grad,(3,dim))
  for i in range(dim):
    inorm=np.linalg.norm(ref[0:2,i])
    for j in range(3):
      d[i]+=np.abs((g[j,i]-ref[j,i])/inorm)
  return d*100/3

def sorting_weight(name):
  # apply custom weight based on gtest ID to allow easy sorting
  g_mode,size=id_file(name)
  weight=g_mode*1e4+size
  return weight

def set_ph_field(ph,val,field):
    # sets val to respective PDB field.
    nat_hierachy=len(ph.hierarchy.atoms())
    nat_numpy=val.shape[0]
    assert (nat_numpy-nat_hierachy==0),'wrong number of atoms!'

    count=0
    if 'occ' in field.lower():
      for a in ph.hierarchy.atoms():
        a.set_occ(val[count])
        count+=1
    elif 'beta' in field.lower():
      for a in ph.hierarchy.atoms():
        a.set_b(val[count])
        count+=1
    else:
      assert mode>0,'wrong value in set_ph_field function'

def atomic_mean_deviation(x):
  nat=int(x.shape[0]/3)
  new=np.zeros(nat)
  for i in range(nat):
    ii = 0 + 3 * i
    new[i]=(x[ii]+x[ii+1]+x[ii+2])/3
  return new

def rmsd(V, W):
    dim = V.shape[0]
    d = 0.0
    for i in range(dim):
        d += (V[i] - W[i])**2.0
    return np.sqrt(d/dim)

#######################################
def run(args,log):
  field='beta'
  ref_name=None
  do_grad=False
  do_delta=False
  do_wdelta=True
  do_expan=False

  args = sys.argv[1:]
  if ('--ref') in args:
      ref_name=args[args.index('--ref')+1]
  if ('-ref') in args:
      ref_name=args[args.index('-ref')+1]
  if ('--occ') in args:
      field='occ'
  if ('--grad') in args:
      do_grad=True
  if ('--delta') in args: # expert debug option
      do_delta=True
  if ('--help') in args:
      get_help()
  if (len(args)<1):
    get_help()
  
  pdbname=args[0]
  print('Reading: %s' %(pdbname))
  pdb_inp = iotbx.pdb.input(file_name=pdbname)
  pdb_obj = iotbx.pdb.hierarchy.input(file_name=pdbname)
  ph = pdb_inp.construct_hierarchy()
  nat=len(pdb_obj.hierarchy.atoms())
  

  # vectors
  ref_grad=None
  # idl=[]
  iname=[]
  igrad=[]

  # find all npy
  files=[]
  for f in os.listdir("."):
      if f.endswith(".npy"):
        files.append(f)

  # sort & label
  sfiles=sorted(files,key=sorting_weight)
  n_files=len(sfiles)
  print('Found the following files:',file=log)
  print(sfiles,file=log) 
  print('',file=log)
  for i in range(0,n_files):
    iname.append(sfiles[i])
    igrad.append(np.load(iname[i]))
  #set or read reference gradient
  if ref_name is None:
    if '0-0.npy' in iname:
      print('Found supercell calc')
      ref_idx  = 0
    else:
      ref_idx  = n_files-1 # last gradient should have best parameters
    ref_name = iname[ref_idx]
    ref_grad = np.array(igrad[ref_idx])
  else:
    ref_grad=np.load(ref_name)
  print('reference gradient taken from: %s \n' %(ref_name),file=log)

  # loop over available gradients and compare
  for i in range(n_files):
    print('     ~g_mode - max_res:', iname[i][:-4],file=log)
    get_deviations(log,ref_grad,igrad[i])
    if do_delta:
      dg=np.abs(ref_grad-igrad[i])
      delta=atomic_mean_deviation(dg)
    elif do_grad:
      delta=atomic_mean_deviation(igrad[i])
    elif do_wdelta:
      delta=get_grad_wdelta(ref_grad,igrad[i])
    set_ph_field(ph=pdb_obj,val=delta,field=field)
    output_pdb = "%s.pdb" %(iname[i][:-4])
    pdb_obj.hierarchy.write_pdb_file(file_name=output_pdb)


  return 0


if (__name__ == "__main__"):
  log = sys.stdout
  run(args=sys.argv[1:], log=log)




