from __future__ import division
from __future__ import print_function
# LIBTBX_SET_DISPATCHER_NAME qr.granalyse
import numpy as np
import os
import sys
from qrefine.utils.mathbox import get_grad_mad, get_grad_angle
from libtbx.utils import Usage


log = sys.stdout

def get_help():
  raise Usage("""
    grana.py analyses gradients obtained from 'qr.refine mode=gtest' runs.

    Examples:
    i)  grana.py  
    ii) grana.py -ref 3-15.npy

    Options:
    qr.refine --ref <npy files> (set reference gradient)
    qr.refine --help  (print this help)
    """)
  sys.exit(0)
  return  


def id_file(filename):
  g_mode=int(filename[0])
  clustersize=int(filename[2:].split('.')[0])
  return g_mode, clustersize

def get_input(cline):
  for arg in range(len(cline)):
    if '-ref' or '--ref' in arg:
      ref_name=cline[arg]
    if '-help' or '--help' or '-h' in arg:
      get_help()
  return 0

# def f_name(x):        
  # return "-".join(map(str,x))

def get_deviations(log,ref_grad,grad):
  '''
  calculates several deviation metrics for two gradient vectors
  '''
  ref_max=max(abs(i) for i in ref_grad)
  ref_min=min(abs(i) for i in ref_grad)
  ref_gnorm=np.linalg.norm(ref_grad)
  # grad=np.array(grad)

  print(' d(angle)  %f'  %(get_grad_angle(grad,ref_grad) ),file=log)
  print(' d(gnorm)  %f'  %(abs(np.linalg.norm(grad)-ref_gnorm) ),file=log)
  print(' d(max_g)  %f'  %(abs(max(abs(i) for i in grad)-ref_max) ),file=log)
  print(' d(min_g)  %f'  %(abs(min(abs(i) for i in grad)-ref_min) ),file=log)
  print(' MAD       %f'  %(get_grad_mad(grad,ref_grad) ),file=log) 
  print(' ',file=log)


def sorting_weight(name):
  # apply custom weight based on gtest ID to allow easy sorting
  g_mode,size=id_file(name)
  weight=g_mode*1e4+size
  return weight

#######################################
def run(args,log):

  ref_name=None
  log = sys.stdout
  args = sys.argv[1:]
  get_input(args)

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
  if ref_grad is None:
    ref_idx= n_files-1 # last gradient should have best parameters
    ref_name=iname[ref_idx]
    ref_grad=np.array(igrad[ref_idx])
  else:
    ref_grad=np.load(ref_name)
  print('reference gradient taken from: %s \n' %(ref_name),file=log)

  # loop over available gradients and compare
  for i in range(n_files):
    print('     ~g_mode - max_res:', iname[i][:-4],file=log)
    get_deviations(log,ref_grad,igrad[i])

  return 0


if (__name__ == "__main__"):
  run(args=sys.argv[1:], log=log)




