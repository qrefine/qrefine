import numpy as np
import math

def get_grad_mad(grad,ref):
  dim=grad.size
  diff=ref-grad
  s=sum(np.absolute(diff))
  return s/dim

def unit_vec(v):
  # much faster than np.linalg.norm because safetly checks
  l=math.sqrt((v[0])**2+(v[1])**2+(v[2])**2)
  return [v[0]/l,v[1]/l,v[2]/l]

def np_angle(v1, v2):
    v1_u=unit_vec(v1)
    v2_u=unit_vec(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))

def get_grad_angle(grad,ref):
  ''' mean angle of all gradient components in deg'''
  import numpy as np
  g=np.reshape(grad,(3,-1))
  r=np.reshape(ref,(3,-1))
  dim=g.shape[1]
  angle=0.0
  for i in range(dim):
    angle+=np.degrees(np_angle(g[:,i],r[:,i]))
  return angle/dim