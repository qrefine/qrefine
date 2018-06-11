# LIBTBX_SET_DISPATCHER_NAME qr.development.status
import os, sys

results = {}

def process_file(filename):
  print 'process_file',filename
  f=file(filename, 'rb')
  lines = f.read()
  f.close()
  done=False
  for line in lines.splitlines():
    pass #rint line
    if line.find('***  Starting incremental Fock matrix formation  ***')>-1:
      done='Starting SCF'
    elif line.find('TOTAL RUN TIME:')>-1:
      done='FINISHED'
    elif line.find('== MOPAC DONE ==')>-1:
      done=line
  print line
  return done

def process_dir(d, engine_name='orca'):
  for filename in os.listdir(d):
    if not filename.endswith('.out'): continue
    rc = process_file(os.path.join(d, filename))
    results.setdefault((d,filename), None)
    results[(d,filename)] = rc

def check_output_pdbs(d):
  for filename in os.listdir(d):
    if filename.find('cycle.pdb')>-1:
      if filename.find('weight')>-1:
        results.setdefault('weight', [])
        results['weight'].append(filename)
      elif filename.find('refine')>-1:
        results.setdefault('refine', [])
        results['refine'].append(filename)

def run(cwd=None):
  if cwd is None:
    cwd=os.getcwd()
  else:
    os.chdir(cwd)
  for root, dirs, files in os.walk(cwd):
    if os.path.basename(root)=='ase_error':
      print 'ASE error directory exists'
    elif root.find('ase_error')>-1:
      pass
    elif root.find('ase')>-1:
      try:
        i = int(os.path.basename(root))
      except ValueError, e:
        i = None
      if i is not None:
        process_dir(root)

    if os.path.basename(root)=='pdb':
      check_output_pdbs(root)

  for key, item in sorted(results.items()):
    print key, item

if __name__=='__main__':
  run(*tuple(sys.argv[1:]))
