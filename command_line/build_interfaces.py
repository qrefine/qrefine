"""

1. This command line tool checks to see that the correct version of java is installed.
Java 1.8 is needed to run the density based descriptor code. This is called Yoink-0.0.1.jar and is
installed in the plugin/yoink folder.

2. This command line tool also checks to see if environment variables are pointing to installed qm engines.

3. A conditional import of python based engines.

A list of installed engines are printed to the screen for your convenience.


"""

# LIBTBX_SET_DISPATCHER_NAME qr.build_interfaces
from __future__ import print_function
import os, sys

from libtbx import easy_run

def draw_box_around_text(msg, width=78, log=None):
  def _draw_lines(msg, terminating_char):
    rc = ''
    for line in msg:
      outl = '%s' % terminating_char
      assert len(line)<width-2
      outl += ' %s ' % line
      outl += '%s' % (' '*(width-len(line)-2))
      outl += '%s' % terminating_char
      rc += '%s\n' % outl
    return rc[:-1]
  ##############################
  if log is None: log=sys.stdout
  ascii_128 = False
  try:
    # assert 0
    print('\n%s' % ('%s%s%s' % (u'\u250C',
                                u'\u2500'*width,
                                u'\u2510',
                               )
                  ),
          file=log,
          )
  except:
    ascii_128 = True
    top_line = '%s' % ('%s%s%s' % ('|',
                                   '-'*width,
                                   '|',
                                  )
                      )
  if ascii_128:
    print(top_line, file=log)
    print(_draw_lines(msg, '|'), file=log)
    print(top_line, file=log)
  else:
    print(_draw_lines(msg, u'\u2502'), file=log)
    print('%s' % ('%s%s%s' % (u'\u2514',
                              u'\u2500'*width,
                              u'\u2518',
                             )
                  ),
         file=log,
         )

def run():
  msg = ['',
         'Welcome to the Q|R interface checker',
         '',
         'For more information: https://github.com/qrefine/qrefine/wiki/Installation',
         ]
  draw_box_around_text(msg, width=78)

  cmd = 'java -version'
  rc = easy_run.go(cmd)
  for line in rc.stdout_lines:
    if line.find('java version')>-1:
      version = line.split('"')[1][:3]
      if float(version)<1.8:
        print('''
  Need at least Java 1.8. Please update your system Java using a JDK bundle
  and try again.
        ''')
        sys.exit()

  java_env_vars = {'JAVA_HOME' : 'absolute_path_of_java_home',
                   'JAVA_LIB_PATH' : 'absolute_path_of_java_lib',
                   'LD_LIBRARY_PATH' : '$LD_LIBRARY_PATH:$JAVA_LIB_PATH/server',
                   }
  count = 0
  for env_var in java_env_vars:
    print('  Set? "%s" "%s"' % (env_var, os.environ.get(env_var, False)))
    if not os.environ.get(env_var, False):
      print('''
      The following environment variables need setting.
      ''')
      for env_var, help in java_env_vars.items():
        print('%s %s : %s' % (' '*10, env_var, help))
      if env_var.startswith('JAVA'):
        print('''
      On OSX use
        /usr/libexec/java_home -v 1.8
      to find the install directory
      ''')
      elif env_var=='LD_LIBRARY_PATH':
        print('''
      The Phenix environment ignores the user set LD_LIBRARY_PATH unless the
      environmental variable PHENIX_TRUST_OTHER_ENV is set.
      ''')
      count+=1
      break
  else:
    print('\n  Java appears to be installed\n')
  if count:
    print('   STOPPING')
    #sys.exit()

  qm_engine_env_vars = {'MOPAC_COMMAND' : 'Mopac executable',
                        'TERACHEM_COMMAND' : 'TeraChem directory',
                        'Orca_COMMAND': 'Orca directory',
                        'XTBHOME':'XTB directory',
                        'g16root':'Gaussian 16 directory',
                        'TURBODIR ':'Turbomole directory',
                        }
  count = []
  for env_var, value in os.environ.items():
    print(' ~> %s : %s' % (env_var, value))
  for env_var in qm_engine_env_vars:
    if os.environ.get(env_var, False):
      print('\n  Environmental variable %s set for "%s" to "%s"\n' % (
        env_var,
        qm_engine_env_vars[env_var],
        os.environ[env_var],
        ))
      if not os.path.exists(os.environ[env_var]):
        print('''
        Environmental variable "%s" : "%s" does not point to anything!

        STOPPING
        ''' % (env_var, os.environ[env_var]))
        sys.exit()
      count.append(env_var)
    else:
      print('  Environmental variable %s for "%s" not found\n' % (
        env_var,
        qm_engine_env_vars[env_var],
        ))

  qm_engines_python =    { 'PyScf':' A collection of electronic structure programs powered by Python',
                           'ANI':'ANI-1 neural net potential with python interface (ASE)',
                           'TorchANI':'Accurate Neural Network Potential on PyTorch'

  }

  qm_engines_python_installed = {}
  for name,description in qm_engines_python.items():

    if name == 'PyScf':
      try:
        import pyscf
        print("  PyScf successfully imported")
        qm_engines_python_installed[name] = description
      except:
        print("  PyScf could not be imported")

    if name == 'TorchANI':
      try:
        import torch
        import torchani
        print("  TorchANI successfully imported")
        qm_engines_python_installed[name] = description
      except:
        print("  TorchANI could not be imported")

    if name == 'ANI':
      try:
        import ani
        from ani.ase_interface import aniensloader
        from ani.ase_interface import ANIENS
        print("  ANI successfully imported")
        qm_engines_python_installed[name] = description
      except:
        print("  ANI could not be imported")

  if count:
    print('''
    QM engines set
    ''')
    for env_var in count:
      print('%s %s : %s' % (' '*10, env_var, os.environ[env_var]))
    for name, description in qm_engines_python_installed.items():
      print('%s %s : %s' % (' '*10, name, description))
  else:
    print('''
    No QM engines found!

    Install and set an environmental variable from the list.
    ''')
    for env_var, help in qm_engine_env_vars.items():
      print('%s %s : %s' % (' '*10, env_var, help))
    print


  draw_box_around_text(['','Done!','','Run a Q|R job using qr.refine','','https://github.com/qrefine/qrefine/wiki/Installation'])

if __name__=="__main__":
  args = sys.argv[1:]
  del sys.argv[1:]
  run(*tuple(args))
