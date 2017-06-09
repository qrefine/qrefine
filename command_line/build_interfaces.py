# LIBTBX_SET_DISPATCHER_NAME qr.build_interfaces
import os, sys

from libtbx import easy_run

def run():
  width=78
  msg = ['',
         'Welcome to the Q|R inteface checker',
         '',
         ]
  print '\n%s' % ('%s%s%s' % (u'\u250C',
                              u'\u2500'*width,
                              u'\u2510',
                             )
                )
  for line in msg:
    outl = '%s' % u'\u2502'
    assert len(line)<width-2
    outl += ' %s ' % line
    outl += '%s' % (' '*(width-len(line)-2))
    outl += '%s' % u'\u2502'
    print outl
  print '%s' % ('%s%s%s' % (u'\u2514',
                            u'\u2500'*width,
                            u'\u2518',
                           )
                )

  cmd = 'java -version'
  rc = easy_run.go(cmd)
  for line in rc.stdout_lines:
    if line.find('java version')>-1:
      version = line.split('"')[1][:3]
      if float(version)<1.8:
        print '''
  Need at least Java 1.8. Please update your system Java using a JDK bundle
  and try again.
        '''
        sys.exit()

  java_env_vars = {'JAVA_HOME' : 'absolute_path_of_java_home',
                   'JAVA_LIB_PATH' : 'absolute_path_of_java_lib',
                   'LD_LIBRARY_PATH' : '$LD_LIBRARY_PATH:$JAVA_LIB_PATH/server',
                   }
  for env_var in java_env_vars:
    if not os.environ.get(env_var, False):
      print '''
      The following environment variables need setting.
      '''
      for env_var, help in java_env_vars.items():
        print '%s %s : %s' % (' '*10, env_var, help)
      print '''
      On OSX use
        /usr/libexec/java_home -v 1.8
      to find the install directory
      '''
      break
  else:
    print '\n  Java appears to be installed'

  qm_engine_env_vars = {'MOPAC_CMD' : 'Mopac executable',
                        }
  count = 0
  for env_var in qm_engine_env_vars:
    if os.environ.get(env_var, False):
      print '\n  Environmental variable %s set for "%s" to %s\n' % (
        env_var,
        qm_engine_env_vars[env_var],
        os.environ[env_var],
        )
      count+=1
    else:
      print '\n  Environmental variable %s for "%s" not found\n' % (
        env_var,
        qm_engine_env_vars[env_var],
        )
  if count:
    pass
  else:
    print '''
    No QM engines found!

    Install and set an environmental variable from the list.
    '''
    for env_var, help in qm_engine_env_vars.items():
      print '%s %s : %s' % (' '*10, env_var, help)
    print

if __name__=="__main__":
  args = sys.argv[1:]
  del sys.argv[1:]
  run(*tuple(args))
