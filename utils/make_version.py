"""
 get current qrefine version and update init file

"""
import os
from subprocess import Popen,PIPE,STDOUT
import libtbx.load_env

def main(init_filename=None):
  if init_filename is None:
    qrefine_path = libtbx.env.find_in_repositories("qrefine")
    init_filename = os.path.join(qrefine_path, '__init__.py')
  git_command='git describe --abbrev=6 --dirty --always --tags'
  shell = Popen([git_command], shell=True, stderr=PIPE,stdout=PIPE)
  out, err = shell.communicate()
  print 'Current version is %s ' % (out.rstrip())
  version='__version__=\" %s  \"'% (out.rstrip())

  with open(init_filename,'w') as outfile:
    outfile.write(version)

if __name__ == '__main__':
  main()
