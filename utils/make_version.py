from __future__ import print_function
"""
 get current qrefine version and update init file

	style: major.minor-patch level

	Manually release a new (full) version: git tag -a v1.0 -m 'version 1.0'
	patch level will be added automatically: git describe --abbrev=4 --dirty --always --tags

	E.g. after 1 commit this gives:
	v1.0-1-gab2f, where 1 indicates how many commits have been made and the rest is the git hash.

  The version string will be updated during the bootstrap install.
  For development updates please execute: phenix.python utils/make_version.py

"""

import os
from subprocess import Popen,PIPE,STDOUT

def main(init_filename=None):
  if init_filename is None:
    qrefine_path = os.path.abspath('.')
    assert os.path.basename(qrefine_path)=='qrefine'
    init_filename = os.path.join(qrefine_path, '__init__.py')
  git_command='git describe --abbrev=6 --dirty --always --tags'
  shell = Popen([git_command], shell=True, stderr=PIPE,stdout=PIPE)
  out, err = shell.communicate()
  print('Current version is %s ' % (out.rstrip()))
  version='__version__=\" %s  \"'% (out.rstrip())

  with open(init_filename,'w') as outfile:
    outfile.write(version)

if __name__ == '__main__':
  main()
