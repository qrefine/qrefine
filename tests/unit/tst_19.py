from __future__ import division

import time, os
import iotbx.pdb
import libtbx
from libtbx import easy_run

qrefine = libtbx.env.find_in_repositories("qrefine")
qr_unit_tests_data = os.path.join(qrefine,"tests","unit","data_files")

stats = {'h_altconf.pdb' :
         [
           [{' N  ': 14, ' O  ': 17, ' H  ': 57, ' C  ': 42}],
           [{' N  ': 14, ' O  ': 17, ' H  ': 2, ' C  ': 42}],
          ],
         'h_altconf_2.pdb' :
         [
           [{' N  ': 14, ' O  ': 17, ' H  ': 57, ' C  ': 42}],
           [{' N  ': 14, ' O  ': 17, ' H  ': 2, ' C  ': 42}],
          ],
         }

def run1():
  for key, item in stats.items():
    print '-'*80
    print key
    print item
  for i, (pdb_filename, results) in enumerate(stats.items()):
    for j, complete_cap in enumerate(['model_completion',
                                      'capping',
                                    ]):
      #if not j: continue
      pdb_filename = os.path.join(qr_unit_tests_data, pdb_filename)
      cmd = 'qr.finalise %(pdb_filename)s action=%(complete_cap)s' % locals()
      cmd += ' keep_alt_loc=True'
      cmd += ' skip_validation=True'
      print cmd
      rc = easy_run.go(cmd)
      print dir(rc)
      rc.show_stdout()
      rc.show_stderr()
      assert rc.return_code==0, 'qr.finalise failed'
      fn = ['complete', 'capping'][j]
      pdb_inp = iotbx.pdb.input(
        file_name=os.path.basename(pdb_filename).replace('.pdb',
                                                         '_%s.pdb' % fn),
      )
      ph = pdb_inp.construct_hierarchy()
      oc = ph.overall_counts()
      oc.show()
      print '-'*80
      assert oc.element_charge_types==results[j][0], '''
      Failed to match %s %s
      ''' % (oc.element_charge_types, results[j][0])

if(__name__ == "__main__"):
  run1()
