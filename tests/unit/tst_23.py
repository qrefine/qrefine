from __future__ import absolute_import
import os, sys
from qrefine.tests.unit import run_tests
from qrefine.charges import charges_class

def get_charge(fn, assert_correct_chain_terminii=True):
  cc=charges_class(fn)
  return cc.get_total_charge(
    assert_correct_chain_terminii=assert_correct_chain_terminii)

pdbs = {
  'ACY' : {
    'ph_7' : '''
CRYST1  100.846  117.207  186.414  90.00  90.00  90.00 P 1
HETATM  305  C   ACY B 202       9.141  -2.833  -6.329  1.00 19.63           C
HETATM  306  O   ACY B 202       9.864  -1.976  -5.776  1.00 19.27           O
HETATM  307  CH3 ACY B 202       9.093  -2.854  -7.828  1.00 17.69           C
HETATM  308  OXT ACY B 202       8.460  -3.679  -5.708  1.00 16.11           O
HETATM  309  H1  ACY B 202       9.837  -2.163  -8.224  0.00 17.69           H
HETATM  310  H2  ACY B 202       9.303  -3.862  -8.186  0.00 17.69           H
HETATM  311  H3  ACY B 202       8.102  -2.547  -8.163  0.00 17.69           H
''',
      },
  'GLU' : {
    'polymer' : '''
CRYST1  155.738  152.261  119.955  90.00  90.00  90.00 P 1
ATOM     64  N   GLU A  66      26.383  32.549   0.663  1.00 28.08           N
ATOM     65  CA  GLU A  66      26.929  32.850  -0.611  1.00 30.46           C
ATOM     66  C   GLU A  66      25.993  32.476  -1.705  1.00 32.90           C
ATOM     67  O   GLU A  66      26.243  32.768  -2.863  1.00 33.11           O
ATOM     68  CB  GLU A  66      28.337  32.264  -0.790  1.00 29.74           C
ATOM     69  CG  GLU A  66      29.263  32.540   0.443  1.00 32.04           C
ATOM     70  CD  GLU A  66      29.956  33.866   0.313  1.00 32.83           C
ATOM     71  OE1 GLU A  66      30.149  34.362  -0.823  1.00 34.77           O
ATOM     72  OE2 GLU A  66      30.257  34.427   1.497  1.00 32.63           O
ATOM     73  H   GLU A  66      26.672  31.858   1.087  1.00 28.08           H
ATOM     74  HA  GLU A  66      27.030  33.814  -0.659  1.00 30.46           H
ATOM     75  HB2 GLU A  66      28.268  31.303  -0.905  1.00 29.74           H
ATOM     76  HB3 GLU A  66      28.750  32.664  -1.571  1.00 29.74           H
ATOM     77  HG2 GLU A  66      28.727  32.554   1.251  1.00 32.04           H
ATOM     78  HG3 GLU A  66      29.939  31.847   0.499  1.00 32.04           H
''',
    'terminii' : '''
CRYST1  155.738  152.261  119.955  90.00  90.00  90.00 P 1
ATOM      1  N   GLU A  66      26.383  32.549   0.663  1.00 28.08           N
ATOM      2  CA  GLU A  66      26.929  32.850  -0.611  1.00 30.46           C
ATOM      3  C   GLU A  66      25.993  32.476  -1.705  1.00 32.90           C
ATOM      4  O   GLU A  66      26.243  32.768  -2.863  1.00 33.11           O
ATOM      5  CB  GLU A  66      28.337  32.264  -0.790  1.00 29.74           C
ATOM      6  CG  GLU A  66      29.263  32.540   0.443  1.00 32.04           C
ATOM      7  CD  GLU A  66      29.956  33.866   0.313  1.00 32.83           C
ATOM      8  OE1 GLU A  66      30.149  34.362  -0.823  1.00 34.77           O
ATOM      9  OE2 GLU A  66      30.257  34.427   1.497  1.00 32.63           O
ATOM         OXT GLU A  66      24.955  31.866  -1.446  1.00 32.90           O
ATOM     10  H   GLU A  66      26.672  31.859   1.087  1.00 28.08           H
ATOM         H2  GLU A  66      26.531  33.237   1.224  1.00 28.08           H
ATOM         H3  GLU A  66      25.497  32.412   0.585  1.00 28.08           H
ATOM     11  HA  GLU A  66      27.030  33.814  -0.659  1.00 30.46           H
ATOM     12  HB2 GLU A  66      28.268  31.303  -0.905  1.00 29.74           H
ATOM     13  HB3 GLU A  66      28.750  32.664  -1.571  1.00 29.74           H
ATOM     14  HG2 GLU A  66      28.727  32.554   1.251  1.00 32.04           H
ATOM     15  HG3 GLU A  66      29.939  31.847   0.499  1.00 32.04           H
''',
    'mixed' : '''
CRYST1  155.738  152.261  119.955  90.00  90.00  90.00 P 1
ATOM      1  N   GLU A  66      26.383  32.549   0.663  1.00 28.08           N
ATOM      2  CA  GLU A  66      26.929  32.850  -0.611  1.00 30.46           C
ATOM      3  C   GLU A  66      25.993  32.476  -1.705  1.00 32.90           C
ATOM      4  O   GLU A  66      26.243  32.768  -2.863  1.00 33.11           O
ATOM      5  CB  GLU A  66      28.337  32.264  -0.790  1.00 29.74           C
ATOM      6  CG  GLU A  66      29.263  32.540   0.443  1.00 32.04           C
ATOM      7  CD  GLU A  66      29.956  33.866   0.313  1.00 32.83           C
ATOM      8  OE1 GLU A  66      30.149  34.362  -0.823  1.00 34.77           O
ATOM      9  OE2 GLU A  66      30.257  34.427   1.497  1.00 32.63           O
ATOM         OXT GLU A  66      24.955  31.866  -1.446  1.00 32.90           O
ATOM     10  H   GLU A  66      26.672  31.859   1.087  1.00 28.08           H
ATOM         H2  GLU A  66      26.531  33.237   1.224  1.00 28.08           H
ATOM     11  HA  GLU A  66      27.030  33.814  -0.659  1.00 30.46           H
ATOM     12  HB2 GLU A  66      28.268  31.303  -0.905  1.00 29.74           H
ATOM     13  HB3 GLU A  66      28.750  32.664  -1.571  1.00 29.74           H
ATOM     14  HG2 GLU A  66      28.727  32.554   1.251  1.00 32.04           H
ATOM     15  HG3 GLU A  66      29.939  31.847   0.499  1.00 32.04           H
''',
    'capping' : '''
CRYST1  155.738  152.261  119.955  90.00  90.00  90.00 P 1
ATOM      1  N   GLU A  66      26.383  32.549   0.663  1.00 28.08           N
ATOM      2  CA  GLU A  66      26.929  32.850  -0.611  1.00 30.46           C
ATOM      3  C   GLU A  66      25.993  32.476  -1.705  1.00 32.90           C
ATOM      4  O   GLU A  66      26.243  32.768  -2.863  1.00 33.11           O
ATOM      5  CB  GLU A  66      28.337  32.264  -0.790  1.00 29.74           C
ATOM      6  CG  GLU A  66      29.263  32.540   0.443  1.00 32.04           C
ATOM      7  CD  GLU A  66      29.956  33.866   0.313  1.00 32.83           C
ATOM      8  OE1 GLU A  66      30.149  34.362  -0.823  1.00 34.77           O
ATOM      9  OE2 GLU A  66      30.257  34.427   1.497  1.00 32.63           O
ATOM     10  H   GLU A  66      26.672  31.858   1.087  1.00 28.08           H
ATOM         H2  GLU A  66      26.532  33.237   1.224  1.00 28.08           H
ATOM     11  HA  GLU A  66      27.030  33.814  -0.659  1.00 30.46           H
ATOM     12  HB2 GLU A  66      28.268  31.303  -0.905  1.00 29.74           H
ATOM     13  HB3 GLU A  66      28.750  32.664  -1.571  1.00 29.74           H
ATOM     14  HG2 GLU A  66      28.727  32.554   1.251  1.00 32.04           H
ATOM     15  HG3 GLU A  66      29.939  31.847   0.499  1.00 32.04           H
ATOM         HC  GLU A  66      25.150  31.981  -1.494  1.00 32.90           H
''',
  },
}

results = {'ACY':{'ph_7': -1},
           'GLU':{'polymer'  : -3,
                  'terminii' : -1,
                  'capping'  : -1,
                  'mixed'    : -2,
                  },
           }

def run(prefix):
  for code, item in pdbs.items():
    for action, lines in item.items():
      #print code, action,
      fn = '%s_%s.pdb' % (code, action)
      f=open(fn, 'w')
      f.write(lines)
      f.close()

      rc = get_charge(fn,
                      assert_correct_chain_terminii=False,
                      #verbose=1,
      )
      ans = None
      level1 = results.get(code, None)
      if level1:
        ans = level1.get(action, None)
      if ans is not None:
        assert ans==rc, 'calculated charge %d not equal to expected %s' % (
            rc,
            ans,
            )
        os.remove(fn)

if(__name__=='__main__'):
  prefix = os.path.basename(__file__).replace(".py","")
  run_tests.runner(function=run, prefix=prefix, disable=False)
