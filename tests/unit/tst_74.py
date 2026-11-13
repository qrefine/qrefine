from __future__ import absolute_import, division, print_function
import os
import iotbx.pdb
from qrefine import conformers

# A residue with two residue names and a shared blank atom (1ejg PRO A /
# SER B, C 22 as deposited, one N) is read as two residue groups (PRO: blank
# and A; SER: B and C). The residue groups of one residue share their labels:
# conformer A takes the N and PRO A alone, not SER B next to it. (With the N
# split per conformer iotbx reads one residue group of three atom groups.)

pdb_str = """\
CRYST1   50.000   50.000   50.000  90.00  90.00  90.00 P 1
ATOM    400  N   THR A  21       4.060   9.406  -4.225  1.00  1.66           N
ATOM    401  CA  THR A  21       3.672  10.572  -3.450  1.00  1.42           C
ATOM    402  C   THR A  21       4.921  11.354  -2.986  1.00  1.58           C
ATOM    403  O   THR A  21       5.842  10.726  -2.484  1.00  2.24           O
ATOM    404  CB  THR A  21       2.930  10.114  -2.187  1.00  1.52           C
ATOM    405  OG1 THR A  21       1.959   9.127  -2.549  1.00  1.89           O
ATOM    406  CG2 THR A  21       2.278  11.275  -1.460  1.00  2.12           C
ATOM    414  N   PRO A  22       4.893  12.687  -3.112  1.00  1.83           N
ATOM    415  CA APRO A  22       6.065  13.428  -2.599  0.57  1.81           C
ATOM    416  C  APRO A  22       6.374  13.133  -1.164  0.57  1.66           C
ATOM    417  O  APRO A  22       5.479  13.001  -0.339  0.57  2.08           O
ATOM    418  CB APRO A  22       5.642  14.893  -2.754  0.57  2.86           C
ATOM    419  CG APRO A  22       4.664  14.854  -4.051  0.57  2.66           C
ATOM    420  CD APRO A  22       3.965  13.504  -3.905  0.57  2.27           C
ATOM    429  CA BSER A  22       5.996  13.423  -2.690  0.21  1.55           C
ATOM    430  C  BSER A  22       6.363  13.088  -1.209  0.21  2.92           C
ATOM    431  O  BSER A  22       5.421  13.041  -0.353  0.21  1.87           O
ATOM    432  CB BSER A  22       5.436  14.841  -2.858  0.21  2.60           C
ATOM    433  OG BSER A  22       4.763  15.242  -1.644  0.21  2.11           O
ATOM    440  CA CSER A  22       6.038  13.629  -2.643  0.22  2.31           C
ATOM    441  C  CSER A  22       6.403  13.244  -1.155  0.22  1.92           C
ATOM    442  O  CSER A  22       5.611  12.721  -0.295  0.22  1.77           O
ATOM    443  CB CSER A  22       5.623  15.088  -2.692  0.22  3.56           C
ATOM    444  OG CSER A  22       6.750  15.783  -2.279  0.22  6.38           O
ATOM    450  N   GLU A  23       7.666  13.200  -0.875  1.00  1.86           N
ATOM    451  CA  GLU A  23       8.098  12.996   0.522  1.00  1.81           C
ATOM    452  C   GLU A  23       7.395  13.972   1.459  1.00  1.66           C
ATOM    453  O   GLU A  23       7.067  13.587   2.585  1.00  1.96           O
ATOM    454  CB  GLU A  23       9.623  13.168   0.673  1.00  2.12           C
ATOM    455  CG AGLU A  23      10.456  12.141  -0.052  0.58  2.73           C
ATOM    456  CD AGLU A  23      11.939  12.222   0.243  0.58  4.60           C
ATOM    457  OE1AGLU A  23      12.418  13.328   0.542  0.58  7.56           O
ATOM    458  OE2AGLU A  23      12.632  11.204   0.128  0.58  7.42           O
ATOM    465  CG BGLU A  23      10.417  11.978   0.093  0.42  1.79           C
ATOM    466  CD BGLU A  23      11.904  12.034   0.444  0.42  2.33           C
ATOM    467  OE1BGLU A  23      12.431  13.116   0.797  0.42  2.70           O
ATOM    468  OE2BGLU A  23      12.517  10.950   0.319  0.42  4.43           O
END
"""

def residue_22(hierarchy, iselection):
  """{(resname, altloc): n atoms} of residue 22 in the selection."""
  result = {}
  for a in hierarchy.atoms().select(iselection):
    ag = a.parent()
    if(ag.parent().resseq.strip() != "22"): continue
    key = (ag.resname, ag.altloc)
    result[key] = result.get(key, 0) + 1
  return result

def test_two_residue_groups():
  ph = iotbx.pdb.input(source_info=None, lines=pdb_str).construct_hierarchy()
  assert len([rg for rg in ph.residue_groups() if rg.resseq.strip()=="22"])==2
  sel = conformers.conformer_selections(ph)
  assert [s.label for s in sel] == ["A", "B", "C"]
  assert residue_22(ph, sel[0].iselection) == {("PRO", ""): 1, ("PRO", "A"): 6}
  assert residue_22(ph, sel[1].iselection) == {("PRO", ""): 1, ("SER", "B"): 5}
  assert residue_22(ph, sel[2].iselection) == {("PRO", ""): 1, ("SER", "C"): 5}
  # residue 22 has every label: only GLU 23 (A, B) falls back, in C
  assert [s.n_fallback for s in sel] == [0, 0, 1]
  # every atom of residue 22 is in exactly one conformer
  mult = conformers.multiplicity(ph)
  for a in ph.atoms():
    if(a.parent().parent().resseq.strip() == "22"):
      expected = 3 if a.parent().altloc == "" else 1
      assert mult[a.i_seq] == expected, (a.id_str(), mult[a.i_seq])

# 1pq7 CYS A 41 (A, B) and CYS A 57 (blank) with the H of reduce2: the
# disulfide SG A 41 - SG 57 takes HG A 57; HG B 57 is an H-only atom group
# on the blank SG, not a conformer of residue 57.
cys_str = """\
CRYST1   40.000   40.000   40.000  90.00  90.00  90.00 P 1
ATOM     25  N  ACYS A  41       5.508   7.360   1.752  0.59  3.97           N
ATOM     26  CA ACYS A  41       4.721   6.164   1.517  0.60  3.70           C
ATOM     27  C  ACYS A  41       3.333   6.308   2.146  0.60  3.80           C
ATOM     28  O  ACYS A  41       2.979   7.330   2.757  0.60  3.74           O
ATOM     29  CB ACYS A  41       5.431   4.898   2.010  0.61  3.50           C
ATOM     30  SG ACYS A  41       6.731   4.335   0.827  0.60  3.82           S
ATOM     31  H  ACYS A  41       5.397   7.692   2.538  0.59  3.97           H
ATOM     32  HA ACYS A  41       4.609   6.046   0.561  0.60  3.70           H
ATOM     33  HB2ACYS A  41       5.855   5.081   2.863  0.61  3.50           H
ATOM     34  HB3ACYS A  41       4.780   4.185   2.109  0.61  3.50           H
ATOM     35  N  BCYS A  41       5.374   7.503   1.660  0.40  4.85           N
ATOM     36  CA BCYS A  41       4.595   6.308   1.311  0.39  5.11           C
ATOM     37  C  BCYS A  41       3.340   6.305   2.209  0.40  4.25           C
ATOM     38  O  BCYS A  41       3.124   7.213   3.049  0.40  3.20           O
ATOM     39  CB BCYS A  41       5.483   5.140   1.450  0.40  7.12           C
ATOM     40  SG BCYS A  41       5.836   4.613   3.064  0.39 10.10           S
ATOM     41  H  BCYS A  41       5.164   7.892   2.398  0.40  4.85           H
ATOM     42  HA BCYS A  41       4.272   6.271   0.397  0.39  5.11           H
ATOM     43  HB2BCYS A  41       5.068   4.392   0.992  0.40  7.12           H
ATOM     44  HB3BCYS A  41       6.331   5.356   1.031  0.40  7.12           H
ATOM     45  HG BCYS A  41       6.022   3.427   3.056  0.39 10.10           H
ATOM     71  N   CYS A  57       9.445   0.516   1.315  1.00  4.44           N
ATOM     72  CA  CYS A  57       9.313   1.815   1.946  1.00  4.72           C
ATOM     73  C   CYS A  57      10.241   1.988   3.148  1.00  4.55           C
ATOM     74  O   CYS A  57      10.553   3.132   3.509  1.00  5.19           O
ATOM     75  CB  CYS A  57       7.875   2.020   2.439  1.00  5.19           C
ATOM     76  SG  CYS A  57       6.674   2.294   1.132  1.00  6.27           S
ATOM     77  H   CYS A  57       8.749   0.015   1.378  1.00  4.44           H
ATOM     78  HA  CYS A  57       9.546   2.481   1.280  1.00  4.72           H
ATOM     79  HB2 CYS A  57       7.600   1.229   2.929  1.00  5.19           H
ATOM     80  HB3 CYS A  57       7.856   2.795   3.022  1.00  5.19           H
ATOM     81  HG BCYS A  57       7.054   3.193   0.434  0.39  6.27           H
END
"""

def test_h_only_group():
  ph = iotbx.pdb.input(source_info=None, lines=cys_str).construct_hierarchy()
  sel = conformers.conformer_selections(ph)
  assert [s.label for s in sel] == ["A", "B"]
  assert [s.n_fallback for s in sel] == [0, 0]
  def names(iselection):
    return sorted(a.name.strip() + a.parent().altloc for a in
      ph.atoms().select(iselection) if a.parent().parent().resseq.strip()=="57")
  blank = ["C", "CA", "CB", "H", "HA", "HB2", "HB3", "N", "O", "SG"]
  assert names(sel[0].iselection) == blank, names(sel[0].iselection)
  assert names(sel[1].iselection) == sorted(blank + ["HGB"])
  mult = conformers.multiplicity(ph)
  for a in ph.atoms():
    if(a.parent().parent().resseq.strip() != "57"): continue
    assert mult[a.i_seq] == (1 if a.parent().altloc == "B" else 2), a.id_str()

def run(prefix = "qrefine_"+os.path.basename(__file__).replace(".py","")):
  os.makedirs(prefix, exist_ok=True)
  os.chdir(prefix)
  test_two_residue_groups()
  test_h_only_group()
  test_split_h_copies()

split_h_str = """\
CRYST1   34.580   30.910   26.860  90.00  90.00  90.00 P 1
ATOM     61  N  AILE A   9      50.042  12.619   9.863  0.60 15.47           N  
ATOM     62  N  BILE A   9      50.064  12.669   9.955  0.40 15.79           N  
ATOM     63  CA AILE A   9      49.760  12.880   8.443  0.60 16.16           C  
ATOM     64  CA BILE A   9      49.894  13.068   8.557  0.40 14.97           C  
ATOM     65  C  AILE A   9      50.740  12.096   7.568  0.60 16.68           C  
ATOM     66  C  BILE A   9      50.739  12.182   7.629  0.40 15.01           C  
ATOM     67  O  AILE A   9      50.344  11.527   6.549  0.60 13.30           O  
ATOM     68  O  BILE A   9      50.215  11.661   6.634  0.40 13.83           O  
ATOM     69  CB AILE A   9      49.712  14.388   8.233  0.60 17.94           C  
ATOM     70  CB BILE A   9      50.201  14.546   8.306  0.40 15.90           C  
ATOM     71  CG1AILE A   9      48.529  15.065   8.881  0.60 19.71           C  
ATOM     72  CG1BILE A   9      49.111  15.466   8.868  0.40 14.83           C  
ATOM     73  CG2AILE A   9      49.644  14.632   6.712  0.60 18.98           C  
ATOM     74  CG2BILE A   9      50.339  14.857   6.814  0.40 13.92           C  
ATOM     75  CD1AILE A   9      48.524  16.563   9.024  0.60 23.19           C  
ATOM     76  CD1BILE A   9      49.568  16.907   8.941  0.40 21.29           C  
ATOM     77  N   ALA A  10      52.006  12.040   7.976  1.00 16.22           N  
ATOM      0  H  AALA A  10      52.291  12.523   8.804  0.60 10.00           H
ATOM      0  H  BALA A  10      52.319  12.486   8.814  0.40 10.00           H
ATOM     78  CA  ALA A  10      52.970  11.271   7.207  1.00 15.86           C  
ATOM     79  C   ALA A  10      52.525   9.824   7.134  1.00 14.08           C  
ATOM     80  O   ALA A  10      52.586   9.200   6.066  1.00 14.53           O  
ATOM     81  CB  ALA A  10      54.346  11.404   7.820  1.00 18.40           C  
ATOM     82  N  AGLU A  11      52.061   9.296   8.275  0.44 15.90           N  
ATOM     83  N  BGLU A  11      52.057   9.197   8.209  0.56 15.09           N  
ATOM     84  CA AGLU A  11      51.546   7.936   8.363  0.44 14.66           C  
ATOM     85  CA BGLU A  11      51.682   7.783   8.120  0.56 15.80           C  
ATOM     86  C  AGLU A  11      50.383   7.727   7.397  0.44 13.66           C  
ATOM     87  C  BGLU A  11      50.422   7.562   7.298  0.56 14.26           C  
ATOM     88  O  AGLU A  11      50.391   6.775   6.606  0.44 12.45           O  
ATOM     89  O  BGLU A  11      50.190   6.519   6.668  0.56 17.23           O  
ATOM     90  CB AGLU A  11      51.112   7.561   9.802  0.44 17.94           C  
ATOM     91  CB BGLU A  11      51.550   7.251   9.567  0.56 16.15           C  
ATOM     92  CG AGLU A  11      52.323   7.209  10.666  0.44 21.16           C  
ATOM     93  CG BGLU A  11      52.838   7.443  10.375  0.56 20.25           C  
ATOM     94  CD AGLU A  11      52.011   7.191  12.138  0.44 25.14           C  
ATOM     95  CD BGLU A  11      52.867   6.691  11.676  0.56 22.28           C  
ATOM     96  OE1AGLU A  11      50.841   7.137  12.571  0.44 31.42           O  
ATOM     97  OE1BGLU A  11      51.880   5.993  12.030  0.56 32.07           O  
ATOM     98  OE2AGLU A  11      52.985   7.223  12.925  0.44 28.65           O  
ATOM     99  OE2BGLU A  11      53.923   6.803  12.349  0.56 26.18           O  
HETATM    0  O  CHOH A 500      63.101  10.466   8.775  0.30 20.00           O
END
"""

def test_split_h_copies():
  """1exr A 9-11 with a water labelled C: ALA 10 is blank with H A and H B
  (split by C A 9 / C B 9). Conformer C takes ILE 9 A and GLU 11 A by the
  fallback and, with them, H A of ALA 10; A and B take their own copy."""
  ph = iotbx.pdb.input(source_info=None,
    lines=split_h_str.splitlines()).construct_hierarchy()
  ph.atoms().reset_i_seq()
  atoms = ph.atoms()
  def h_copies(s):
    return sorted([atoms[i].parent().altloc for i in s.iselection
      if atoms[i].name.strip() == "H" and
         atoms[i].parent().parent().resseq.strip() == "10"])
  sel = dict([(s.label, s) for s in conformers.conformer_selections(ph)])
  assert sorted(sel.keys()) == ["A", "B", "C"], sorted(sel.keys())
  assert h_copies(sel["A"]) == ["A"], h_copies(sel["A"])
  assert h_copies(sel["B"]) == ["B"], h_copies(sel["B"])
  assert h_copies(sel["C"]) == ["A"], h_copies(sel["C"])
  # the heavy atoms of ILE 9 in C are the A ones
  for s in (sel["A"], sel["C"]):
    assert sorted(set([atoms[i].parent().altloc for i in s.iselection
      if not atoms[i].element_is_hydrogen() and
         atoms[i].parent().parent().resseq.strip() == "9"])) == ["A"]
  # the water labelled C is absent from A and B (unbonded, C only)
  assert sel["A"].n_absent == 1 and sel["C"].n_absent == 0
  m = conformers.multiplicity(ph)
  h_a = [a.i_seq for a in atoms if a.name.strip() == "H" and
         a.parent().altloc == "A"][0]
  h_b = [a.i_seq for a in atoms if a.name.strip() == "H" and
         a.parent().altloc == "B"][0]
  assert m[h_a] == 2 and m[h_b] == 1, (m[h_a], m[h_b])

if(__name__ == "__main__"):
  run()
  print("OK")
