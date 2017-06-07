import os
import time

#get pure terachem gradients for comparison
pdbdir = "/home/xuyanting/qr-tests-p1/01_finalised"
coords="/home/xuyanting/qr-tests-p1/01_finalised/"
work_dir ="/home/xuyanting/sp/"
pdbs =[]
for pdb_file in os.listdir(pdbdir):
 # pdbs.append(os.path.join(pdbdir,pdb_file))
  pdbs.append(pdb_file[:-4])

template = """run gradient
$multibasis
Se lanl2dz_ecp
Cl lanl2dz_ecp
Cd lanl2dz_ecp
Zn lanl2dz_ecp
Mg lanl2dz_ecp
$end
basis sto-3g
scf diis+a
coordinates """

template2 ="""
gpumem 512
charge 6
seed 1351351
maxit 200
threall 1e-12
pointcharges no
gpus 4
watcheindiis yes
method rhf
dftd yes
end
 """

def pbs():
  for pdb in pdbs:
     with open(work_dir+pdb+".sp","w") as f:
        f.write(template + coords + pdb + ".pdb " +template2)
        
pbs()
