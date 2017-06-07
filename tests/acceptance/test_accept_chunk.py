import os
import libtbx
from libtbx import easy_run
from libtbx.command_line import easy_qsub

# swith to either ENV variable or command line argument
rootdir = "/home/pdb/mirror/pub/pdb/data/structures/divided/pdb"
qsub_str= "qsub  -N test_clusters -M xuyanting@i.shu.edu.cn -m ae -q qr  -l nodes=1:ppn=4"
work_dir="/home/xuyanting/work"
phenix_source="/home/xuyanting/phenix/phenix-1.11.1-2575/build/setpaths_all.sh"

def acceptance_tests():
  return  {
    "hierarchy"      ："phenix.pdb.hierarchy"
    " finalise"      : "qr.final" ,
    " chunk"         : "qr.chunk",
    " restraints"    : "qr.restrain",
    " refinement"    : "qr.refine",
    }

commands =[]
for subdir, dirs, files in os.walk(rootdir):
    for file in files:
      for name , command in acceptance_tests().items():
        cmd = "%s  %s > %s.log" % (os.path.join(command,subdir, file),file)
      #easy_run.call(cmd)
        commands.append(cmd)


#TODO compare with  /home/xuyanting/test_prime/elbow.py
easy_qsub.run(
    phenix_source=phenix_source,
    where=work_dir,
    size_of_chunks = 50,
    qsub_cmd=qsub_str,
    commands= commands )
