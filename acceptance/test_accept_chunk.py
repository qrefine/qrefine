import os
import libtbx
from libtbx import easy_run
from libtbx.command_line import easy_qsub

# swith to either ENV variable or command line argument
rootdir = "/home/pdb/mirror/pub/pdb/data/structures/divided/pdb"
qsub_str= "qsub  -N test_clusters -M xuyanting@i.shu.edu.cn -m ae -q qr  -l nodes=1:ppn=4"
work_dir="/home/xuyanting/work"
phenix_source="/home/xuyanting/phenix/phenix-1.11.1-2575/build/setpaths_all.sh"

commands =[]
for subdir, dirs, files in os.walk(rootdir):
    for file in files:
      cmd = " phenix.pdb.hierarchy %s >> %s.log" % (os.path.join(subdir, file),file)
      #easy_run.call(cmd)
      commands.append(cmd)

easy_qsub.run(
    phenix_source=phenix_source,
    where=work_dir,
    size_of_chunks = 50,
    qsub_cmd=qsub_str,
    commands= commands )
