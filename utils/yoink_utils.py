from __future__ import division

def write_yoink_infiles(cluster_file_name,qmmm_file_name,pdb_hierarchy,dat_dir):

  def read_pdb_hierarchy(pdb_hierarchy):
    positions = []
    symbols = []
    for chain in pdb_hierarchy.chains():
      for residue_group in chain.residue_groups():
        symbols_res = []
        positions_res = []
        for atom in residue_group.atoms():
          element = atom.element.strip()
          if(len(element) == 2):
            element = element[0]+ element[1].lower()
          symbols_res.append(element)
          positions_res.append(list(atom.xyz))
        positions.append(positions_res)
        symbols.append(symbols_res)
    return symbols,positions

  def write_xml(cluster_file_name,qmmm_file_name,pdb_hierarchy,dat_dir):
    symbols,positions = read_pdb_hierarchy(pdb_hierarchy)
    cluster_file = open(cluster_file_name,"w+")
    cluster_file.write("""<?xml version="1.0" encoding="UTF-8" standalone="yes"?> \n """)
    cluster_file.write("""<cml xmlns="http://www.xml-cml.org/schema"> \n """)
    cluster_file.write("<moleculeList> \n")
    for resid in range(len(positions)):
      cluster_file.write("""<molecule id="MM"> \n""")
      cluster_file.write("<atomArray> \n")
      for j in range(len(positions[resid])):
        atom = "<atom id=\""+str(j)+"\" "+" elementType=\""+symbols[resid][j]+"\" "+" x3=\""+ str(positions[resid][j][0])+"\"  "+ "y3=\""+str(positions[resid][j][1])+"\""+" z3=\""+str(positions[resid][j][2])+"\""+"></atom> \n"
        cluster_file.write(atom)
      cluster_file.write("</atomArray> \n")
      cluster_file.write("</molecule> \n")
    cluster_file.write("</moleculeList> \n")
    cluster_file.write("""<parameterList title="parameters"> \n""")
    cluster_file.write(""" <parameter name="DENSITY_BUFFER" value="1"/> \n""")
    cluster_file.write("""  <parameter name="DENSITY_QM" value="0.001"/> \n""")
    cluster_file.write("""  <parameter name="DENSITY_ASR_QM" value="0.000000001"/>\n""")
    cluster_file.write("""  <parameter name="DENSITY_ASR_QMCORE" value="0.1"/>  \n """)
    cluster_file.write("""  <parameter name="DENSITY_DORI" value="0.001"/>  \n """)
    cluster_file.write("""  <parameter name="DENSITY_SEDD" value="0.1"/>   \n""")
    cluster_file.write("""  <parameter name="SEDD" value="5"/> \n """)
    cluster_file.write("""  <parameter name="DORI" value="0.8"/>  \n""")
    cluster_file.write("""  <parameter name="DENSITY_RATIO_MIN" value="0.06"/> \n """)
    cluster_file.write("""  <parameter name="DENSITY_RATIO_MAX" value="16"/>  \n""")
    cluster_file.write("""  <parameter name="SEDD_STEPSIZE" value="0.3 0.3 0.3"/>  \n """)
    cluster_file.write("""  <parameter name="DORI_STEPSIZE" value="0.5 0.5 0.5"/>  \n """)
    cluster_file.write("""  <parameter name="SMOOTHNER" value="ABRUPT"/>  \n """)
    cluster_file.write("""  <parameter name="DGRID" value="true"/>  \n""")
    cluster_file.write("""  <parameter name="INTERACTION_WEIGHT" value="false"/> \n """)
    cluster_file.write("""  <parameter name="WFC_PATH" value="%(dat_dir)s" /> \n """%locals())
    cluster_file.write("""  <parameter name="REGION_CUBE" value="SYSTEM"/> \n """)
    cluster_file.write("""  <parameter name="PARTITIONER" value="INTERACTION"/>  \n""")
    cluster_file.write("""  </parameterList> \n""")
    cluster_file.write(""" </cml> \n""")
    cluster_file.close()
    lines = open(cluster_file_name,"r").readlines()
    qmmm_file = open(qmmm_file_name,"w+")
    lines[-4] = """  <parameter name="REGION_CUBE" value="ADAPTIVE_SEARCH"/> \n """
    lines[-3] = """   <parameter name="PARTITIONER" value="DORI"/>  \n"""
    for line in lines:
      qmmm_file.write(line)
    qmmm_file.close()

  write_xml(cluster_file_name,qmmm_file_name,pdb_hierarchy,dat_dir)
