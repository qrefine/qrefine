from libtbx.utils import Sorry
import numpy as np
try:
  from jpype import shutdownJVM
  from jpype import JClass
  from jpype import JString
  from jpype import getDefaultJVMPath
  from jpype import java
  from jpype import startJVM
  import jpype
except ImportError, e:
  print str(e)
  jpype = None
#  raise Sorry(str(e))

class PYoink(object):

  def __init__(self, yoink_jar_path, yoink_dat_path=None, input_file=None,out_file=None,system=None):
    if jpype is None:
      raise Sorry('module jpype not loaded. Failed to load with the following:\n\n%s\n' % e)
    if not jpype.isJVMStarted():
      startJVM(getDefaultJVMPath(), "-Djava.class.path="+yoink_jar_path)
    Yoink=JClass("org.wallerlab.yoink.Yoink")
    self.yoink_dat_path = yoink_dat_path
    yoink=Yoink()
    JavaApplicationContext=JClass(
      "org.springframework.context.annotation.AnnotationConfigApplicationContext")
    javaApplicationContext=JavaApplicationContext()
    yoink.getBeans(javaApplicationContext)
    FileAdaptiveQMMMProcessor=JClass(
      "org.wallerlab.yoink.service.processor.FileYoinkProcessor")
    self.adaptiveQMMM=javaApplicationContext.getBean(FileAdaptiveQMMMProcessor)
    JaxbFileReader=JClass("org.wallerlab.yoink.molecule.data.JaxbFileReader")
    self.jaxbFileReader=javaApplicationContext.getBean(JaxbFileReader)
    self.Cml=JClass("org.xml_cml.schema.Cml")
    JaxbFileWriter=JClass("org.wallerlab.yoink.molecule.data.JaxbFileWriter")
    self.jaxbFileWriter=javaApplicationContext.getBean(JaxbFileWriter)
    self.result=None
    self.input_file=input_file
    self.out_file=out_file
    self.system=system
    self.set_up()

  def set_up(self):
    if(self.input_file is not None):
      self.jaxb_cml=self.jaxbFileReader.read(JString(self.input_file), self.Cml())
      self.atoms,self.molecules=self.get_atoms_molecules()

  def shutdown(self):
    shutdownJVM()

  def get_interactions_list(self):
    #print self.input_file
    result=self.adaptiveQMMM.process(self.input_file)
    interactions_temp=result.getGraph().getEdges()
    interactions_list=[]
    weights_temp=result.getGraph().getWeights()
    weights=[]
    for i in  range (interactions_temp.size()):
      temp=interactions_temp.get(i)
      plist=[temp.get(0).intValue(),temp.get(1).intValue()]
      interactions_list.append(plist)
      temp=weights_temp.get(i).doubleValue()
      weights.append(temp)
    return interactions_list,weights

  def get_qm_indices(self):
    self.result=self.adaptiveQMMM.process(self.input_file)
    #Region=JClass("org.wallerlab.yoink.api.model.regionizer.Region")
    RegionName=JClass("org.wallerlab.yoink.api.model.region.Region$Name")
    qm_atoms=self.result.getRegions().get(RegionName.valueOf("QM")).getAtoms()
    qm_atom_indices=[]
    qm_size= qm_atoms.size()
    for i in range(qm_size):
      atom=qm_atoms.get(i)
      qm_atom_indices.append(atom.getIndex())
    qm_atom_indices=np.sort(qm_atom_indices)
    qm_molecules_temp=self.result.getRegions().get(
      RegionName.valueOf("QM")).getMolecules()
    qm_molecules=java.util.ArrayList()
    qm_molecules.addAll(qm_molecules_temp)
    qm_molecule_indices=[]
    qm_size= qm_molecules.size()
    for i in range(qm_size):
      molecule=qm_molecules.get(i)
      qm_molecule_indices.append(molecule.getIndex())
    qm_molecule_indices=np.sort(qm_molecule_indices)
    return qm_atom_indices, qm_molecule_indices

  def write_result(self):
    if(self.out_file==None):
      self.out_file=self.input_file
    print "self.out_file",self.out_file
    self.jaxbFileWriter.write(JString(self.out_file), \
                              self.result.getInput().getValue())

  def update_input_file(self, positions=None,qm_core_fixed=None):
    Double=JClass("java.lang.Double")
    if(positions is not  None):
      for iposition, position in enumerate(positions):
        atom=self.atoms[iposition]
        atom.setX3(Double(float(position[0])))
        atom.setY3(Double(float(position[1])))
        atom.setZ3(Double(float(position[2])))
    if(qm_core_fixed != None):
      for m in self.molecules:
        m.setId(JString("MM"))
      for q in qm_core_fixed:
        self.molecules[q-1].setId(JString("QM_CORE_FIXED"))
    self.jaxbFileWriter.write(JString(self.input_file),self.jaxb_cml.getValue())

  def get_atoms_molecules(self):
    MoleculeList=JClass("org.xml_cml.schema.MoleculeList")
    Molecule=JClass("org.xml_cml.schema.Molecule")
    AtomArray=JClass("org.xml_cml.schema.AtomArray")
    Atom=JClass("org.xml_cml.schema.Atom")
    cmlSystem=self.jaxb_cml.getValue().getAnyCmlOrAnyOrAny()
    atoms=[]
    molecules=[]
    self.qm_core_fixed_indices=[]
    counter=0
    for l in range(cmlSystem.size()):
      if(cmlSystem.get(l).getValue().getClass()==MoleculeList):
        moleculeList=cmlSystem.get(l).getValue().getAnyCmlOrAnyOrAny()
        for i in range(moleculeList.size()):
          if(moleculeList.get(i).getValue().getClass()==Molecule):
            molecules.append(moleculeList.get(i).getValue())
            molecule_id=str( moleculeList.get(i).getValue().getId())
            qm_core=( str( moleculeList.get(i).getValue().getId())=="QM_CORE_FIXED")
            molecule=moleculeList.get(i).getValue().getAnyCmlOrAnyOrAny()
            for j in range(molecule.size()):
              if(molecule.get(j).getValue().getClass()==AtomArray):
                atomArray=molecule.get(j).getValue().getAnyCmlOrAnyOrAny()
                for k in range(atomArray.size()):
                  if(atomArray.get(k).getValue().getClass()==Atom):
                    atom=atomArray.get(k).getValue()
                    atoms.append(atom)
                    counter=counter+1
                    if(qm_core):
                      self.qm_core_fixed_indices.append(counter)
    return atoms,molecules

  def update(self,qm_core_fixed=None,positions=None):
    self.set_up()
    if(self.system!=None):
      positions=self.system.get_positions()
    if(positions !=None  and  qm_core_fixed!=None):
      self.update_input_file(positions,qm_core_fixed)
    elif(qm_core_fixed!=None and positions ==None):
      self.update_input_file(qm_core_fixed=qm_core_fixed)
    elif(qm_core_fixed ==None and positions!=None):
      self.update_input_file(positions)
    self.set_up()

  def write_clustering_qmmm_xmls(self, clustering_file_name, positions, qmmm_file_name, symbols):
    yoink_dat_path = self.yoink_dat_path
    cluster_file = open(clustering_file_name, "w+")
    cluster_file.write(
      """<?xml version="1.0" encoding="UTF-8" standalone="yes"?> \n """)
    cluster_file.write("""<cml xmlns="http://www.xml-cml.org/schema"> \n """)
    cluster_file.write("<moleculeList> \n")
    for resid in range(len(positions)):
      cluster_file.write("""<molecule id="MM"> \n""")
      cluster_file.write("<atomArray> \n")
      for j in range(len(positions[resid])):
        atom = "<atom id=\"" + str(j) + "\" " + " elementType=\"" \
          + symbols[resid][j] + "\" " + " x3=\"" + str(positions[resid][j][0]) \
          + "\"  " + "y3=\"" + str(positions[resid][j][1]) + "\"" + " z3=\"" \
          + str(positions[resid][j][2]) + "\"" + "></atom> \n"
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
    cluster_file.write(
      """  <parameter name="WFC_PATH" value="%(yoink_dat_path)s" /> \n """ % locals())
    cluster_file.write("""  <parameter name="REGION_CUBE" value="SYSTEM"/> \n """)
    cluster_file.write("""  <parameter name="PARTITIONER" value="INTERACTION"/>  \n""")
    cluster_file.write("""  </parameterList> \n""")
    cluster_file.write(""" </cml> \n""")
    cluster_file.close()
    lines = open(clustering_file_name, "r").readlines()
    qmmm_file = open(qmmm_file_name, "w+")
    lines[-4] = """  <parameter name="REGION_CUBE" value="ADAPTIVE_SEARCH"/> \n """
    lines[-3] = """   <parameter name="PARTITIONER" value="DORI"/>  \n"""
    for line in lines:
      qmmm_file.write(line)
    qmmm_file.close()
