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
  pass
#  raise Sorry(str(e))

class PYoink(object):

  def __init__(self, yoink_jar_path, input_file=None,out_file=None,system=None):
    if not jpype.isJVMStarted():
      startJVM(getDefaultJVMPath(), "-Djava.class.path="+yoink_jar_path)	
    Yoink=JClass("org.wallerlab.yoink.Yoink")
    yoink=Yoink()
    JavaApplicationContext=JClass("org.springframework.context.annotation.AnnotationConfigApplicationContext")
    javaApplicationContext=JavaApplicationContext()
    yoink.getBeans(javaApplicationContext)
    FileAdaptiveQMMMProcessor=JClass("org.wallerlab.yoink.service.processor.FileAdaptiveQMMMProcessor")
    self.adaptiveQMMM=javaApplicationContext.getBean(FileAdaptiveQMMMProcessor)
    JaxbFileReader=JClass("org.wallerlab.yoink.molecular.data.JaxbFileReader")
    self.jaxbFileReader=javaApplicationContext.getBean(JaxbFileReader)
    self.Cml=JClass("org.xml_cml.schema.Cml")
    JaxbFileWriter=JClass("org.wallerlab.yoink.molecular.data.JaxbFileWriter")
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

  def get_interaction_list(self):
    result=self.adaptiveQMMM.process(self.input_file)
    interaction_temp= result.getInteractionList()
    interaction_list=[]
    weight_temp=result.getInteractionWeight()
    weight=[]
    for i in  range (interaction_temp.size()):
      temp=interaction_temp.get(i)
      plist=[temp.get(0).intValue(),temp.get(1).intValue()]
      interaction_list.append(plist)
      temp=weight_temp.get(i).doubleValue()
      weight.append(temp)
    return interaction_list,weight


  def get_qm_indices(self):
    self.result=self.adaptiveQMMM.process(self.input_file)
    #Region=JClass("org.wallerlab.yoink.api.model.regionizer.Region")
    RegionName=JClass("org.wallerlab.yoink.api.model.regionizer.Region$Name")
    qm_atoms=self.result.getRegions().get(RegionName.valueOf("QM")).getAtoms()
    qm_atom_indices=[]
    qm_size= qm_atoms.size()
    for i in range(qm_size):
      atom=qm_atoms.get(i)
      qm_atom_indices.append(atom.getIndex())
    qm_atom_indices=np.sort(qm_atom_indices)
    qm_molecules_temp=self.result.getRegions().get(RegionName.valueOf("QM")).getMolecules()
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
    self.jaxbFileWriter.write(JString(self.out_file),self.result.getInput().getValue())

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
