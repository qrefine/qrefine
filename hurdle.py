# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function
#
import os
import iotbx
from libtbx.utils import Sorry
from libtbx import Auto
from phenix.program_template import ProgramTemplate

# =============================================================================
class Program(ProgramTemplate):

  description = '''
Program for testing models for Q|R

Minimum required data:
  Model file

'''

  datatypes = ['model', 'phil']

  master_phil_str = '''
inputs
{
  qm_method = xtb
    .type = str
  max_atoms = 5000
    .type = int
}
'''
  results={}

  def validate(self):
    assert self.params.inputs.qm_method in ['xtb']

  def run(self):
    print('Using model: %s' % self.data_manager.get_default_model_name())
    prefix = '%s' % os.path.splitext(self.data_manager.get_model_names()[0])[0]
    #
    model = self.data_manager.get_model()
    self.results['has_hd']=model.has_hd()
    if not self.results['has_hd']:
      print('\n\tNeed to add Hydrogen atoms\n')
    #
    ph = model.get_hierarchy()
    # print(dir(ph))
    self.results['max_atoms']=len(ph.atoms())
    self.results['altloc_indices']=list(ph.altloc_indices())
    del self.results['altloc_indices'][0]
    self.results['is_ca_only'] = ph.is_ca_only()

    sd = model.scattering_dictionary()
    self.results['elements']=list(sd.keys())

    xrs = model.get_xray_structure()
    self.results['special_position_indices']=list(xrs.special_position_indices())

    print(self.results)
    must_be_false = ['is_ca_only',
                     'special_position_indices',
                     'altloc_indices',
                     ]
    for attr in must_be_false:
      if self.results[attr]:
        raise Sorry('\nModel has a True value for "%s" of "%s"' % (attr,
                                                                 self.results[attr]))

    max_atoms = self.params.inputs.max_atoms
    if not self.results['has_hd']:
      max_atoms*=2
    if self.results['max_atoms']>max_atoms:
      raise Sorry('\nToo many atoms : "%d"' % (self.results['max_atoms']))

    from elbow.chemistry.AtomClass import AtomClass
    for element in self.results['elements']:
      a = AtomClass(element[:2])
      if a.isMetal():
        raise Sorry('\nModel contains metal element "%s"' % element)

    return self.results

  def get_results(self):
    return self.results
