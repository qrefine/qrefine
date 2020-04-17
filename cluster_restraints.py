from __future__ import division
from __future__ import print_function
from __future__ import absolute_import

import os
from libtbx import Auto
from libtbx.utils import Sorry
from libtbx import adopt_init_args
from libtbx.easy_mp import parallel_map
from scitbx.array_family import flex
from .fragment import fragment_extracts, write_cluster_and_fragments_pdbs
from .restraints import from_qm
from libtbx import group_args

def check_no_altlocs(h, file_name):
  altlocs = []
  for m in h.models():
    for c in m.chains():
      for rg in c.residue_groups():
        for ag in rg.atom_groups():
          altlocs.append(ag.altloc)
  altlocs = list(set(altlocs))
  assert len(altlocs)<=2, [file_name, altlocs]
  cntr = 0
  for a in altlocs:
    if(len(a)==0): cntr+=1
  assert cntr==1, [file_name, altlocs]

class from_cluster(object):
  def __init__(self, restraints_manager, fragment_manager, parallel_params):
    adopt_init_args(self, locals())

  def energies_sites(self, sites_cart, compute_gradients=True):
    tg = self.target_and_gradients(sites_cart=sites_cart)
    return group_args(
      target    = tg[0],
      gradients = tg[1])

  def target_and_gradients(self, sites_cart):
    # update the pdb hierarchy of the center
    system_size = sites_cart.size()
    self.fragment_manager.update_xyz(sites_cart)
    sites_cart = self.fragment_manager.pdb_hierarchy_super.atoms().extract_xyz()
    #
    self.fragment_manager.pdb_hierarchy_super.write_pdb_file(
      file_name=self.restraints_manager.file_name,
      crystal_symmetry=self.fragment_manager.expansion.cs_box)
    fragment_extracts_obj = fragment_extracts(self.fragment_manager)
    # super_sphere_geometry_restraints_manager will cause qusb submits
    # a single job instead of batch jobs
    if(isinstance(self.restraints_manager, from_qm)):
      fragment_extracts_obj.super_sphere_geometry_restraints_manager=None
    self.restraints_manager.fragment_extracts = fragment_extracts_obj
    selection_and_sites_cart=[]
    write_cluster_and_fragments_pdbs(fragments=fragment_extracts_obj,directory='frag_pdbs')
    for index, selection_fragment in enumerate(
                       self.fragment_manager.fragment_selections):
      selection_and_sites_cart.append([selection_fragment, sites_cart,index])
       ## DEBUG begin
       #try:
       #  super_selection = self.restraints_manager.\
       #    fragment_extracts.fragment_super_selections[index]
       #  tmp_h = self.fragment_manager.pdb_hierarchy_super.select(super_selection)
       #  fn = self.restraints_manager.file_name.replace(".pdb","_%s.pdb"%str(index))
       #  tmp_h.write_pdb_file(
       #    file_name        = fn,
       #    crystal_symmetry = self.fragment_manager.expansion.cs_box)
       #  check_no_altlocs(h=tmp_h, file_name=fn)
       #  self.restraints_manager.target_and_gradients(sites_cart=sites_cart,
       #    selection=selection_fragment, index=index)
       #except Exception as e:
       #  print "Failed before entering lbfgs minimization"
       #  print str(e)
       #  STOP()
       ## DEBUG end
    if(self.parallel_params.nproc is None):
      self.parallel_params.nproc = Auto
    ncount=0
    energy_gradients=None
    while(ncount<5 and energy_gradients is None):
      try:
        energy_gradients = parallel_map(
          func                       = self.restraints_manager,
          iterable                   = selection_and_sites_cart,
          method                     = self.parallel_params.method,
          preserve_exception_message = True,
          processes                  = self.parallel_params.nproc,
          qsub_command               = self.parallel_params.qsub_command,
          use_manager                = True)
      except Exception as e:
        import shutil
        if os.path.exists('ase_error'):
          shutil.rmtree('ase_error')
        try:
          shutil.copytree('ase', 'ase_error')
        # Directories are the same
        except shutil.Error as e:
          print('Directory not copied. Error: %s' % e)
        # Any error saying that the directory doesn't exist
        except OSError as e:
          print('Directory not copied. Error: %s' % e)
        ncount=ncount+1
        energy_gradients=None
        print("check independent QM jobs")
        print(e)
        raise Sorry('process finished with error: %s' % e)
    target=0
    gradients=flex.vec3_double(system_size)
    for index, item in enumerate(energy_gradients):
      t = item[0]
      g = item[1]
      selection_fragment = self.fragment_manager.fragment_selections[index]
      selection_buffer =  self.fragment_manager.buffer_selections[index]
      gradients_i = flex.vec3_double(system_size)
      gradients_i = gradients_i.set_selected(selection_fragment, g)
      gradients_i = gradients_i.set_selected(selection_buffer,[0,0,0])
      gradients += gradients_i
      target += t
    ### for debugging parallel_map, remove later.
    if(0):
      for selection_fragment, selection_buffer in zip(
                                      self.fragment_manager.fragment_selections,
                                      self.fragment_manager.buffer_selections):
        t, g =  self.restraints_manager.target_and_gradients(
                           sites_cart=sites_cart, selection=selection_fragment)
        # XXX INEFFICIENT, MOVE TO C++
        gradients_i = flex.vec3_double(system_size)
        gradients_i = gradients_i.set_selected(selection_fragment, g)
        gradients_i = gradients_i.set_selected(selection_buffer,[0,0,0])
        gradients += gradients_i
        target += t
    return target, gradients
