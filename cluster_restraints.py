from __future__ import division
from __future__ import print_function
from __future__ import absolute_import

import os
from libtbx import Auto
from libtbx.utils import Sorry
from libtbx import adopt_init_args
from libtbx.easy_mp import parallel_map
from scitbx.array_family import flex
from .restraints import from_qm
from qrefine.utils.hierarchy_utils import write_hierarchy
from libtbx import group_args

class from_cluster(object):
  """
  Restraints from fragments. keep_error_copy: when the fragment jobs fail,
  copy the working folder (that of restraints_manager.file_name) into its
  sub-folder error (replacing an earlier one) before the Sorry.
  """
  def __init__(self, restraints_manager, fragment_manager, parallel_params,
               keep_error_copy=False):
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
    file_name = getattr(self.restraints_manager, "file_name", None)
    if(file_name is not None and isinstance(self.restraints_manager, from_qm)):
      # the current super sphere, for inspection next to the fragments
      write_hierarchy(self.fragment_manager.pdb_hierarchy_super, file_name,
        crystal_symmetry=self.fragment_manager.expansion.cs_box)
    # The cctbx restraints of the super sphere are needed by from_cctbx only;
    # from_qm must not carry them (pickled with every parallel_map job).
    fragment_extracts_obj = self.fragment_manager.get_fragment_extracts(
      with_super_sphere_grm = not isinstance(self.restraints_manager, from_qm))
    self.restraints_manager.fragment_extracts = fragment_extracts_obj
    selection_and_sites_cart=[]
    for index, selection_fragment in enumerate(
                       self.fragment_manager.fragment_selections):
      selection_and_sites_cart.append([selection_fragment, sites_cart,index])
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
        import sys, traceback
        import shutil
        # Keep a copy of the working folder (<folder>/error) for inspection
        error_folder = None
        if(self.keep_error_copy and file_name is not None):
          folder = os.path.dirname(file_name) or "."
          error_folder = os.path.join(folder, "error")
          if os.path.exists(error_folder):
            shutil.rmtree(error_folder)
          try:
            shutil.copytree(folder, error_folder,
              ignore=shutil.ignore_patterns("error"))
          except (shutil.Error, OSError) as copy_error:
            print('Directory not copied. Error: %s' % copy_error)
            error_folder = None
        ncount=ncount+1
        energy_gradients=None
        print("check independent QM jobs")
        # Sometimes the 'standard' traceback is not available.
        # Below sort of forces the same information at the risk of doing things twice.
        # It was needed to find some bugs, but perhaps needs to revisisted later.
        exc_type, exc_value, exc_traceback = sys.exc_info()
        traceback_template = ''' ** qrefine exception handler: **
        %(type)s => File "%(filename)s" \n line %(lineno)s, in %(name)s: \n %(message)s
         \n'''
        traceback_details = {
                         'filename': exc_traceback.tb_frame.f_code.co_filename,
                         'lineno'  : exc_traceback.tb_lineno,
                         'name'    : exc_traceback.tb_frame.f_code.co_name,
                         'type'    : exc_type.__name__,
                         'message' : exc_value.args[0]
                        }
        del(exc_type, exc_value, exc_traceback)
        print(traceback.format_exc())
        print(traceback_template % traceback_details)
        msg = 'process finished with error: %s' % e
        if(error_folder is not None):
          msg += ' (working folder copied to %s)' % error_folder
        raise Sorry(msg)
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
    return target, gradients
