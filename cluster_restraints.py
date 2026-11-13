from __future__ import division
from __future__ import print_function
from __future__ import absolute_import

import os
from libtbx import Auto
from libtbx.utils import Sorry
from libtbx import adopt_init_args
from libtbx.easy_mp import parallel_map
from scitbx.array_family import flex
from .restraints import from_qm, from_cctbx
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
    if(isinstance(self.restraints_manager, from_cctbx)):
      # cctbx restraints of the fragments: in this process. A process per
      # fragment costs far more than the fragment's restraints (1PWM, 373
      # jobs: about 19 s per evaluation with nproc=4, 0.2 s here), and the
      # fragment restraints selected at the first evaluation are kept for
      # the next ones only in this process (from_cctbx._fragment_grm).
      # An exception propagates as it is (no Sorry wrapper, no keep_error_copy
      # folder) and parallel.nproc is not used for cctbx fragments.
      energy_gradients = [self.restraints_manager(item)
        for item in selection_and_sites_cart]
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
    return self._combine(energy_gradients, system_size)

  def _combine(self, energy_gradients, system_size):
    """
    Weighted combination of the job (conformer variant) results: per atom,
    num = sum of w_j*g_j and den = sum of w_j over the jobs having it as a
    cluster atom (buffer atoms are zeroed per job), gradient num/den; den is
    the conformer multiplicity of the atom (1 without conformers, where
    num/den is the plain sum). Target: sum_j w_j*e_j / n_conformers, with
    n_conformers the number of conformer labels (the weight sum of every
    fragment, less for a fragment whose cluster is absent from some
    conformers: fragment._variants): the mean over conformers of the sum of
    fragment energies, not a model energy (gradient_only is forced with
    clustering).
    Sets self.n_jobs and self.multiplicity (den), available after the first
    call.
    """
    fm = self.fragment_manager
    weights = fm.variant_weights
    self.n_jobs = len(weights)
    assert len(energy_gradients) == self.n_jobs
    target = 0
    num = flex.vec3_double(system_size)
    den = flex.double(system_size, 0)
    for index, item in enumerate(energy_gradients):
      t, g = item[0], item[1]
      w = weights[index]
      selection_fragment = fm.fragment_selections[index]
      selection_buffer = fm.buffer_selections[index]
      gradients_i = flex.vec3_double(system_size)
      gradients_i = gradients_i.set_selected(selection_fragment, g)
      gradients_i = gradients_i.set_selected(selection_buffer, [0,0,0])
      if(w != 1): gradients_i = gradients_i * float(w)
      num += gradients_i
      core = selection_fragment & ~selection_buffer
      den.set_selected(core, den.select(core) + w)
      target += w * t
    cluster = flex.bool(system_size, False)
    for s in fm.fragment_cluster_selections:
      cluster = cluster | s
    assert (den.select(cluster) > 0).all_eq(True)
    den.set_selected(den == 0, 1)
    self.multiplicity = den
    weight_sums = [0] * fm.n_fragments
    for index, i in enumerate(fm.variant_fragment_index):
      weight_sums[i] += weights[index]
    labels = getattr(fm, "labels", [])
    if(len(labels) == 0): assert len(set(weight_sums)) <= 1, weight_sums
    n_conformers = max([1, len(labels)] + weight_sums)
    # division by 1 is exact: the single-conformer sum is unchanged
    return target / n_conformers, num / den
