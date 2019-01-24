import iotbx.pdb
from scitbx.array_family import flex
from cctbx import uctbx
from cctbx import crystal
from libtbx.utils import null_out
import mmtbx
import mmtbx.monomer_library.server
import mmtbx.monomer_library.pdb_interpretation
import mmtbx.restraints
from mmtbx import monomer_library
from libtbx import group_args
import math
import iotbx.pdb.utils

mon_lib_srv = mmtbx.monomer_library.server.server()
ener_lib    = mmtbx.monomer_library.server.ener_lib(
  use_neutron_distances = True)

def create_super_sphere(pdb_hierarchy,
                        crystal_symmetry,
                        select_within_radius,
                        link_min=1.0,
                        link_max=2.0,
                        r=None):
  if(r is None):
    # This is equivalent to (but much faster):
    #
    #r = grm.geometry.pair_proxies().nonbonded_proxies.\
    #  get_symmetry_interacting_indices_unique(
    #    sites_cart = pdb_hierarchy.atoms().extract_xyz())
    #
    sites_cart = pdb_hierarchy.atoms().extract_xyz()
    sst = crystal_symmetry.special_position_settings().site_symmetry_table(
      sites_cart = sites_cart)
    r = {}
    from cctbx import crystal
    cutoff = select_within_radius+1 # +1 is nonbonded buffer
    conn_asu_mappings = crystal_symmetry.special_position_settings().\
      asu_mappings(buffer_thickness=cutoff)
    conn_asu_mappings.process_sites_cart(
      original_sites      = sites_cart,
      site_symmetry_table = sst)
    conn_pair_asu_table = crystal.pair_asu_table(
        asu_mappings=conn_asu_mappings)
    conn_pair_asu_table.add_all_pairs(
        distance_cutoff=cutoff)
    pair_generator = crystal.neighbors_fast_pair_generator(
       conn_asu_mappings,
       distance_cutoff=cutoff)
    for pair in pair_generator:
      rt_mx_i = conn_asu_mappings.get_rt_mx_i(pair)
      rt_mx_j = conn_asu_mappings.get_rt_mx_j(pair)
      rt_mx_ji = rt_mx_i.inverse().multiply(rt_mx_j)
      if str(rt_mx_ji)=="x,y,z": continue
      r.setdefault(pair.j_seq, []).append(rt_mx_ji)
    for k,v in zip(r.keys(), r.values()): # remove duplicates!
      r[k] = list(set(v))
  c = iotbx.pdb.hierarchy.chain(id = "SS") # all symmetry related full residues
  fm = crystal_symmetry.unit_cell().fractionalization_matrix()
  om = crystal_symmetry.unit_cell().orthogonalization_matrix()
  selection = r.keys()
  for rg in pdb_hierarchy.residue_groups():
    keep=False
    for i in rg.atoms().extract_i_seq():
      if(i in selection):
        keep=True
        break
    if(keep):
      ops = r[i]
      for op in ops:
        rg_ = rg.detached_copy()
        xyz = rg_.atoms().extract_xyz()
        new_xyz = flex.vec3_double()
        for xyz_ in xyz:
          t1 = fm*flex.vec3_double([xyz_])
          t2 = op*t1[0]
          t3 = om*flex.vec3_double([t2])
          new_xyz.append(t3[0])
        rg_.atoms().set_xyz(new_xyz)
        rg_.link_to_previous=True
        c.append_residue_group(rg_)
  resseq = 0
  for rg in c.residue_groups():
    rg.resseq = "%4d" % resseq
    resseq+=1
  # Now we need to re-order residues and group by chains so that they can be
  # linked by pdb interpretation.
  super_sphere_hierarchy = pdb_hierarchy.deep_copy()
  #
  all_chains = iotbx.pdb.utils.all_chain_ids()
  all_chains = [chid.strip() for chid in all_chains]
  for cid in list(set([i.id for i in pdb_hierarchy.chains()])):
    all_chains.remove(cid)
  all_chains = iter(all_chains)
  #
  def get_atom(atoms, name):
    for a in atoms:
      if(a.name.strip().lower()==name.strip().lower()):
        return a.xyz
  residue_groups_i = list(c.residue_groups())
  residue_groups_j = list(c.residue_groups())
  residue_groups_k = list() # standard aa residues only
  #
  def dist(r1,r2):
    return math.sqrt((r1[0]-r2[0])**2+(r1[1]-r2[1])**2+(r1[2]-r2[2])**2)
  #
  def grow_chain(residue_groups_j, chunk, ci):
    n_linked = 0
    for rgj in residue_groups_j:
      nj = get_atom(atoms=rgj.atoms(), name="N")
      if(nj is None): continue
      d_ci_nj = dist(ci,nj)
      if(d_ci_nj<link_max and d_ci_nj>link_min):
        n_linked += 1
        chunk.append(rgj)
        residue_groups_j.remove(rgj)
        break
    return n_linked, rgj
  # Find all isolated single residues, and also save starters
  starters = []
  for rgi in residue_groups_i:
    ci = get_atom(atoms=rgi.atoms(), name="C")
    ni = get_atom(atoms=rgi.atoms(), name="N")
    if(ci is None or ni is None): # collect non-proteins
      c = iotbx.pdb.hierarchy.chain(id = all_chains.next())
      c.append_residue_group(rgi)
      super_sphere_hierarchy.models()[0].append_chain(c)
      continue
    residue_groups_k.append(rgi)
    c_linked = 0
    n_linked = 0
    for rgj in residue_groups_j:
      cj = get_atom(atoms=rgj.atoms(), name="C")
      nj = get_atom(atoms=rgj.atoms(), name="N")
      if(cj is None or nj is None): continue
      d_ci_nj = dist(ci,nj)
      d_ni_cj = dist(ni,cj)
      if(d_ci_nj<link_max and d_ci_nj>link_min):
        n_linked += 1
      if(d_ni_cj<link_max and d_ni_cj>link_min):
        c_linked += 1
    assert c_linked in [1, 0], c_linked # Either linked or not!
    assert n_linked in [1, 0], n_linked # Either linked or not!
    if(c_linked==0 and n_linked==0): # isolated single residue
      c = iotbx.pdb.hierarchy.chain(id = all_chains.next())
      rgi.link_to_previous=True
      c.append_residue_group(rgi)
      super_sphere_hierarchy.models()[0].append_chain(c)
    elif(c_linked==0): # collect c-terminal residues
      starters.append(rgi)
  # Grow continuous chains from remaining residues starting from c-terminals
  for rgi in starters:
    ci = get_atom(atoms=rgi.atoms(), name="C")
    chunk = [rgi]
    n_linked=1
    while n_linked==1:
      n_linked, rgj = grow_chain(residue_groups_k, chunk, ci)
      if(n_linked==1):
        ci = get_atom(atoms=rgj.atoms(), name="C")
    c = iotbx.pdb.hierarchy.chain(id = all_chains.next())
    for i, rg in enumerate(chunk):
      rg.resseq = "%4d" % i
      c.append_residue_group(rg)
    super_sphere_hierarchy.models()[0].append_chain(c)
  #
  box = uctbx.non_crystallographic_unit_cell_with_the_sites_in_its_center(
    sites_cart   = super_sphere_hierarchy.atoms().extract_xyz(),
    buffer_layer = 10)
  cs_super_sphere = box.crystal_symmetry()
  return group_args(
    hierarchy        = super_sphere_hierarchy,
    crystal_symmetry = cs_super_sphere,
    selection_r                = r)

class expand(object):
  def __init__(self, pdb_hierarchy, crystal_symmetry, select_within_radius=15,
               create_restraints_manager=True):
    # fixed members
    self.create_restraints_manager = create_restraints_manager
    self.crystal_symmetry = crystal_symmetry
    self.select_within_radius = select_within_radius
    self.cs_p1 = crystal.symmetry(self.crystal_symmetry.unit_cell(), "P1")
    self.sel_str_focus = " or ".join([
      "chain %s"%c.id for c in pdb_hierarchy.chains()])
    # variable members (can change upon self.update() call) or by end of
    # constructor execution
    self.pdb_hierarchy = pdb_hierarchy
    self.super_sphere_geometry_restraints_manager = None
    self.ph_p1           = None
    self.selection_keep  = None
    self.ph_super_cell   = None
    self.cs_box          = None
    self.ph_super_sphere = None
    #
    self.update()

  def update(self, sites_cart=None, selection_r=None):
    if(sites_cart is not None):
      self.pdb_hierarchy.atoms().set_xyz(sites_cart)
    o = create_super_sphere(
      pdb_hierarchy        = self.pdb_hierarchy,
      crystal_symmetry     = self.crystal_symmetry,
      select_within_radius = self.select_within_radius,
      r                    = selection_r)
    self.ph_super_sphere = o.hierarchy
    self.ph_super_cell = self.ph_super_sphere.deep_copy() # nonsense, we don't have super-cell anymore
    self.cs_box = o.crystal_symmetry
    self.selection_r = o.selection_r
    if(self.create_restraints_manager):
      if(selection_r is None): self.update_super_sphere_geometry_restraints_manager()
    return self

  def update_xyz(self, sites_cart=None):
    # keep the selection of super_sphere, update its coordinates
    self.update(sites_cart=sites_cart, selection_r=self.selection_r)
    return self

  def update_super_sphere_geometry_restraints_manager(self):
    # XXX Unify with process_model_file of qr.py
    params = mmtbx.monomer_library.pdb_interpretation.master_params.extract()
    params.use_neutron_distances = True
    params.restraints_library.cdl = False
    params.sort_atoms = False
    processed_pdb_file = mmtbx.monomer_library.pdb_interpretation.process(
      mon_lib_srv              = mon_lib_srv,
      ener_lib                 = ener_lib,
      params                   = params,
      pdb_hierarchy            = self.ph_super_sphere,
      strict_conflict_handling = False,
      crystal_symmetry         = self.cs_box,
      force_symmetry           = True,
      log                      = null_out())
    xrs = processed_pdb_file.xray_structure()
    sctr_keys = xrs.scattering_type_registry().type_count_dict().keys()
    has_hd = "H" in sctr_keys or "D" in sctr_keys
    geometry = processed_pdb_file.geometry_restraints_manager(
      show_energies                = False,
      assume_hydrogens_all_missing = not has_hd,
      plain_pairs_radius           = 5.0)
    self.super_sphere_geometry_restraints_manager = mmtbx.restraints.manager(
       geometry = geometry, normalization = False)

  def write_super_cell_selected_in_sphere(self, file_name="super_sphere.pdb"):
    self.ph_super_sphere.write_pdb_file(file_name = file_name,
      crystal_symmetry = self.cs_box) # box around super-sphere, not original cs
                                      # or cs in p1.
