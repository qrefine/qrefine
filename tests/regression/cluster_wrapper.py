from __future__ import division
import iotbx.pdb
from qrefine.restraints import from_qm

class Clusterer(object):

  def create(self,pdb):
    self.pdb = pdb
    self.pdb_inp = iotbx.pdb.input(self.pdb)
    self.ph = self.pdb_inp.construct_hierarchy()
    self.cs = self.pdb_inp.crystal_symmetry()
    self.sites_cart = self.ph.atoms().extract_xyz()
    self.maxnum_residues_in_cluster = 12
    self.fq = from_qm(
             use_cluster_qm             = True,
             pdb_hierarchy              = self.ph,
             crystal_symmetry           = self.cs,
             maxnum_residues_in_cluster = int(self.maxnum_residues_in_cluster)
             )

  def process(self):
    chunks = []
    chunk_sizes = []
    for chunk in self.fq.fragments.qm_pdb_hierarchies:
      res_in_chunk = []
      atom_tot_per_residue = 0
      for chain in chunk.only_model().chains():
        for residue_group in chain.residue_groups():
          res_in_chunk.append(residue_group.resid())
          for atom_group in residue_group.atom_groups():
            atom_tot_per_residue += atom_group.atoms_size()
            chunk_sizes.append(atom_tot_per_residue)
      chunks.append(res_in_chunk)
    return Result(self.pdb_file,
                     self.fq.fragments.clusters,
                     chunks,
                     chunk_sizes)

  def __call__(self,pdb):
    self.create(pdb)
    return self.process()

class Result(object):
  def __init__(self, pdb_code, clusters, chunks, chunk_sizes):
    self.pdb_code        = pdb_code
    self.clusters        = clusters
    self.chunks          = chunks
    self.chunk_sizes     = chunk_sizes
    self.num_of_clusters = len(clusters)
    self.num_of_chunks   = len(chunks)
