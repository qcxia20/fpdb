import numpy as np
from copy import copy
from fpdb import fCHEMO
from fpdb import fTOPOLOGY
from fpdb import standard_protein_residues
from _frag_type import FRAG_TYPE


class fFRAGTOPO(fTOPOLOGY):
  def __init__(self, lines):
    fTOPOLOGY.__init__(self, lines)
    self.frags = []
    self.res2frag()

  def find_frags(self, name):
    frags = []
    for frag in self.frags:
      if frag.name == name:
        frags.append(frag)
    return frags
    
  def res2frag(self):
    # rename HIS to HID, HIE, HIP
    self.renameHIS(self.residues)
    self.add_bck()
    self.add_frag()

  def add_bck(self):
    # generate backbone fragment
    default_atom_dict = [
      (1,'N'),(1,'H'),(0,'C'),(0,'O'),
      (0,'CA'),(0,'HA'),(0,('N','CA',0.9)),(0,('CB','CA',0.9)),
      (1,'CA'),(1,'HA'),(1,('C','CA',0.9)),(1,('CB','CA',0.9))]
    pro_atom_dict = [
      (1,'N'),(1,('CD','N',0.9)),(0,'C'),(0,'O'),
      (0,'CA'),(0,'HA'),(0,('N','CA',0.9)),(0,('CB','CA',0.9)),
      (1,'CA'),(1,'HA'),(1,('C','CA',0.9)),(1,('CB','CA',0.9))]
    gly_atom_dict_0 = [
      (1,'N'),(1,'H'),(0,'C'),(0,'O'),
      (0,'CA'),(0,'HA1'),(0,('N','CA',0.9)),(0,'HA2'),
      (1,'CA'),(1,'HA'),(1,('C','CA',0.9)),(1,('CB','CA',0.9))]
    gly_atom_dict_1 = [
      (1,'N'),(1,'H'),(0,'C'),(0,'O'),
      (0,'CA'),(0,'HA'),(0,('N','CA',0.9)),(0,('CB','CA',0.9)),
      (1,'CA'),(1,'HA1'),(1,('C','CA',0.9)),(1,'HA2')]
    gly_atom_dict_01 = [
      (1,'N'),(1,'H'),(0,'C'),(0,'O'),
      (0,'CA'),(0,'HA1'),(0,('N','CA',0.9)),(0,'HA2'),
      (1,'CA'),(1,'HA1'),(1,('C','CA',0.9)),(1,'HA2')]

    n_res = len(self.residues)
    for i in range(n_res):
      if i == n_res-1: break 
      resi = self.residues[i]
      next_resi = self.residues[i+1]
      if (resi.name not in standard_protein_residues
        or next_resi.name not in standard_protein_residues):
        continue
      if next_resi.index - resi.index != 1: 
        continue
      tmp_res = [resi, next_resi]
      f = fCHEMO()
      f.name = 'bck'
      atom_dict = default_atom_dict
      if next_resi.name == 'PRO':
        atom_dict = pro_atom_dict
      if resi.name == 'GLY':
        atom_dict = gly_atom_dict_0
      if next_resi.name == 'GLY':
        atom_dict = gly_atom_dict_1
      if resi.name == 'GLY' and next_resi.name == 'GLY':
        atom_dict = gly_atom_dict_01
      for idx, atom_name in atom_dict:
        if type(atom_name) is tuple:
          atom = self.createH(tmp_res[idx], atom_name)
        else:
          atom = tmp_res[idx].atoms_d[atom_name]
        f.add_atom(atom)
      f.index = len(self.frags)+1
      self.frags.append(f)

  def add_frag(self):
    # generate non-backbone fragment
    for frag_name, res2frag_dict in FRAG_TYPE.items():
      for resi in self.residues:
        if resi.name in res2frag_dict:
          f = fCHEMO()
          f.name = frag_name
          atom_dict = res2frag_dict[resi.name]
          for atom_name in atom_dict:
            if type(atom_name) is tuple:
              atom = self.createH(resi, atom_name)
            else:
              atom = resi.atoms_d[atom_name]
            f.add_atom(atom)
          f.index = len(self.frags)+1
          self.frags.append(f)

  @staticmethod 
  def createH(resi, name_tuple):
    H, C, d = name_tuple
    H = resi.atoms_d[H]
    C = resi.atoms_d[C]
    h = copy(H)
    vec = np.subtract(H.posi, C.posi)
    vec /= np.linalg.norm(vec)
    posi = C.posi + vec * d
    h.posi = list(posi)
    h.element = 'H'
    h.name = 'H' + h.name
    return h

  @staticmethod
  def renameHIS(residues):
    for resi in residues:
      if resi.name == 'HIS':
        if 'HE2' in resi.atoms_d:
          resi.name = 'HIE'
        if 'HD1' in resi.atoms_d:
          if resi.name == 'HIE':
            resi.name = 'HIP'
          else:
            resi.name = 'HID'
      
      
      
      
