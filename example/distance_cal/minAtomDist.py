import fpdb
pdb = fpdb.fPDB('waterDimer.pdb')
residues = pdb.topology.residues
minDist = 10e8
keeper = None
for i, resi in enumerate(residues):
  for resj in residues[:i]:
    for atomi in resi.atoms:
      for atomj in resj.atoms:
        d_2 = fpdb.dist_2(atomi, atomj)
        if d_2 < minDist:
          minDist = d_2
          keeper = (resi, atomi, resj, atomj)

resi, atomi, resj, atomj = keeper
print("\nmin distance atom pair are atom %s%d:%s%d and %s%d:%s%d\n"%(
  resi.name, resi.index, atomi.name, atomi.index, 
  resj.name, resj.index, atomj.name, atomj.index,))
  
