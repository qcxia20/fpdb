import fpdb
pdb = fpdb.fPDB('waterDimer.pdb')
residues = pdb.topology.residues
minDist = 10e8
keeper = None
for i, resi in enumerate(residues):
  for j, resj in enumerate(residues):
    if j >= i: continue
    for k, atomk in enumerate(resi.atoms):
      for l, atoml in enumerate(resj.atoms):
        if l >= k: continue
        d = fpdb.dist(atomk, atoml)
        if d < minDist:
          minDist = d
          keeper = (resi, atomk, resj, atoml)

resi, atomk, resj, atoml = keeper
print("min distance atom pair are atom %s%d:%s%d and %s%d:%s%d"%(
  resi.name, resi.index, atomk.name, atomk.index, 
  resj.name, resj.index, atoml.name, atoml.index,))
  
