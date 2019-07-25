#!/usr/bin/python
import sys,os
import fpdb
import gmx_top

nonbonded = 'amber99sb.ffnonbonded.itp'
rtp = 'amber99sb.aminoacids.rtp'
waterfile = 'tip3p.itp'

gmxtop = gmx_top.gmxtop(nonbonded,rtp,waterfile)

pdb = fpdb.fPDB(sys.argv[1])
r1 = fpdb.topology.residues[0]
r2 = fpdb.topology.residues[1]

atom1 = r1.atoms['CA']

vs =cs = 0
for atom2 in r2.atoms:
    v,c = fpdb.potential_atom_atom(atom1,atom2)
    vs += v
    cs += c

vdw, chg = fpdb.potential_resi(r1,r2)

vdw,chg = fpdb.potential_atom_atom(atom1,atom2)

print ">>>>> vDW:",vdw
print ">>>>> Charge:",chg


