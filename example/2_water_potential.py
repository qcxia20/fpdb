#!/usr/bin/python
import sys,os
import fpdb
import gmx_top

pdb = fpdb.fPDB(sys.argv[1])

gmxtop = gmx_top.gmxtop("ffnonbonded.itp","aminoacids.rtp","tip3p.itp")

pdb.load_ff_params(gmxtop)


vdw,chg = fpdb.potential_resi( pdb.topology.residues[0],pdb.topology.residues[1] )

print "vDW:",vdw
print "Charge:",chg

