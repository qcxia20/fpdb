#!/usr/bin/env python
import sys,os
import argparse
import fpdb
from pathlib import Path

from .config import config
from .utils.utils import fmtPDBname, getsdfchg

for key,value in config["all"].items():
    os.environ[key] = str(value)
OBABEL_EXE = os.environ["OBABEL_EXE"]
PLOP_EXE = os.environ["PLOP_EXE"]
plop_data = os.environ["plop_data"]
DOCKBASE = os.environ["DOCKBASE"]
SCHUTILS = os.environ["SCHUTILS"]


def gen_lig_parm(ligresi, ligname, insdf):
    nc = getsdfchg(insdf)
    ligresi.generate_OPLS_parameters( nc = nc, version = '2005', rename_atom = True, removeH=False )

    if not os.path.isfile(f"{ligname.lower()}_parameter/{ligname.lower()}") :
        print( f"Net charge is {nc}, can not be used to generate proper ligand parameters")
        sys.exit()
    print( f"Net charge is {nc}.")
    
    os.system(f"cp {ligname.lower()}_parameter/{ligname.lower()} ./")


def prep_lig(insdf,hasligpdb=False):
    ligname = 'LIG'
    if not hasligpdb:
        os.system(f"{OBABEL_EXE} -isdf {insdf} -opdb -O xtal-lig.pdb")
        # Since obabel will identify amino acid in a and assign resname, we need to reformat
        fmtPDBname("xtal-lig.pdb")

    element_index = dict()
    newlines = list()
    number = 1
    for line in open('xtal-lig.pdb'):
        if "ATOM" not in line and "HETATM" not in line:
            newlines.append(line)
        else:
            name = line[12:16].strip()
            element = name[0]
            if element not in element_index:
                element_index[element] = 0
                newname = element+"%d"%element_index[element]
            else:
                element_index[element] += 1
                newname = element+"%d"%element_index[element]
            newline = "%s%5d%s %-3s%s"%(line[:6],number,line[11:12],newname,line[16:])
            number += 1
            newline = newline[:17] + ligname + newline[20:] # to deal with ATOM/HETATM PART, to change "VWW" PDB identifier into "LIG"
            newlines.append(newline)

    ligresi = fpdb.fPDB(newlines).topology.residues[0]
    # topology with default(fragmentation=False) will use fTOPOLOGY
    # fTOPOLOGY will only take line starting with "ATOM " or "HETATM" as input
    # fTOPOLOGY use fRESIDUE(fCHEMO) method for residues,
    # type(ligresi) = <fpdb.fpdb.fRESIDUE object at 0x7ffff0904fd0>

    gen_lig_parm(ligresi, ligname, insdf)
    return ligresi, ligname


def prep_cmx_by_icda(recpdb,insdf,hasligpdb=False,notcheckrecpdb=False):
    #### NOTE: Ligand Part #####
    ligresi, ligname = prep_lig(insdf,hasligpdb=hasligpdb)

    #### NOTE: Protein Part ####
    if not notcheckrecpdb:
        os.system(f"{SCHUTILS}/pdbconvert -ipdb {recpdb} -omae {recpdb}.mae")
        os.system(f"cp {recpdb} {recpdb}.bak")
        os.system(f"{SCHUTILS}/pdbconvert -imae {recpdb}.mae -opdb {recpdb}")
    rec = fpdb.fPDB(recpdb)

    #### NOTE: complex part ####
    with open('TMPCOM.pdb','w') as ofp: # no water, but hydrogens are kept (some labeled as HETATM)
        icda_residues = list()
        for resi in rec.topology.residues:
            if resi.name in fpdb.standard_protein_residues:
                resi.write_pdb_plop(ofp)
                if resi.name == 'HIS':
                    icda_residues.append(resi)
        ligresi.write_pdb_plop(ofp)
    # TMPCOM.pdb contains the ligand part info
    
    with open("res.list",'w') as ofp:
        # HIS
        for resi in icda_residues:
            if resi.chain !="":
                ofp.write("%s:%d\n"%(resi.chain,resi.index))
            else:
                ofp.write("_:%d\n"%(resi.index,))
    
    #### NOTE: icda part to prepare the complex ####
    with open("icda.input",'w') as ofp:
        # ofp.write("file datadir /home/qyfu/Software/plop21.0/data \n")
        # ofp.write("file datadir /home/soft/plop/25.6/data \n")
        ofp.write(f"file datadir {plop_data} \n")
        ofp.write("file logfile ICDA.log \n")
        ofp.write("load pdb TMPCOM.pdb het yes ions no wat yes opt yes seqres no \n") # no ions
        # ofp.write("load pdb TMPCOM.pdb het yes ions no wat yes opt no seqres no \n") # no ions
        ofp.write("pka clust file res.list fastpka no adjust yes \n") # use fastpka to set his protonation state
        ofp.write("minim hyd \n") # minimize hydrogen positions
        ofp.write("energy calc \n") # calculate energy
        ofp.write("structure write ICDAOUT.pdb \n")

    os.system(f"{PLOP_EXE} icda.input")
    return ligname


def modify_ligand_sphere():
    from .utils.ligsph import removeligsph
    # remove ligand sphere directly
    newlines, ligcount = removeligsph("./blastermaster/dockfiles/matching_spheres.sph")
    with open("./blastermaster/dockfiles/matching_spheres_noLIGSPH.sph", 'w') as f:
        f.write(''.join(newlines))
    # without ligand sphere, but add more receptor spheres
    os.chdir("blastermaster/working")
    os.system(f"{DOCKBASE}/proteins/makespheres3/makespheres3.cli.pl 1.5 0.8 {45+ligcount} xtal-lig.match.sph all_spheres.sph rec.crg.pdb matching_spheres_{45+ligcount}.sph >& matching_spheres_{45+ligcount}.log")
    os.chdir("../..")
    newlines, _ = removeligsph(f"./blastermaster/working/matching_spheres_{45+ligcount}.sph")
    with open("./blastermaster/dockfiles/matching_spheres_noLIGSPH_45.sph", 'w') as f:
        f.write(''.join(newlines))


def run_blastermaster(ligname,modligsph=False):
    # make blastermaster input
    # Simply use PLOP output
    os.system("mkdir -p blastermaster")
    assert os.path.isfile("ICDAOUT.pdb")
    os.system(f"cat ICDAOUT.pdb|grep -v {ligname.upper()} > blastermaster/rec.pdb")
    os.system("cp xtal-lig.pdb blastermaster")
    
    # run blastermaster
    os.chdir("blastermaster")
    with open("tmp.sh",'w') as ofp:
        ofp.write("source ~/.bashrc\n")
        ofp.write("conda activate dock37\n")
        ofp.write(f"python2 {DOCKBASE}/proteins/blastermaster/blastermaster.py --addhOptions=\" -HIS -FLIPs \" -v\n")
    os.system("bash tmp.sh &> blastermaster.log")
    os.system("rm tmp.sh")
    os.chdir("..")
    if modligsph: modify_ligand_sphere()

def main():
    parser = argparse.ArgumentParser("Prepare protein-ligand complex with PLOP ICDA and blastermaster, HIS will be predicted and hydrogens will be minimized.")
    parser.add_argument("--fixedsdf", type=str, required=True, help="fixed ligand sdf used in complex building")
    parser.add_argument("--recpdb", type=str, required=True, help="receptor pdb used in complex building")
    parser.add_argument("--hasligpdb", action="store_true", default=False, help="whether already has well-modified xtal-lig.pdb, default=False, which means will be generated")
    parser.add_argument("--notcheckrecpdb", action="store_true", default=False, help="wheter to check pdb using schrodinger pdbconvert, default=False, which means check=True")
    parser.add_argument("--modligsph", action="store_true", default=False, help="whether to modify ligand spheres, default=False")
    args = parser.parse_args()

    #### NOTE: Prepare complex using ICDA (a sub module in PLOP) ####
    ligname = prep_cmx_by_icda(
        recpdb=args.recpdb,
        insdf=args.fixedsdf,
        hasligpdb=args.hasligpdb,
        notcheckrecpdb=args.modligsph
        )

    #### NOTE: Blastermaster to prepare dockfiles ####
    run_blastermaster(ligname=ligname, modligsph=args.modligsph)

    # clean
    # remove temparoy files
    garbage_list =( "lig",   # gen OPLS parameters for lig
                    "lig_parameter",
                    "xtal-lig.pdb",
                    ################
                    # f"{args.recpdb}.bak" # backuped 6arv_protein.pdb before pdbconvert
                    f"{args.recpdb}.mae"
                    ################
                    "TMPCOM.pdb", # protein-ligand complex
                    "icda.input", # input for PLOP
                    "res.list", # HIS
                    "plop_job.pka_restart",
                    "ICDA.log",
                    "ICDAOUT.pdb",
    )
    
    for tmp in garbage_list:
        os.system(f"rm -r {tmp}")


if __name__ == "__main__":
    main()
