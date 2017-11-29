#!/usr/bin/python
import sys,os
import fpdb

class fSPA_water(fpdb.fCHEMO):
    def __init__(self,three_lines):
        first_line = three_lines[0]
        items = first_line.split()
        x = float(items[1])
        y = float(items[2])
        z = float(items[3])
        self.posi = (x,y,z)
        self.occ = float(items[5])
        self.index = int(items[4])
        self.residence_t = float(items[6])
        self.vdw_sol = float(items[7])
        self.ele_sol = float(items[8])
        self.vdw_rec = float(items[9])
        self.ele_rec = float(items[10])
        self.vdw_lig = float(items[11])
        self.ele_lig = float(items[12])
        self.trans_entropy = float(items[13])
        self.orient_entropy = float(items[14])
        self.spa_energy = float(items[15])
        
        ow = fpdb.fATOM()
        ow.posi = (x,y,z)

        h1 = fpdb.fATOM()
        second_line = three_lines[1]
        items = second_line.split()
        h1.posi = ( float(items[1]) , float(items[2]) , float(items[3]) )

        h2 = fpdb.fATOM()
        third_line = three_lines[2]
        items = third_line.split()
        h2.posi = ( float(items[1]) , float(items[2]) , float(items[3]) )

        self.atoms = [ow,h1,h2]
        
class fSPA_summary:
    @staticmethod
    def next_three_lines(summaryfile):
        three_lines = list()
        ifp = open(summaryfile)
        for line in ifp:
            if 'OW' in line:
                three_lines = [line,]
            elif 'H2' in line:
                three_lines.append(line)
                yield three_lines
            elif 'H1' in line:
                three_lines.append(line)

    def __init__(self,summaryfile):
        self.waters = list()
        self.waters_d = dict()
        for three_lines in fSPA_summary.next_three_lines(summaryfile):
            try:
                water = fSPA_water(three_lines)
            except Exception as e:
                print e
                sys.stderr.write("###### Error in reading SPA summary line:\n")
                sys.stderr.write("-----# %s"%three_lines[0])
                continue
            self.waters.append(water)
            self.waters_d[water.index] = water
        
class fSPA():
    def __init__():
        pass

def prepare_md(dirpath):
    dirpath = "%s/%s"%("/home/fuqy/work/SPA_database/www/submit_jobs",dirpath)
    print dirpath
    assert os.path.isdir(dirpath)
    assert os.path.isfile("%s/rec.pdb"%dirpath)
    assert os.path.isfile("%s/lig.pdb"%dirpath)

    os.chdir(dirpath)

    if True:
    ## receptor

        # assert(1==0)
        # load receptor
        rec = fpdb.fPDB("%s/rec.pdb"%dirpath)

        # fix his name 
        if os.path.isfile("%s/his.list"%dirpath):
            # if his.list exist 
            # load his.list
            histype=dict()
            for line in open("%s/his.list"%dirpath):
                histype[int(line.split()[0])] = line.split()[1]
            for resi in rec.topology.residues:
                if resi.index in histype.keys():
                    resi.name = histype[resi.index]
        else:
            # else
            for resi in rec.topology.residues:
                if resi.name in ("HIS","HID","HIP","HIE"):
                    #use HID by default
                    resi.name = "HID"

        print ">>>> TEST"
        rec.write_pdb("ftmp.pdb")
        # tleap addH
        # ofp = open("leap.in",'w')
        # ofp.write()

        # amberH 2 gmx H

        # pdb2gmx generate topology

        # add box

    ## ligand
        # fix atom coordinate, nothing to do 
   

    ## minimize receptor 
    
    ## add solvent 

    ## minimize solvent 

    ## check rmsd
        ## if too large of some atoms, warning 

    ## generate sge files

    ## return to cgi, wait for submit 

def run_spa_md():
    pass
    ## submit job
    ## wait 10 seconds
    ## check status
    ## submit to database

