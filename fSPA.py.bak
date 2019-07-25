#!/usr/bin/python
import sys,os
import fpdb
from subprocess import Popen as sP 
from subprocess import PIPE 


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

def cgi_system(cmd,output):
    i,o = sP(cmd,shell=True,stdout=PIPE,stderr=PIPE).communicate()
    output.write(i)
    output.write(o)
    
def prepare_md(dirpath,output=sys.stdout):
    dirpath = "%s/%s"%("/home/fuqy/work/SPA_database/www/submit_jobs",dirpath)
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

        for resi in rec.topology.residues:
            if resi.name in ("CYS","CYZ"):
                #use HID by default
                resi.name = "CYX"

        print ">>>>> OUTPUT of fSPA.prepare_md(): "
        rec.write_pdb("ftmp.pdb")
        cgi_system("addter.py ftmp.pdb ftmp_ter.pdb",output)
        
        #cgi_system("addter.py ftmp.pdb ftmp_ter.pdb").read()
        with open("leap.in",'w') as ofp:
            ofp.write("""
                source leaprc.ff14SB
                rec = loadpdb ftmp_ter.pdb
                savepdb rec tleap_out.pdb
                quit
            """)
        cgi_system("tleap -f leap.in",output)
        cgi_system("sed -i s/CYX/CYS/g tleap_out.pdb",output)
        cgi_system("atomtypeconvert.py a2g tleap_out.pdb tleap_out_trans.pdb",output)
        # cgi_system("gmx pdb2gmx -f tleap_out_trans.pdb -o rec.gro -merge all -water tip3p -ff amber99sb",output)

        # add box
        # cgi_system("gmx editconf -f rec.gro -o box.gro -d 1.0",output)

    ## ligand
    if False:
        # fix atom coordinate, nothing to do 
        # compute shift
        line_o = open('rec.gro').readlines()[2]
        line_f = open('box.gro').readlines()[2]
        xo = float(line_o[20:28]) * 10
        yo = float(line_o[28:36]) * 10
        zo = float(line_o[36:44]) * 10
        xf = float(line_f[20:28]) * 10
        yf = float(line_f[28:36]) * 10
        zf = float(line_f[36:44]) * 10
        xs,ys,zs = xf-xo,yf-yo,zf-zo 
        
        # adjust ligand coordinate 
        with open("lig_shift.pdb",'w') as ofp:
            for line in open("lig.pdb"):
                if len(line)>=6 and line[:6] in ("ATOM  ","HETATM"):
                    x = float(line[30:38]) + xs
                    y = float(line[38:46]) + ys
                    z = float(line[46:54]) + zs
                    ofp.write("%s%8.3f%8.3f%8.3f%s"%(line[:30],x,y,z,line[54:]))
                else:
                    ofp.write(line)

    ## cp mdp files
    # cgi_system("cp -r %s/mdp ./"%fpdb.__path__[0][:-4],output)
   
    ## minimize receptor 
    # cgi_system("gmx grompp -f mdp/em.mdp -o em.tpr -c box.gro -p topol.top -maxwarn 10",output)
    # cgi_system("gmx mdrun -deffnm em",output)
    
    ## add solvent 
    # cgi_system("gmx solvate -cp em.gro -cs spc216 -o sys.gro -p topol.top",output)

    ## minimize solvent 
    # cgi_system("gmx grompp -f mdp/em.mdp -o em_sys.tpr -c sys.gro -p topol.top -r box.gro -maxwarn 10",output)
    #cgi_system("gmx mdrun -deffnm em_sys",output)


    ## check rmsd
        ## if too large of some atoms, warning 

    ## generate sge files

    ## return to cgi, wait for submit 
    return 0

def submit_job(index,stype='AMOEBA',queue='gpu',output=sys.stdout):
    dir_path = "%s/www/submit_jobs/%s"%(os.environ['SPA_DATABASE_HOME'],index)
    name = open("%s/name"%dir_path).read().strip()
    assert name[0] not in "0123456789"
    if stype == "AMOEBA" and queue == 'gpu':
        cgi_system("scp -r %s fuqy@k209:/pubhome/fuqy/%s_%s"%(dir_path,name,index), output)
        # cgi_system("ssh fuqy@k209 cp /pubhome/fuqy/run.sh /pubhome/fuqy/Auto_%s/"%(name,),output )
        cgi_system("ssh fuqy@k209 bash /pubhome/fuqy/submit.sh %s %s"%(name,index),output)
    return

def run_spa_md(dirpath,output=sys.stdout):
    pass

    ## submit job
    ## SA solvent

    ## reminimize
    ## wait 10 seconds
    ## check status
    ## submit to database

