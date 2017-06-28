#!/bin/env python
import math
#import simtk.openmm.app as soa
#import simtk.unit as su
import sys,os

MAX = 99999
TMPFILE = "FPDB.SOA.TMPPDBFILE.PDB"
kJ_to_kcal = 1/4.184
BIG_NUM = 1.e10

if True: ### residue names 
    standard_protein_residues = ('ARG','LYS','HIS','HIP','HIE',
         'HID','ASP','GLU',
         'ASN','GLN','SER','THR','CYS','CYX','GLY','PRO','ALA',
         'VAL','LEU','ILE','MET','PHE','TYR','TRP' )

    standard_ion_redsidues = ('FE','MN','NI','CA','ZN','MG')
    standard_water_redsidues = ('HOH','SOL','WAT')
    standard_redsidues = standard_protein_residues + standard_water_redsidues + standard_ion_redsidues

    protein_hbond_donors = (
        ('*','N','H'),
        ('*','N','HN'),
        ('ARG', 'NE','HE'),
        ('ARG','NH1','HH11'),
        ('ARG','NH1','HH12'),
        ('ARG','NH2','HH21'),
        ('ARG','NH2','HH22'),
        ('LYS', 'NZ', 'HZ1'),
        ('LYS', 'NZ', 'HZ2'),
        ('LYS', 'NZ', 'HZ3'),
        ('HIS','ND1', 'HD1'),
        ('HIS','NE2', 'HE2'),
        ('HIP','ND1', 'HD1'),
        ('HIP','NE2', 'HE2'),
        ('HID','ND1', 'HD1'),
        ('HIE','NE2', 'HE2'),
        ('ASN','ND2','HD21'),
        ('ASN','ND2','HD22'),
        ('GLN','NE2','HE22'),
        ('GLN','NE2','HE21'),
        ('SER','OG','HG'),
        ('THR','OG1','HG1'),
        ('CYS','SG','HG'),
        ('TYR','OH','HH'),
        ('TRP','NE1','HE1'),
    )

    one_letter_code = {   
               'ala'.upper() : 'A' ,  
               'arg'.upper() : 'R' ,
               'asn'.upper() : 'N' ,
               'asp'.upper() : 'D' ,
               'asx'.upper() : 'B' ,
               'cys'.upper() : 'C' ,
               'glu'.upper() : 'E' ,
               'gln'.upper() : 'Q' ,
               'glx'.upper() : 'Z' ,
               'gly'.upper() : 'G' ,
               'his'.upper() : 'H' ,
               'ile'.upper() : 'I' ,
               'leu'.upper() : 'L' ,
               'lys'.upper() : 'K' ,
               'met'.upper() : 'M' ,
               'phe'.upper() : 'F' ,
               'pro'.upper() : 'P' ,
               'ser'.upper() : 'S' ,
               'thr'.upper() : 'T' ,
               'trp'.upper() : 'W' ,
               'tyr'.upper() : 'Y' ,
               'val'.upper() : 'V' ,
    }
    three_letter_code = dict()
    for key,value in one_letter_code.iteritems():
        three_letter_code[value] = key


if True: ### Global varieties 
    newnew = '/home/fuqiuyu/.ffallrain/newnew'
    TMPFILE = 'FIND_RING.PDB'
    TMPFILE_MOL2 = 'FIND_RING.mol2'
    babel = '/usr/bin/babel'

class fATOM():
    def __init__(self,atom_line = None):

        if atom_line == None:
            atom_line = "ATOM      1  X   DEF     1       0.000   0.000   0.000  1.00  0.00           X"        

        tmpname = atom_line[12:16]
        if tmpname[0] in "0123456789" and tmpname[3]!=" ":
            tmpname = tmpname[1:] + tmpname[0]
        tmpname = tmpname.strip()
        self.name = tmpname
        self.element = tmpname[0]  ### !!!! NOT FINISHED
        self.charge = None
        self.sig = None
        self.eps = None
        self.conf = atom_line[16]
        self.index = int(atom_line[6:11])
        x = float(atom_line[30:38])
        y = float(atom_line[38:46])
        z = float(atom_line[46:54])
        self.posi = [x,y,z]
        try:
            self.bf =  float(atom_line[60:66])
            self.occ =  float(atom_line[54:60])
        except:
            self.bf = 0
            self.occ = 0
        self.connect_count = 0 ## for bond
        self.bond = list()
        self.var1 = None
        self.var2 = None
        self.var3 = None

    def addparm(self,charge,sig,eps):
        self.charge = charge
        self.sig = sig
        self.eps = eps

    def _addconnect(self):
        self.connect_count += 1

    def _delconnect(self):
        self.connect_count -= 1
        if self.connect_count < 0 :
            sys.stderr.write("WARNING,BOND CONNECTION LESS THAN ZERO.\n")
    def addbond(self,bond):
        self.bond.append(bond)
        self._addconnect()
    def delbond(self,bond):
        self.bond.remove(bond)
        self._delconnect()

class fCHEMO():
    @staticmethod
    def _next_atom_line(resi_lines):
        for atom_line in resi_lines:
            if len(atom_line)>=6 and atom_line[:6] in ('HETATM','ATOM  '):
                yield atom_line
    def __init__(self,resi_lines=None):
        if resi_lines == None:
            resi_lines = list()
        try:
            self.name = resi_lines[0][17:20].strip()
            self.index = int(resi_lines[0][22:26])
            self.chain = resi_lines[0][21]
            self.insertion = resi_lines[0][26]
        except:
            self.name = "UND"
            self.index = 0
            self.chain = " "
            self.insertion = " "

        self.atoms = list()
        self.index_shift = 0
            
        for atom_line in fRESIDUE._next_atom_line(resi_lines):
            self.atoms.append( fATOM(atom_line) )
        if len(self.atoms)>0:
            self.index_shift = self.atoms[0].index - 1
        else:
            self.index_shift = 0
        self.atoms_d = dict()
        for atom in self.atoms:
            self.atoms_d[atom.name] = atom
        self.var1 = None
        self.var2 = None
        self.var3 = None

    def add_atom(self,atom):
        if hasattr(atom,'name'):
            self.atoms.append(atom)
            self.atoms_d[atom.name] = atom
        else :
            try:
                a = fATOM(atom)
                self.atoms.append( a )
                self.atoms_d[a.name] = a
            except:
                sys.stderr.write("##### Error, Error loading atom %s"%atom)
                

    def find_atoms(self, name = None , index = None ):
        result = list()
        ignore_name = False
        ignore_index = False
        if name == None:
            ignore_name = True
        if index == None:
            ignore_index = True
        for atom in self.atoms:
            if ignore_name or atom.name in name:
                if ignore_index or atom.index in index:
                        result.append(atom)
        return result
    def write_pdb(self,ofile):
        if hasattr(ofile,'write'):
            ofp = ofile
            for atom in self.atoms:
                x,y,z = atom.posi
                line = 'ATOM  %5d %4s%1s%3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n'%(atom.index,atom.name,atom.conf,self.name,self.chain,self.index,x,y,z,atom.occ,atom.bf)
                ofp.write(line)
        else:
            ofp = open(ofile,'w')
            for atom in self.atoms:
                x,y,z = atom.posi
                line = 'ATOM  %5d %4s %3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n'%(atom.index,atom.name,self.name,self.chain,self.index,x,y,z,atom.occ,atom.bf)
                ofp.write(line)
            ofp.close()
    
    def write_pdb_plop(self,ofile):
        if self.name in standard_protein_residues:
            atomline='ATOM  %5d %4s%1s%3s %1s%4d    %8.3f%8.3f%8.3f  1.00  0.00\n'
        else:
            atomline='HETATM%5d %4s%1s%3s %1s%4d    %8.3f%8.3f%8.3f  1.00  0.00\n'
        if self.name in standard_ion_redsidues:
            atomline='HETATM%5d %-4s%1s%3s %1s%4d    %8.3f%8.3f%8.3f  1.00  0.00\n'
        if hasattr(ofile,'write'):
            ofp = ofile
            if self.name not in standard_protein_residues:ofp.write('TER\n')
            for atom in self.atoms:
                x,y,z = atom.posi
                line = atomline%(atom.index,atom.name,atom.conf,self.name,self.chain,self.index,x,y,z)
                ofp.write(line)
            if self.name not in standard_protein_residues:ofp.write('TER\n')
        else:
            ofp = open(ofile,'a')
            if self.name not in standard_protein_residues:ofp.write('TER\n')
            for atom in self.atoms:
                x,y,z = atom.posi
                line = atomline%(atom.index,atom.name,atom.conf,self.name,self.chain,self.index,x,y,z)
                ofp.write(line)
            if self.name not in standard_protein_residues:ofp.write('TER\n')
            ofp.close()

    def debug(self):
        print 'name',self.name
        print 'index',self.index
        print 'insertion',self.insertion
        print 'atoms:'
        for atom in self.atoms:
            print atom.name,
        print

    def find_h(self):
        h_list = list()
        for atom in self.atoms:
            if atom.element == 'H':
                h_list.append(atom)
        return h_list

    def _find_bond_atoms(self,atom,cutoff):
        cutoff_2 = cutoff ** 2
        a_list = list()
        for a in self.atoms:
            if a in self.atoms:
                continue
            if dist_2(atom,a) < cutoff_2:
                a_list.append(a)
        return a_list

    def find_polar_h(self):
        polar_h_list = list()
        for atom in self.h_list():
            for a in self._find_bond_atoms(atom, 1.2): # bond length with H is usually < 1.0 A
                if a.element in ("N","O","S"): # currently only take O,N and S into account.
                    polar_h_list.append(atom)
        return polar_h_list

    def find_hbond_acceptor(self):
        l = list()
        for atom in self.atoms:
            if  atom.element in ("N","O"):
                l.append(atom)
        return l
                    

class fRESIDUE(fCHEMO):
    pass

class fRING(fCHEMO):
    def __init__(self,atom_list):
        fCHEMO.__init__(self,[])
        for atom in atom_list:
            self.atom.append(atom)

class fBOND:
    def __init__(self,atoma,atomb):
        # print "DEBUG ADD BOND BETWEEN",atoma.name,'and',atomb.name
        if atoma.name<atomb.name:
            atoma.addbond(self)
            atomb.addbond(self)
            self._ = (atoma,atomb)
        elif atoma.name>atomb.name:
            atoma.addbond(self)
            atomb.addbond(self)
            self._ = (atomb,atoma)
        else:
            pass
            sys.stderr.write("WARNING,INTENDED TO CONNECT BOND BETWEEN IDENTICAL ATOMS\n")
    def del_me(self):
        self._[0].delbond(self)
        self._[1].delbond(self)

class fCOMPOUND(fCHEMO):
    BOND_CUTOFF_H = 1.5 # Angstrom
    BOND_CUTOFF_H_2 = 2.25 # Angstrom
    BOND_CUTOFF = 2.0 # Angstrom
    BOND_CUTOFF_2 = 4.0 # Angstrom
    def __init__(self,compound_lines):
        fCHEMO.__init__(self,compound_lines)
        self._atom_classification()
        self._makebond()
    def _atom_classification(self):
        self.h_list = list()
        self.heavy_list = list()
        for atom in self.atoms:
            if len(atom.name)>=1 and atom.name[0] == 'H' :
                self.h_list.append(atom)
            else:
                self.heavy_list.append(atom)
    def _makebond(self):
        self.bonds = set()
        self._make_h_bond()
        self._make_heavy_bond()
    def _make_h_bond(self):
        for h in self.h_list:
            for heavy in self.heavy_list:
                if dist_2(h,heavy) <= self.BOND_CUTOFF_H_2:
                    self.bonds.add(fBOND(h,heavy))
                    break
    def _make_heavy_bond(self):
        for heavy_1 in self.heavy_list:
            for heavy_2 in self.heavy_list:
                if heavy_1.name < heavy_2.name :
                    if dist_2(heavy_1,heavy_2) <= self.BOND_CUTOFF_2:
                        self.bonds.add(fBOND(heavy_1,heavy_2))
    def truncate_leaf(self):
        while True:
            to_be_deleted = list()
            for atom in self.atoms:
                assert atom.connect_count > 0
                if atom.connect_count == 1:
                    to_be_deleted.append(atom)
            if len(to_be_deleted) == 0 :
                break
            else:
                for atom in to_be_deleted:
                    bonds = atom.bond
                    for bond in bonds:
                        bond.del_me()
                        self.bonds.remove(bond)
                    self.atoms.remove(atom)
        self._atom_classification()
    def find_ring(self):
        pass
        self.write_mol2(TMPFILE_MOL2)
        answer = os.popen("%s %s "%(newnew,TMPFILE_MOL2)).readlines()
        os.remove(TMPFILE_MOL2)
        rings = list()
        for line in answer:
            tmpring = list()
            numbers = [ int(x)+self.index_shift for x in line.replace(',' , '').split()[2:] ]
            tmpring = self.find_atoms(index = numbers)
            rings.append(tmpring)
        self.rings = self._merge_ring(rings)
        return self.rings
        #analyze answer and tranlate to my class  ( fRING class )
    def _merge_ring(self,rings):
        new_rings = list()
        for ring in rings:
            flag_keep = True
            for i in range(len(new_rings)):
                cp_ring = new_rings[i]
                if len(ring)+len(cp_ring)-len( set( ring + cp_ring )) >=2 :
                    new_rings[i] = list( set(ring+cp_ring) )
                    flag_keep = False
                    break
            if flag_keep:
                new_rings.append(ring)
        return new_rings
    def write_mol2(self,ofile):
        self.write_pdb(TMPFILE)
        os.system("%s -ipdb %s -omol2 %s"%(babel,TMPFILE,ofile))
        os.system("echo '@<TRIPOS>ATOM' >> %s"%ofile)
        os.remove(TMPFILE)
    def list_connect_count(self):
        for atom in self.atoms:
            print atom.name,atom.connect_count
    def debug(self):
        fCHEMO.debug(self)
        print "Hydrogen atoms:"
        for h in self.h_list:
            print h.name,
        print "Heavy atoms:"
        for heavy in self.heavy_list:
            print heavy.name,
        print "Num of Bonds:",len(self.bonds)
        print "Atom index:"
        for atom in self.atoms:
            print atom.index,
        print
   
class fTOPOLOGY():
    @staticmethod
    def _next_resi_lines(lines):
        resi_lines = list()
        oldresindex = None
        for line in lines:
            if len(line)<6 or line[:6] not in ('ATOM  ','HETATM'):
                continue
            else:
                resindex = line[22:26].strip()+line[21]
                if resindex == oldresindex or oldresindex == None:
                    resi_lines.append(line)
                else:
                    yield resi_lines
                    resi_lines = list()
                    resi_lines.append(line)
            oldresindex = resindex
        yield resi_lines

    def __init__(self,lines):
        self.residues = list()
        for resi_lines in fTOPOLOGY._next_resi_lines(lines):
            self.residues.append( fRESIDUE(resi_lines) )
        self.residues_d = dict()
        for resi in self.residues:
            self.residues_d[resi.index] = resi

    def get_protein_residues(self):
        prot_residues = list()
        for resi in self.residues:
            if resi.name in standard_protein_residues:
                prot_residues.append(resi)
        return prot_residues

    def get_water_residues(self):
        water_residues = list()
        for resi in self.residues:
            if resi.name in standard_water_redsidues:
                water_residues.append(resi)
        return water_residues

    def find_residues(self, name = None , index = None ,atom_index = None ):
        result = list()
        ignore_name = False
        ignore_index = False
        ignore_atom_index = False
        if name == None:
            ignore_name = True
        if index == None:
            ignore_index = True
        if atom_index == None:
            ignore_atom_index = True
        for residue in self.residues:
            if ignore_name or residue.name in name:
                if ignore_index or residue.index in index:
                    if ignore_atom_index or sum( [ atom_index.count(x.index) for x in residue.atoms] ) > 0:
                        result.append(residue)
        return result

    def write_model(self,ofp):
        for residue in self.residues:
            residue.write_pdb(ofp)
    
class fPDB:
    def __init__(self,frame):
        lines = None
        if hasattr(frame,'isalpha'):
            lines = open(frame).readlines()
        else:
            lines = frame

        for line in lines:
            if len(line)>=6 and line[:6]=="CRYST1":
                x = float(line.split()[1])
                y = float(line.split()[2])
                z = float(line.split()[3])
                self.box = x,y,z
                break

        if not hasattr(self,'box'):
            self.box = 0,0,0

        self.model_n = 1
        for line in lines:
            if len(line)>=6 and line[:6]=="MODEL ":
                try:
                    if len(line.split())==2:
                        self.model_n = int(line.split()[-1])    
                    else:
                        self.model_n = int(line[5:14])
                except:
                    pass
        
        self.topology = fTOPOLOGY(lines)

    @staticmethod
    def load_ff_param_resi(resi,gmxtop):
        resi_atoms = set( [ x.name for x in resi.atoms ] )
        for gmx_resi in gmxtop.get_resilist():
            gmx_atoms = set([ x[0] for x in gmxtop.get_resi(gmx_resi) ] )
            if resi_atoms == gmx_atoms:
                # print ">>>>> Loading GMX parameters : %s"%gmx_resi
                if resi.name != gmx_resi :
                    pass
                    # print "===== Warning : different residue names while loading parameters %s and %s"%(resi.name,gmx_resi)
                for atom in resi.atoms:
                    for gmx_atom in gmxtop.get_resi(gmx_resi):
                        if atom.name == gmx_atom[0]:
                            atom.addparm( gmx_atom[3],gmx_atom[1],gmx_atom[2] )
                return 
        for gmx_resi in gmxtop.get_resilist_amber():
            gmx_atoms = set( [ x[0] for x in gmxtop.get_resi_amber(gmx_resi) ] )

            if resi_atoms == gmx_atoms:
                pass
                # print ">>>>> Loading AMBER parameters : %s"%gmx_resi
                if resi.name != gmx_resi :
                    pass
                    # print "##### Warning : Different residue names while loading parameters %s and %s"%(resi.name,gmx_resi)
                for atom in resi.atoms:
                    for gmx_atom in gmxtop.get_resi_amber(gmx_resi):
                        if atom.name == gmx_atom[0]:
                            atom.addparm( gmx_atom[3],gmx_atom[1],gmx_atom[2] )
                return 
        # print "##### Warning : Error in loading residue parameters ",resi.name,resi.index, " set all parameter within this residue to zero "
        for atom in resi.atoms:
            atom.addparm(0,0,0)
        # print "Error loading residue parameters",resi.name,resi.index
        return

    def load_ff_params(self, gmxtop ):
        residues = list(self.topology.residues)
        for resi in residues :
            fPDB.load_ff_param_resi(resi,gmxtop)

    def check_params(self):
        for resi in self.topology.residues:
            for atom in resi.atoms:
                print atom.name,atom.charge,atom.sig,atom.eps

    def find_protein_center(self):
        n = 0
        xs = 0.
        ys = 0.
        zs = 0.
        for residue in self.topology.residues:
            if residue.name in standard_protein_residues:
                for atom in residue.atoms:
                    x,y,z = atom.posi
                    xs += x
                    ys += y
                    zs += z
                    n += 1
        if n == 0 :
            return 0,0,0
        else:
            return xs/n,ys/n,zs/n

    def center(self, center_posi = (0,0,0) ):
        for resi in self.topology.residues:
            for atom in resi.atoms:
                for i in range(3):
                    if atom.posi[i] - center_posi[i] > 0:
                        atom.posi[i] = atom.posi[i] - self.box[i]*int( ( atom.posi[i]-center_posi[i] )/self.box[i] + 0.5 )
                    else :
                        atom.posi[i] = atom.posi[i] + self.box[i]*int( ( - atom.posi[i] + center_posi[i] )/self.box[i] + 0.5 )

    def write_pdb(self,ofp):
        if hasattr(ofp,'write'):
            flag_is_handle = True
            flag_is_str = False
        else:
            ofp = open(ofp,'w')
            flag_is_handle = False
            flag_is_str = True

        ofp.write("MODEL     %4d\n"%self.model_n)
        ofp.write("CRYST1%9.3f%9.3f%9.3f  90.00  90.00  90.00 P 65 2 2     12\n"%self.box)
        self.topology.write_model(ofp)
        ofp.write("ENDMDL\n")

        if flag_is_str:
            ofp.close()

###### Function to compute vdw
def calc_vdw(a,b,dist_2):
    sig1,eps1=a.sig,a.eps
    sig2,eps2=b.sig,b.eps
    sig = 0.5 * ( sig1 + sig2 )
    sig = sig * sig
    eps = math.sqrt(eps1*eps2)
    _dist = dist_2 / 100.
    _ = ( sig/_dist ) ** 3
    return  4 * eps * ( _ * _ - _ ) 

###### Function to compute charge
def calc_chg(a,b,dist_2):
    #f_charge = 138.935485
    charge1 = a.charge
    charge2 = b.charge
    return 138.935485*charge1*charge2/(math.sqrt(dist_2)/10.)

###### dist_2 atom atom
def dist_2(a,b):
    if hasattr(a,'posi'):
        x1,y1,z1 = a.posi
    else:
        x1,y1,z1 = a
    if hasattr(b,'posi'):
        x2,y2,z2 = b.posi
    else:
        x2,y2,z2 = b
    return (x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2 

def dist(a,b):
    return math.sqrt( dist_2(a,b) )
def dist_resi_resi_2(a,b):
    assert len(a.atoms)>0 and len(b.atoms)>0
    dist = BIG_NUM
    for i in a.atoms:
        for j in b.atoms:
            tmp = dist_2(i,j)
            if tmp < dist:
                dist = tmp
    return dist
def dist_resi_resi(a,b):
    return math.sqrt(dist_resi_resi_2(a,b))
def dist_atom_resi_2(a,b):
    dist = BIG_NUM 
    if hasattr(b,'atoms'):
        assert len(b.atoms)>0
    else:
        assert len(b) > 0
    if hasattr(b,'atoms'):
        for j in b.atoms:
            tmp = dist_2(a,j)
            if tmp< dist:
                dist = tmp
    else:
        for j in b:
            tmp = dist_2(a,j)
            if tmp< dist:
                dist = tmp
    return dist

def dist_atom_resi(a,b):
    return math.sqrt( dist_atom_resi_2(a,b) )

def angle( a1,b,a2):
    dist_a1_a2_2 = dist_2(a1,a2)
    dist_a1_b = dist(a1,b)
    dist_a2_b = dist(a2,b)
    cos_angle = ( - dist_a1_a2_2 + dist_a1_b**2 + dist_a2_b**2 )/( 2*dist_a1_b*dist_a2_b )
    angle = math.acos(cos_angle)
    return angle

def potential_atom_atom( atoma,atomb ):
    d_2 = dist_2(atoma,atomb)
    vdw =  calc_vdw(atoma,atomb,d_2) 
    chg =  calc_chg(atoma,atomb,d_2)
    return vdw*kJ_to_kcal,chg*kJ_to_kcal

def potential_resi( resia,resib ):
    if resia == resib:
        # sys.stderr.write("WARNING, YOU ARE ATTEMPING TO COMPUTE THE POTENTIAL BETWEEN IDENTICAL RESIDUES. THE FUNCTION WILL IGNORE IT AND RETURN ZERO.\n")
        return 0,0
    vdw = 0
    chg = 0
    for atoma in resia.atoms:
        for atomb in resib.atoms:
            v,c = potential_atom_atom(atoma,atomb)
            vdw += v
            chg += c 
    return vdw, chg

def next_frame(filename):
    frame = list()
    for line in open(filename):
        if len(line)>=6 and line[:6] == "MODEL ":
            frame = list()
            frame.append(line)
        elif len(line)>=6 and line[:3] == "END":
            frame.append(line)
            yield frame
            frame = list()
        else:
            frame.append(line)

if False:
    def dist_atom(a,b):
        return math.sqrt( (a[0]-b[0])**2 + (a[1]-b[1])**2 + (a[2]-b[2])**2 )

    def dist_atom_resi(atom, resi):
        assert len(atom)>=3
        assert len(resi)>0
        dist = MAX
        node = None
        for atomb in resi:
            tmp = dist_atom(atom,atomb)
            if tmp < dist:
                dist = tmp
                node = atomb
        assert node != None
        return dist,node

# Test
if __name__ == '__main__':
    pass
        
