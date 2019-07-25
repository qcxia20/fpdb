#!/usr/bin/python
import fpdb

def next_sdflines(infile):
    if hasattr(infile,'readlines'):
        lines = infile.readlines()
    elif infile != None :
        lines = open(infile).readlines()
    else:
        lines = list()
    this_frame = list()
    for line in lines:
        if "$$$$" in line:
            yield this_frame
            this_frame = list()
        else:
            this_frame.append(line)
        

class fSDF_ARCHIVE:
    def __init__(self,infile = None):
        self.molecules = list()
        for frame in next_sdflines(infile):
            self.molecules.append(fSDF(frame))


class fSDF:
    def __init__(self,infile = None):
        if hasattr(infile,'readlines'):
            lines = infile.readlines()
        elif infile != None :
            lines = infile
        else:
            infile = list()

        try :
            self.titleline = lines.pop(0)
            self.commentline = lines.pop(0)
            self.blankline = lines.pop(0)
            countsline = lines.pop(0)
            items = countsline.split()
            assert items[-1].strip()=='V2000'
            N_atomline = int(items[0])
            N_bondline = int(items[1])

            atomlines = list()
            for i in range(N_atomline):
                atomlines.append(lines.pop(0))
            bondlines = list()
            for i in range(N_bondline):
                bondlines.append(lines.pop(0))
            
            self.atoms = list()
            index = 1
            for line in atomlines:
                tmp = fpdb.fATOM()
                tmp.index  = index
                x = float(line.split()[0])
                y = float(line.split()[1])
                z = float(line.split()[2])
                tmp.posi = (x,y,z)
                tmp.element = line.split()[3]
                index += 1
                self.atoms.append(tmp)
            
        except Exception as e:
            print("##### Error Loading SDF file %s "%str(infile))
            print(e.message)
    def save_xyz(self,out = None,name=None):
        print( out )
        if hasattr(out,'write'):
            ofp = out
        elif hasattr(out,'upper'):
            ofp = open(out,'w')
        elif out == 'None':
            ofp = open("TMP_SDF_2_XYZ.TMP",'w')
        else:
            pass
        # N atom
        ofp.write("%d\n"%len(self.atoms))
        # Title
        if name == None:
            ofp.write(self.titleline)
        else:
            ofp.write("%s\n"%name)
        # atom
        for atom in self.atoms:
            ofp.write("%-2s %15.5f%15.5f%15.5f\n"%(atom.element,atom.posi[0],atom.posi[1],atom.posi[2]))
        # End
        ofp.flush()
        
        
