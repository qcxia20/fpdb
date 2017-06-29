#!/usr/bin/python
import sys,os

class fXYZ:
    def __init__(self,frame):
        try:
            self.n_atom = int(frame[0].split()[0])
            self.n_lines = len(frame)
            self.lines = frame
            if self.lines[1].split()[0] == '1' :
                self.box = 0,0,0
            else:
                x = float(self.lines[1].split()[0])
                y = float(self.lines[1].split()[1])
                z = float(self.lines[1].split()[2])
                self.box = x,y,z
        except:
            print "ERROR format of *.xyz file"
            sys.exit()
    def write_xyz(self,ofp):
        for line in self.lines:
            ofp.write(line)
    def write_pdb(self,ofp):
        TMP_XYZ = 'fTMP.xyz'
        pass ###### alpha 


def next_frame(filename):
    frame = list()
    for line in open(filename):
        if len(line.split())==1 :
            if len(frame) == 0 :
                pass
            else:
                yield frame
                frame = list()
            frame.append(line)
        else:
            frame.append(line)
    yield frame
