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
            except:
                sys.stderr.write("###### Error in reading SPA summary line:\n")
                sys.stderr.write("-----# %s"%three_lines[0])
                continue
            self.waters.append(water)
            self.waters_d[water.index] = water
        
