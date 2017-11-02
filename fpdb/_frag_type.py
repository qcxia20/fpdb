"""residues to fragment type
"""
# bck = {1:['N','H','CA','HA',
#           ('C','CA',0.9),('CB','CA',0.9)],
#        0:['C','O','CA','HA',
#           ('C','CA',0.9),('CB','CA',0.9)],
bck = [(1,'N'),(1,'H'),(0,'C'),(0,'O'),
       (0,'CA'),(0,'HA'),(0,'C'),(0,'CB'),
       (1,'CA'),(1,'HA'),(1,'C'),(1,'CB')]

wtr = {'HOH':['OW','HW1','HW2']}
alc = {'THR':['CB', ('CA','CB', 1.09), ('CG2','CB',1.09), 'HB',
              'OG1', 'HG1'],  
       'SER':['CB', ('CA','CB', 1.09), 'HB1', 'HB2', 'OG', 'HG']}

arg = {'ARG':['CZ', 'NE', 'NH1', 'NH2', ('CD','NE', 1.0), 'HE', 
              'HH11', 'HH12', 'HH21', 'HH22']}
ben = {'PHE':['CG', 'CD1', 'CE1', 'CZ', 'CE2', 'CD2', 
             ('CB','CG',1.09), 'HD1', 'HE1', 'HZ', 'HE2', 'HD2']}

co2 = {'ASP':['CB',('CA', 'CB', 1.09),'HB1','HB2','CG','OD1','OD2'],
       'GLU':['CG',('CB', 'CG', 1.09),'HG1','HG2','CD','OE1','OE2']}

cys = {'CYS':['CB',('CA','CB', 1.08),'HB1','HB2','SG','HG']}

gln = {'ASN':['CB','CG','OD1','ND2',('CA','CB',1.09),
              'HB1','HB2','HD21','HD22'],
       'GLN':['CG','CD','OE1','NE2',('CB','CG',1.09),
              'HG1','HG2','HE21','HE22']}

hid = {'HID':['CE1','ND1','CG','CD2','NE2',
              'HE1','HD1',('CB','CG',1.07),'HD2'],
       'HIE':['CE1','NE2','CD2','CG','ND1',
              'HE1','HE2','HD2',('CB','CG',1.07)]} 

hip = {'HIP':['CE1','NE2','CD2','CG','ND1',
              'HE1','HE2','HD2',('CB','CG',1.07),'HD1']}

met = {'MET':['CG',('CB','CG',1.09),'HG1','HG2','SD',
              'CE','HE1','HE2','HE3']}

nh4 = {'LYS':['NZ','HZ1','HZ2','HZ3',
              'CE','HE1','HE2',('CD','CE', 1.08)]}

prp = {'ARG':['CB',('CA','CB',1.09),'HB1','HB2','CG','HG1','HG2',
              'CD','HD1','HD2',('NE','CD',1.09)],
       'ILE':['CB',('CA','CB',1.09),('CG2','CB',1.09),'HB',
              'CG1','HG11','HG12','CD','HD1','HD2','HD3'],
       'LEU':['CD1','HD11','HD12','HD13','CG',('CB','CG',1.09),'HG',
              'CD2','HD21','HD22','HD23'],
       'LYS':['CB',('CA','CB',1.09),'HB1','HB2','CG','HG1','HG2',
              'CD','HD1','HD2',('CE','CD',1.09)],
       'PRO':['CB',('CA','CB',1.09),'HB1','HB2','CG','HG1','HG2',
              'CD','HD1','HD2',('N','CD',1.09)],
       'VAL':['CG1','HG11','HG12','HG13','CB',('CA','CB',1.09),'HB',
              'CG2','HG21','HG22','HG23']}

trp = {'TRP':['CG','CD1','CD2','NE1','CE2','CE3','CZ2','CZ3','CH2',
              'HD1','HE1','HE3','HZ2','HZ3','HH2',('CB','CG',1.08)]}
tyr = {'TYR':['CE1','CD1','CG','CD2','CE2','CZ','OH',
              'HE1','HD1',('CB','CG',1.08),'HD2','HE2','HH']}

lnk = {
  'ALA':['CB','1HB','2HB','3HB'],
  'GLU':['CB','1HB','2HB'],
  'GLN':['CB','1HB','2HB'],
  'HID':['CB','1HB','2HB'],
  'HIE':['CB','1HB','2HB'],
  'HIP':['CB','1HB','2HB'],
  'ILE':['CG2','1HG2','2HG2','3HG2'],
  'LEU':['CB','1HB','2HB'],
  'MET':['CB','1HB','2HB'],
  'PHE':['CB','1HB','2HB'],
  'THR':['CG2','1HG2','2HG2','3HG2'],
  'TRP':['CB','1HB','2HB'],
  'TYR':['CB','1HB','2HB']}


FRAG_TYPE = {
  'alc':alc,
  'arg':arg,
  'ben':ben,
  'co2':co2,
  'cys':cys,
  'gln':gln,
  'hid':hid,
  'hip':hip,
  'met':met,
  'nh4':nh4,
  'prp':prp,
  'trp':trp,
  'tyr':tyr,
  'wtr':wtr}
