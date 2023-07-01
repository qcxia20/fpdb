from pathlib import Path

def fmtPDBname(infile):
    lines = Path(infile).read_text().split("\n")
    newlines = []
    for line in lines:
        if line.startswith("ATOM") or line.startswith("HETATM"):
            newlines.append(line[:17] + "LIG" + line[20:25] + "1" + line[26:])
        else:
            newlines.append(line)

    with open(infile, 'w') as f:
        f.write("\n".join(newlines))

def getsdfchg(sdffile):
    totchg = 0
    lines = Path(sdffile).read_text().split("\n")
    for line in lines:
        if line.startswith('M  CHG'):
            line = line.strip()
            line = line.split()
            atom_num = int(line[2])
            counter = 4
            charge = []
            for i in range(atom_num):
                charge.append(int(line[counter]))
                counter +=2
            assert len(charge) == atom_num, 'something wrong'
            totchg = sum(charge)

    return totchg
