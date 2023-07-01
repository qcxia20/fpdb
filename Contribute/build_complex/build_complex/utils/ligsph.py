def removeligsph(matching_sph):
    with open(matching_sph, 'r') as f:
        lines = f.readlines()
    newlines = []   
    ligcount = 0
    for line in lines:
        if len(line) == 64:
            if line.strip().split()[4] == "0.000": # lig sphere
                ligcount += 1
            else:
                newlines.append(line)
        else:
            newlines.append(line)
    newlines[13] = newlines[13][:-3] + str(int(newlines[13][-3:-1]) - ligcount) + "\n"

    return newlines, ligcount

