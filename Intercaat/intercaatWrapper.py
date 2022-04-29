import intercaat.intercaat_functions as icaat
import sys

def intercaat(pdb: str, qc: str, ic: str, mi: int = 4,di: str = "yes",cc: str = "yes", sr: float = 1.4, vi = [], fp: str = "./", qhull = False):
    arg1 = pdb
    arg2 = qc.split(',')
    arg3 = ic.split(',')
    arg4 = int(mi)
    arg5 = di
    arg6 = cc
    arg7 = float(sr)
    arg8 = vi
    arg9 = fp
    if arg8 != []:
        arg8 = arg8.split(',')
    try:
        pdb  = icaat.parse(arg1, (arg2+arg3 + arg8), arg9)
    except FileNotFoundError:
        sys.exit("filename not found")
    coordinates = []
    match = []

    # Parse PDB file into list of PDB lines containing atoms and hetatms excluding hydrogen and water.
    for line in pdb:
        coordinates.append([line[8], line[9], line[10]])

    # Creates 3D voronoi diagram and returns indices of neighboring atoms
    if qhull:
        contacts = icaat.voroC(coordinates)
    else:
        contacts = icaat.voroPyhull(coordinates)
    # Creates a list (pdbAtomClass) that contains classes for all atoms analayzed
    pdbAtomClass = icaat.appendAtomClasses(pdb)

    # Creates list of all atomic interactions
    count1 = 0
    count2 = 0
    for buddies in contacts:
        while count2 < len(buddies):
            # XYZ = atom, X, Y, Z
            XYZ1 = [pdb[count1][2][0:2], float(pdb[count1][8]), float(pdb[count1][9]), float(pdb[count1][10])]
            XYZ2 = [pdb[buddies[count2]][2][0:2], \
            float(pdb[buddies[count2]][8]), float(pdb[buddies[count2]][9]), float(pdb[buddies[count2]][10])]
            # Returns the distance between two atoms and min distance reqeuired for solvent molecule
            Ad, Vd = icaat.inter(XYZ1,XYZ2,arg7)
            # Finds class of atom1 and atom2
            class1 = pdbAtomClass[count1]
            class2 = pdbAtomClass[buddies[count2]]
            # atom = residue, residue #, chain, atom
            atom1 = '{0:<3} {1:>5} {2} {3:<4}'.format(pdb[count1][4], pdb[count1][6], pdb[count1][5], pdb[count1][2])
            atom2 = '{0:<3} {1:>5} {2} {3:<4}'.format(pdb[buddies[count2]][4], pdb[buddies[count2]][6],\
                     pdb[buddies[count2]][5], pdb[buddies[count2]][2])
            Line = '{0} | {1} | {2:<4} |    {3}   {4}'.format(atom1, atom2, str(round(Ad,2)), str(class1), str(class2))
            # Only appends if classes are compatible or the user inputs that class compatibility does not matter
            # if class compatibility is unknown, the atomic interaction will be shown
            if icaat.compatible(class1,class2) == True or arg6 == 'no':
                # line 1: appends if distance between atoms is within bound and accounts for occupancy < 1
                # line 2: appends only for specific chain
                if Ad < Vd and any((atom1 + atom2) in sub for sub in match) == False \
                and pdb[count1][5] == arg2[0]:
                    # appends all residues from specified chain
                    if len(arg2) == 1:
                        # appends only against queried neighbor chains
                        if pdb[buddies[count2]][5] in arg3:
                            match.append(Line)
                    #appends only specified residues from specified chain
                    elif pdb[count1][6] in arg2:
                        # appends only against queried neighbor chains
                        if pdb[buddies[count2]][5] in arg3:
                            match.append(Line)
            count2 += 1
        count2 = 0
        count1 += 1
    # Filters the list generated above (match) based on arg2 and arg4
    # Also displays interactions matrix based on arg5
    newMatch, newInteractionRes, newInteractions = icaat.filterMatch(match, pdb, arg2, arg4, arg5)
    matches = matches_to_dict(newMatch)
    interactions = interactions_to_dict(newInteractionRes, newInteractions)
    return matches, interactions


def interactions_to_dict(newInteractionRes, newInteraction):
    return  {f"{i}{l[0]}": [i,l[0],l[1]] for i, l in zip(newInteractionRes, newInteraction)}
    
def matches_to_dict(newMatch):
    # d=  {i.split("|")[0]: i.split("|")[1:] for i in newMatch}
    d = {}
    for i in newMatch:
        key  = i.split("|")[0]
        key = "".join(key.split())
        value = i.split("|")[1:]
        value = ["".join(s.split()) for s in value]
        d[key] = value
    return d
