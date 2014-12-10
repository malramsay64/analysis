#!/usr/bin/python

import math
import molecule
import sys

PATH='{path}/data'
thePath='/home/malcolm/Honours/make/states'

def molFile(molecule):
    string = ''
    string += '# Defining the {0} molecule\n\n'.format(molecule.getName())
    string += '{0}  atoms\n'.format(molecule.numAtoms())
    string += '{0}  bonds\n'.format(molecule.numBonds())
    string += '{0}  angles\n\n'.format(molecule.numAngles())
    string += 'Coords\n\n'
    atomID = 1
    for atom in molecule:
        x,y = atom.getPos()
        string += '{0} {1} {2} 0\n'.format(atomID, x, y)
        atomID += 1
    string += '\nTypes\n\n'
    atomID = 1
    for atom in molecule:
        string += '{0} {1}\n'.format(atomID, atom.getType())
        atomID += 1
    if molecule.getBonds():
        string += '\nBonds\n\n'
        bondID = 1
        for a,b,dist in molecule.getBonds():
            string += "{0} {type} {atom1} {atom2}\n".format(\
                    bondID, type=1, atom1=a, atom2=b)
            bondID += 1
    if molecule.getAngles():
        string += '\nAngles\n\n'
        angleID = 1
        for a,b,c,theta in molecule.getAngles():
            string += '{0} {type} {atom1} {atom2} {atom3}\n'.format(\
                    angleID, type=1, atom1=a, atom2=b, atom3=c)
            angleID += 1
    if molecule.getBonds():
        string += "\nSpecial Bond Counts\n\n"
        for atom in molecule: 
            string += "{0!s} {1} {2} 0\n".format(\
                    atom.getID(), molecule.bondCount(atom.getID()), \
                    molecule.angleCount(atom.getID()))
        string += "\nSpecial Bonds\n\n"
        for atom in molecule:
            string += "{0!s} {1} {2}\n".format(atom.getID(), \
                    ' '.join(str(v) for v in molecule.get12(atom.getID())),\
                    ' '.join(str(v) for v in molecule.get13(atom.getID())))
    # Write to file
    if molecule.getAngles():
        filename="{shape}-{radius}-{dist}-{theta}".format(shape=molecule.getName(),radius=molecule.radius, dist=molecule.dist, theta=molecule.theta)
    else:
        filename="{shape}-{radius}-{dist}".format(shape=molecule.getName(),radius=molecule.radius, dist=molecule.dist)
    f = open('{PATH}/{filename}.mol'.format(PATH=PATH,filename=filename), 'w')
    f.write(string)
    f.close()

def dataFile(mol):
    string = ''
    string += '#Atom data for {} molecule\n\n'.format(mol.getName())
    string += '{}   atom types\n'.format(mol.numAtomTypes())
    string += '{}   extra bond per atom\n'.format(mol.maxBondCount())
    string += '{}   extra angle per atom\n\n'.format(mol.maxAngleCount())
    string += '0 {xdim:.0f}     xlo xhi\n'.format(xdim=math.sqrt(n_mol/dens))
    string += '0 {ydim:.0f}     ylo yhi\n\n'.format(ydim=math.sqrt(n_mol/dens))
    if mol.getBonds():
        string += 'Bond Coeffs\n\n'
        a,b,d = mol.getBonds()[0]
        string += '{bondID} {coeff} {dist}\n'.format(\
                bondID=1, coeff=500, dist=d)
    if mol.getAngles():
        string += '\nAngle Coeffs\n\n'
        a,b,c,theta = mol.getAngles()[0]
        string += '{angleID} {coeff} {theta}\n'.format(\
                angleID=1, coeff=5000, theta=theta*180/math.pi)
    string += '\nMasses\n\n'
    for t in mol.getAtomTypes():
        string += '{type} {mass}\n'.format(type=t.getType(), mass=1)
    string += '\nPair Coeffs\n\n'
    for t in mol.getAtomTypes():
        string += '{0} {strength} {dist}\n'.format(t.getType(), strength=1, dist=2*t.getSize())
    if mol.getAngles():
        filename="{shape}-{radius}-{dist}-{theta}".format(shape=mol.getName(),radius=mol.radius, dist=mol.dist, theta=mol.theta)
    else:
        filename="{shape}-{radius}-{dist}".format(shape=mol.getName(),radius=mol.radius, dist=mol.dist)
    f = open('{PATH}/{filename}.dat'.format(PATH=PATH,filename=filename),'w')
    f.write(string)
    f.close()

if __name__=='__main__':
    args = sys.argv
    nargs = len(args)

    if nargs > 1:
        PATH=PATH.format(path=args[1])
        if sys.argv[2] == 'Snowman':
            mol = molecule.Snowman(*args[3:-2])
            dens = float(args[-1])
            n_mol = int(args[-2])
        elif sys.argv[2] == 'Trimer':
            mol = molecule.Trimer(*args[3:-2])
            dens = float(args[-1])
            n_mol = int(args[-2])
        elif sys.argv[2] == 'Disc':
            mol = molecule.Disc()
            dens = float(args[-1])
            n_mol = int(args[-2])
        dataFile(mol)
        molFile(mol)



