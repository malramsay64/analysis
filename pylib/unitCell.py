#! /usr/bin/python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2014 malcolm <malcolm@asaph-VirtualBox>
#
# Distributed under terms of the MIT license.

"""

"""
import molecule
from copy import deepcopy as copy
from math import *

class cell:
    def __init__(self, a, b, theta, molecule):
        self.a = a
        self.b = b
        self.theta = theta*pi/180
        self.mol = molecule
        self.mols = []
        self.crys = ""

    def addParticle(self,x,y,phi):
        x *= self.a
        y *= self.b
        xy = abs(y/tan(self.theta))
        x += xy
        m = copy(self.mol)
        m.setAngle(phi)
        m.setPos(x,y)
        self.mols.append(m)

    def addMol(self, x, y, phi):
        self.addParticle(x,y,phi)

    def getMol(self):
        return self.mol

    def getMols(self):
        return self.mols

    def numMols(self):
        return len(self.mols)

    def numAtoms(self):
        return len(self.mols)*self.mol.numAtoms()

    def getShape(self):
        return self.a, self.b, self.theta

    def getCrys(self):
        return self.crys
    
    def getA(self):
        return self.a
    
    def replicate(self, nx=1, ny=1):
        base = self.getMols()
        new = []
        a,b,theta = self.getShape()
        for i in xrange(nx):
            for j in xrange(ny):
                for mol in base:
                    mnew = copy(mol)
                    if theta != 0:
                        mnew.translate(i*a+j*b*tan(pi/2-theta), j*b)
                    else:
                        mnew.translate(i*a, j*b)
                    new.append(mnew)

        self.a = nx*a
        self.b = ny*b
        self.mols = new


    def __str__(self):
        mid = 1
        aid = 1
        s = ""
        for mol in self.getMols():
            for atom in mol:
                x,y = atom.getPos()
                # Wrap x coordinates
                #x = wrap(x,self.a)
                s += "{aid} {mid} {tid} {taid} {atype} {x} {y} {z}\n"\
                        .format(aid=aid, mid=mid, tid=1,\
                        taid=atom.getType(), atype=atom.getType(),\
                        x=x, y=y, z = 0)
                aid += 1
            mid += 1
        return s

    def __repr__(self):
        xy = self.b/tan(self.theta)
        s = ""
        s += "0,0\n"
        s += "{a},0\n\n".format(a=self.a)
        s += "0,0\n"
        s += "{xy},{b}\n\n".format(xy=xy, b=self.b)
        s += "{xy},{b}\n".format(xy=xy, b=self.b)
        s += "{a},{b}\n\n".format(a=self.a+xy, b=self.b)
        s += "{a},0\n".format(a=self.a)
        s += "{a},{b}\n\n".format(a=self.a+xy, b=self.b)
        for mol in self.mols:
            for atom in mol:
                x,y = atom.getPos()
                #x = wrap(x,self.a)
                s += "{x}, {y}, {size}\n".format(x=x,y=y,size=atom.getSize())
            s += "\n"
        return s

    def rotation(self, pos, degree, n=1):
        x,y = self.mols[pos].atoms[0].getPos()
        xy = self.b/tan(self.theta)
        x = self.a - x + xy
        y = self.b - y
        self.mols[pos].rotate((360/degree)*n)
        self.mols[pos].setPos(x,y)

class p2gg(cell):

    def __init__(self, a, b, theta, x, y, phi, mol):
        self.a = a
        self.b = b
        self.theta = theta
        self.mol = mol
        self.mols = []
        self.addMol(x,y,phi)
        self.addMol(x,y,phi)
        self.rotation(1,2,1)
        self.addMol(0.5-x, 0.5+y, phi)
        self.addMol(0.5+x, 0.5-y, phi+180)
        self.crys = "p2gg"

class p2(cell):
    def __init__(self, a, b, theta, x, y, phi, mol):
        self.a = a
        self.b = b
        self.theta = theta
        self.mol = mol
        self.mols = []
        self.addMol(x,y,phi)
        self.addMol(x,y,phi)
        self.rotation(1,2,1)
        self.crys = "p2"

def lammpsFile(cell,path='.'):
    mol = cell.getMol()
    a,b,theta = cell.getShape()
    xy = 0
    if theta != 0:
        xy = b*tan(pi/2-theta)
    string = ""
    string += "ITEM: TIMESTEP\n"
    string += "0\n"
    string += "ITEM: NUMBER OF ATOMS\n"
    string += "{}\n".format(cell.numAtoms())
    string += "ITEM: BOX BOUNDS\n"
    string += "{xmin} {xmax} {xy}\n".format(xmin=0, xmax=a+xy, xy=xy)
    string += "{ymin} {ymax} {yz}\n".format(ymin=0, ymax=b, yz=0)
    string += "{zmin} {zmax} {zx}\n".format(zmin=-0.5, zmax=0.5, zx=0)
    string += "ITEM: ATOMS\n"
    atomID = 1
    molID = 1
    for m in cell.getMols():
        for atom in m:
            x,y = atom.getPos()
            x = wrap(x, cell.getA())
            string += "{id} {molID} {type} {diam} {x} {y} {z}\n".format(\
                    id=atomID, molID=molID, type=atom.getType(),\
                    diam=2*atom.getSize(), x=x, y=y, z=0)
            atomID += 1
        molID += 1
    if mol.getAngles():
        filename="{shape}-{radius}-{dist}-{theta}".format(shape=mol.getName(),radius=mol.radius, dist=mol.dist, theta=mol.theta)
    else:
        filename="{shape}-{radius}-{dist}".format(shape=mol.getName(),radius=mol.radius, dist=mol.dist)
    f = open('{path}/{filename}.lammpstrj'.format(path=path, filename=filename),'w')
    f.write(string)
    f.close()

def wrap(x, a):
    x = 2*pi*x/a
    x = (atan2(sin(x),cos(x))+pi)*(a/(2*pi))
    return x

def cellFile(cell,path='.'):
    mol = cell.getMol()
    a,b,theta = cell.getShape()
    string = ""
    string += "# Unit cell data for {} molecule\n\n".format(mol.getName())
    string += '{}   atoms\n'.format(cell.numAtoms())
    string += '{}   atom types\n'.format(mol.numAtomTypes())
    string += '{}   extra bond per atom\n'.format(mol.maxBondCount())
    string += '{}   extra angle per atom\n\n'.format(mol.maxAngleCount())
    string += '0 {xdim}     xlo xhi\n'.format(xdim=a)
    string += '0 {ydim}     ylo yhi\n'.format(ydim=b)
    string += '{xy} {xz} {yz}   xy xz yz\n\n'.format(xy=0,xz=0,yz=0)
    
    string += '\nAtoms\n\n'
    string += str(cell)
    if mol.getBonds():
        string += '\nBond Coeffs\n\n'
        a,b,d = mol.getBonds()[0]
        string += '{bondID} {coeff} {dist}\n'.format(\
                bondID=1, coeff=500, dist=d)
    if mol.getAngles():
        string += '\nAngle Coeffs\n\n'
        a,b,c,theta = mol.getAngles()[0]
        string += '{angleID} {coeff} {theta}\n'.format(\
                angleID=1, coeff=5000, theta=theta*180/pi)
    string += '\nMasses\n\n'
    for t in mol.getAtomTypes():
        string += '{type} {mass}\n'.format(type=t.getType(), mass=1)
    string += '\nPair Coeffs\n\n'
    for t in mol.getAtomTypes():
        string += '{0} {strength} {dist}\n'.format(t.getType(), strength=1, dist=2*t.getSize())
    string += "\nAtoms\n"
    if mol.getAngles():
        filename="{shape}-{radius}-{dist}-{theta}".format(shape=mol.getName(),radius=mol.radius, dist=mol.dist, theta=mol.theta)
    elif cell.getCrys():
        filename="{shape}-{radius}-{dist}-{crys}".format(shape=mol.getName(),radius=mol.radius, dist=mol.dist, crys=cell.getCrys())
    else:
        filename="{shape}-{radius}-{dist}".format(shape=mol.getName(),radius=mol.radius, dist=mol.dist)
    f = open('{path}/{filename}.dat'.format(path=path, filename=filename),'w')
    f.write(string)
    f.close()

def molFile(molecule, path='.', crys = ""):
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
    elif crys:
        filename="{shape}-{radius}-{dist}-{crys}".format(shape=molecule.getName(),radius=molecule.radius, dist=molecule.dist, crys=crys)
    else:
        filename="{shape}-{radius}-{dist}".format(shape=molecule.getName(),radius=molecule.radius, dist=molecule.dist)

    f = open('{path}/{filename}.mol'.format(path=path,filename=filename), 'w')
    f.write(string)
    f.close()

if __name__ == "__main__":
    c = cell(3.43,3.13, 40.5, molecule.Trimer(0.7,1.0,180))
    c.addParticle(0.734,0.25,0)
    c.addParticle(1.177,0.75,0)
    print repr(c)
    c.replicate(2,2)
    cellFile(c)
    molFile(c.getMol())

