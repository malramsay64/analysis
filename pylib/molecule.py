#!/usr/bin/python

import math

class Atom:
    """Creates a atom with size, mass and position"""
    def __init__(self, atomID, atomType, size, mass=-1):
        self.atomID = atomID
        self.atomType = atomType
        size = float(size)
        if size > 0:
            self.size = size
        else:
            raise Exception("The size of the particle must be positive")
        if mass > 0:
            self.mass = mass
        else:
            # default density = 1/pi
            self.mass = size**2
        self.pos = (0,0)

    def setPos(self, x, y):
        self.pos = (x,y)

    def getPos(self):
        return self.pos

    def getMass(self):
        return self.mass

    def getSize(self):
        return self.size
    
    def getType(self):
        return self.atomType

    def getID(self):
        return self.atomID

    def area(self):
        return math.pi*size**2

class Molecule:
    """Creates a molecule, a collection of atoms"""
    def __init__(self):
        atoms = []
        bonds = []
        angles = []

    def getAtoms(self):
        return self.atoms
    
    def __iter__(self):
        return iter(self.atoms)
    
    def getName(self):
        return self.name

    def getAngles(self):
        return self.angles

    def getBonds(self):
        return self.bonds

    def numBonds(self):
        return len(self.getBonds())

    def numAngles(self):
        return len(self.getAngles())

    def maxBondCount(self):
        m = 0
        for atom in self:
            count = self.bondCount(atom.getID())
            if count > m:
                m = count
        return m
    
    def maxAngleCount(self):
        m = 0
        for atom in self:
            count = self.angleCount(atom.getID())
            if count > m:
                m = count
        return m


    def bondCount(self, ID):
        count = 0
        for a,b,dist in self.getBonds():
            if a == ID or b == ID:
                count += 1
        return count

    def getBondLength(self, atom1, atom2):
        for a,b,dist in self.getBonds():
            if atom1 in [a,b] and atom2 in [a,b]:
                return dist
        return 0

    def get12(self, atomID):
        bonded = []
        for a,b,dist in self.getBonds():
            if a == atomID:
                bonded.append(b)
            elif b == atomID:
                bonded.append(a)
        return bonded

    def get13(self, atomID):
        bonded = []
        for a,b,c,theta in self.getAngles():
            if a == atomID:
                bonded.append(c)
            elif c == atomID:
                bonded.append(a)
        return bonded

    def angleCount(self, atomID):
        count = 0
        for a,b,c,theta in self.getAngles():
            if a == atomID or c == atomID:
                count += 1
        return count

    def getAtomTypes(self):
        seen = set()
        types = []
        for atom in self:
            if atom.getType() not in seen:
                types.append(atom)
                seen.add(atom.getType())
        return types

    def numAtoms(self):
        return len(self.getAtoms())

    def numAtomTypes(self):
        seen = set()
        [seen.add(atom.getType()) for atom in self]
        return len(seen)

    def getAtom(self, name):
        for atom in self:
            if atom.getID() == name:
                return atom
        return 0

    def posRange(self):
        xmin = 0
        ymin = 0
        xmax = 0
        ymax = 0
        for atom in self:
            x,y = atom.getPos()
            if x < xmin:
                xmin = x
            if x > xmax:
                xmax = x
            if y < ymin:
                ymin = y
            if y > ymax:
                maxy = y
        return xmin, xmax, ymin, ymax

    def dimensions(self):
        xmin, xmax, ymin,ymax = self.valRange()
        return xmax - xmin, ymax - ymin

    def valRange(self):
        xmin = 0
        xmax = 0
        ymin = 0
        ymax = 0
        for atom in self:
            r = atom.getSize()
            x,y = atom.getPos()
            if x-r < xmin:
                xmin = x-r
            if x+r > xmax:
                xmax = x+r
            if y-r < ymin:
                ymin = y-r
            if y+r > ymin:
                ymax = y+r
        return xmin, xmax, ymin, ymax

    def area(self):
        return sum([a.area() for a in self])  

    def COM(self):
        x = 0.
        y = 0.
        mass = 0
        for atom in self:
            dx, dy = atom.getPos()
            x += dx*atom.getMass()
            y += dy*atom.getMass()
            mass += atom.getMass()
        return (x/float(mass), y/float(mass))

    def orientation(self):
        x1,y1 = self.COM()
        x2,y2 = self.getAtom(1).getPos()
        v1 = x1-x2
        v2 = y1-y2
        l = math.sqrt(v1*v1 + v2*v2)
        if self.numAtoms() == 2 or l < 0.05:
            x2,y2 = self.getAtom(2).getPos()
            v2 = x1-x2
            v1 = -(y1-y2)
            l = math.sqrt(v1*v1 + v2*v2)
        return (v1 / l, v2/l)

    def angle(self):
        v = self.orientation()
        return math.atan2(v[0], v[1])*180/math.pi + 90


    def setAngle(self, theta):
        """ The molecule is rotated around the first/largest molecule"""
        dTheta = (self.angle() - theta)*math.pi/180
        x1,y1 = self.getAtom(1).getPos()
        for atom in self.getAtoms():
            if atom.getID() != 1:
                dist = self.getBondLength(1,atom.getID())
                x,y = atom.getPos()
                bondAngle = math.atan2(x-x1,y-y1)
                x = x1+dist*math.sin(bondAngle + dTheta)
                y = y1+dist*math.cos(bondAngle + dTheta)
                atom.setPos(x,y)


    def translate(self, dx, dy):
        for atom in self.getAtoms():
            x,y = atom.getPos()
            atom.setPos(x+dx, y+dy)

    def rotate(self, dtheta):
        self.setAngle(self.angle()+dtheta)

    def setPos(self, x, y):
        dx, dy = self.getAtoms()[0].getPos()
        dx = x - dx
        dy = y - dy
        self.translate(dx,dy)


    def __str__(self):
        str = ''
        for atom in self.getAtoms():
            x,y = atom.getPos()
            str += "{x} {y} {size}\n".format(x=x,y=y,size=atom.getSize())
        return str

class Snowman(Molecule):
    """Creates a snowman molecule with a single degree of freedom
       in the size of the second atom. The first atom has size 1 
       and is placed at the origin. The other has its center
       at a distance 1, on the radius of the first circle """
    def __init__(self, radius=1, dist=1):
        self.radius = radius
        self.dist = dist
        radius = float(radius)
        dist = float(dist)
        self.atoms = []
        self.bonds = []
        self.angles = []
        self.atoms.append(Atom(1,1,1))
        if radius == 1:
            self.atoms.append(Atom(2,1,1))
        else:
            self.atoms.append(Atom(2,2,radius))
        self.atoms[1].setPos(0,dist)
        self.bonds.append((1,2,dist))
        self.name = "Snowman"
    
    def area(self):
        return math.pi*(sum([a.size()**2 for a in self]))

class Trimer(Molecule):
    """ Creates a trimer with 2 degrees of freedom. The first atom
        has size 1 and placed at the origin. The remaining 2 atoms
        of size radius are placed on the circumference of the first
        circle with an angle of theta (degrees) separating them."""
    def __init__(self, radius=1, dist=1, theta=120, ratio=1):
        self.radius = radius
        self.dist = dist
        self.theta = theta
        radius = float(radius)
        dist = float(dist)
        theta = float(theta)
        ratio = float(ratio)
        self.atoms = []
        self.bonds = []
        self.angles = []
        self.atoms.append(Atom(1,1,1))
        theta *= (math.pi/180) # convert to radians
        if ratio == 1 and radius == 1:
            self.atoms.append(Atom(2,1,1))
            self.atoms.append(Atom(3,1,1))
        elif ratio == 1:
            self.atoms.append(Atom(2,2,radius))
            self.atoms.append(Atom(3,2,radius))
        else:
            self.atoms.append(Atom(2,2,radius))
            self.atoms.append(Atom(3,3,radius*ratio))
        self.atoms[1].setPos(dist*math.cos(theta/2), dist*math.sin(theta/2))
        self.atoms[2].setPos(dist*math.cos(theta/2), dist*-math.sin(theta/2))
        self.angles.append((2,1,3,theta))
        self.bonds.append((1,2,dist))
        self.bonds.append((1,3,dist))
        self.name = "Trimer"

    def area(self):
        d2 = math.sqrt(sum([d[-1]**2 for d in self.getBonds()]) + reduce(mul, [x[-1] for x in self.getBonds()])*math.cos(self.getAngles()[-1]))
        return math.pi*([ d[-1]**2 for d in self.getBonds()])

class Disc(Molecule):
    """ Creates a single disc as a molecule, with radius 1"""
    def __init__(self):
        self.atoms = []
        self.bonds = []
        self.angles = []
        self.atoms.append(Atom(1,1,1))
        self.name = "Disc"

def combinations(iterable, r):
    pool = tuple(iterable)
    n = len(pool)
    for indices in permutations(range(n), r):
        if sorted(indices) == list(indices):
            yield tuple(pool[i] for i in indices)

def overlap(r1,r2,d):
    if (d >= r1+r2):
        return 0
    try:
        height = math.sqrt(r1**2 - (dist**2 - r2**2 + r1**2)/(2*dist))
    except:
        return 0
    theta1 = 2*math.asin(height/r1)
    theta2 = 2*math.asin(height/r2)
    sec1 = 0.5*theta1*r1**2
    sec2 = 0.5*theta2*r2**2
    quad = dist*height
    return sec1 + sec2 - quad

