#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2014 malcolm <malcolm@macbook.local>

import sys
import math

class thermo():

    def __init__(self, filename):
        self.data = []
        f = open(filename, 'r')
        l = f.readline()
        while l:
            if l.startswith('Step'):
                if len(self.data) == 0:
                    [self.data.append([prop]) for prop in l.split()]
                l = f.readline()
                while not l.startswith('Loop'):
                    d = l.split()
                    for i in xrange(len(d)):
                        self.data[i].append(float(d[i]))
                    l = f.readline()
            l = f.readline()

    def last(self, prop):
        d = [x[-1] for x in self.data if x[0] == prop]
        return d[0]

    def output(self, fname, *args):
        f = open(fname, 'w')
        px = [x[0] for x in self.data].index("Step")
        py = []
        for prop in args:
            py.append([x[0] for x in self.data].index(prop))
        for i in xrange(len(self.data[px])):
            f.write("{0},{1}\n".format(self.data[px][i],",".join([str(self.data[h][i]) for h in py])))

        

def area(r, dist):
    return math.pi*(1+r**2) - overlap(1, r, dist)
    
def area3(r, dist, theta):
    d2 = math.sqrt(2*dist**2 + 2*dist*math.cos(theta))
    return math.pi*(1+2*r**2) - 2*overlap(1,r,dist) - overlap(r,r,d2)

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
         

def coverage(r,dist,density):
    a = area(r, dist)
    a /= 2 
    return density*a
    
def coverage3(r,dist,theta,density):
    a = area3(r,dist,theta)
    a /= 3
    return density*a


if __name__ == "__main__" :
    if len(sys.argv) > 1:
        filename  = sys.argv[1]
    else:
        filename = '1.log'
    f = thermo(filename)
    print 'Density: {0:.4f}'.format(f.last('Density'))
    print 'Total Energy: {0:.4f}'.format(f.last('TotEng'))
    if len(sys.argv) == 4:
        r = float(sys.argv[2])
        dist = float(sys.argv[3])
        print 'Packing Fraction: {0:.4f}'.format(coverage(r, dist, float(f.last('Density'))))
    elif len(sys.argv) == 5:
        r = float(sys.argv[2])
        dist = float(sys.argv[3])
        try:
            theta = float(sys.argv[4])
            print 'Packing Fraction: {0:.4f}'.format(coverage3(r, dist, theta, float(f.last('Density'))))
        except ValueError:
            print 'Packing Fraction: {0:.4f}'.format(coverage(r, dist, float(f.last('Density'))))
    f.output('energy.csv', "TotEng") 
    
