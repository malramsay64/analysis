#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2014 malcolm <malcolm@asaph-VirtualBox>
#
# Distributed under terms of the MIT license.

"""

"""
import unitCell
import molecule
from math import *
import sys

def create(a, b, theta, x, y, phi, molecule, crys, mols=2500, path='.'):
    s = crys(a,b,theta, x, y, phi, molecule)
    if mols:
        mols /= s.numMols()
        na = sqrt(mols*b/a)
        nb = mols/na
        na = int(na)
        nb = int(nb)
        s.replicate(na, nb)
    print repr(s)
    unitCell.cellFile(s, path)
    unitCell.molFile(s.getMol(), path, s.getCrys())
    unitCell.lammpsFile(s, path)

if __name__ == "__main__":
    args = sys.argv
    line = sys.stdin.readline().split()
    path = '.'
    mols = 2500

    if len(args) == 6 and len(line) == 7:
        theta = pi- float(line[0])
        a = float(line[1])
        b = float(line[2])
        x = float(line[3])
        y = float(line[4])
        phi = float(line[5])
        m = float(line[6])
        
        path = args[1]
        mols = int(args[2])
        r = float(args[3])
        d = float(args[4])
        crys = getattr(unitCell,args[5])

        s = molecule.Snowman(r,d)
       
        # Convert to my coordinates
        ab = -b/tan(theta)
        #a += ab
        
        #phi += theta
        #phi = pi-phi
        x = a*x + b*cos(theta)
        y = b*sin(theta)
       

        conv = -(1**2 +d**2 - r**2)/(2*d)
        x += conv*cos(phi)
        y += conv*sin(phi)
        
        #phi += pi/2
        #phi = phi*180/pi
        #b *= sin(theta)
        #ab = b*sin(theta)
        create(a, b, theta, x/a, y/b, phi, s, crys, mols, path)
         
