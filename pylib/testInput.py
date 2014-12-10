#!/usr/bin/python

import molecule
import unitCell
import sys

prefix=sys.argv[1]


m = molecule.Snowman(0.637556, 1.637556)
positions = [(0.399, 0.559, 203.5), (1.00, 0.559, 23.5)]
s = unitCell.cell(4.74, 1.79, 37, m)
for x,y,t, in positions:
    s.addParticle(x,y,t)

s.replicate(5,5)
unitCell.lammpsFile(s,prefix)
