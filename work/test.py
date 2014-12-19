
from math import *

class cell:
    def __init__(self,a,b,theta,x,y,phi,na,nb):
        self.a = a
        self.b = b
        self.na = na
        self.nb = nb
        self.x = x
        self.y = y
        self.phi = phi
        self.theta = theta

    def __str__(self):
        na = self.na
        nb = self.nb
        a = self.a
        b = self.b
        theta = self.theta
        x = self.x
        y = self.y
        s = ''
        s = ''
        s += "{0},{1}\n".format(0,0)
        s += "{0},{1}\n\n".format(nb*b*cos(theta),nb*b*sin(theta))
        s += "{0},{1}\n".format(0,0)
        s += "{0},{1}\n\n".format(na*a,0)
        s += "{0},{1}\n".format(na*a,0)
        s += "{0},{1}\n\n".format(na*a + nb*b*cos(theta), nb*b*sin(theta))
        s += "{0},{1}\n".format(nb*b*cos(theta), nb*b*sin(theta))
        s += "{0},{1}\n\n".format(na*a+nb*b*cos(theta), nb*b*sin(theta))
        xi = 1-x
        yi = 1-y
        x = x*a + y*b*cos(theta)
        y = y*b*sin(theta)
        x0 = x
        y0 = y
        xi = xi*a + yi*b*cos(theta)
        yi = yi*b*sin(theta)
        d = 1.0
        r = 0.7
        const = (1**2 + d**2 - r**2)/(2*d*1)
        #x = x - const*cos(self.phi)
        #y = y - const*sin(self.phi)
        #xi = xi - const*cos(self.phi+pi)
        #yi = yi - const*sin(self.phi+pi)
        #xi = a+b*cos(theta)-x
        #yi = b*sin(theta)-y
        
        for i in xrange(na):
            for j in xrange(nb):
                xn = x + i*a + j*b*cos(theta)
                yn = y + j*b*sin(theta)
                xin = xi + i*a + j*b*cos(theta)
                yin = yi + j*b*sin(theta)
                x0n = x0 + i*a + j*b*cos(theta)
                y0n = y0 + j*b*sin(theta)
                s += "{0},{1},1\n".format(xn, yn)
                s += "{0},{1},{2}\n\n".format(xn+d*cos(self.phi), yn+d*sin(self.phi),r)
                s += "{0},{1},1\n".format(xin, yin)
                s += "{0},{1},{2}\n\n".format(xin+d*cos(self.phi+pi), yin+d*sin(self.phi+pi),r)
                s += "{0},{1},0.1\n\n".format(x0n, y0n)
        return s

# Perfect Packing p2
c = cell(3.9,3.183253,0.926517,0.205,0.1675,2*pi-2.5,4,4)

# d =1 r = 0.637556
c = cell(3.234679,2.675387,1.674524,0.772394,0.506326,2*pi-1.290112,4,4)

c = cell(4.152241,2.000398,1.570796,0.864893,0.704416,0.618550,1,1)

#c = cell(4,4,pi/3,0.2,0.2,0,4,4)

print c
