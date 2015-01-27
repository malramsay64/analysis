//
//  particle.h
//  analysis
//
//  Created by Malcolm Ramsay on 7/12/2014.
//  Copyright (c) 2014 Malcolm Ramsay. All rights reserved.
//

#ifndef PARTICLE_H
#define PARTICLE_H

#include <iostream>
#include <vector>
#include <list>
#include <algorithm>
#include <math.h>
#include <fstream>
#include "functions.h"
#include <pthread.h>
#include "vect.h"
#include "constants.h"

// Particle class
class particle {
public:
    vect pos;
    int id;
    int molid;
    int type;
    double radius;
    double mass;
    std::vector<particle *> my_neighbours;
    
    particle();
    void append(particle *);
    int numn();
    int set_xpos(double);
    int set_ypos(double);
    int set_pos(vect);
    vect pos_vect();
    int index();
    int m_i();
    
    bool operator> (const particle &b);
    bool operator>= (const particle &b);
    
    bool operator< (const particle &b);
    bool operator<= (const particle &b);
    
    bool operator== (const particle &b);
    bool operator!= (const particle &b);
};

bool operator < (const particle &a, const particle &b);


#endif

