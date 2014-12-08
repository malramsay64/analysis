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
    bool traversed;
public:
    int id;
    int molid;
    int type;
    double radius;
    double mass;
    double pos[2];
    std::vector<particle *> my_neighbours;
    
    particle();
    //~particle();
    void append(particle *);
    int numn();
    int n_large();
    int n_small();
    double xpos();
    double ypos();
    int set_xpos(double);
    int set_ypos(double);
    int set_pos(vect);
    vect pos_vect();
    int index();
    int m_i();
    bool get_traversed();
    void traverse();
    void reset_traverse();
    
    bool operator> (const particle &b);
    bool operator>= (const particle &b);
    
    bool operator< (const particle &b);
    bool operator<= (const particle &b);
    
    bool operator== (const particle &b);
    bool operator!= (const particle &b);
};

bool operator < (const particle &a, const particle &b);


#endif

