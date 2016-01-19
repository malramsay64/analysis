//
//  Particle.h
//  analysis
//
//  Created by Malcolm Ramsay on 7/12/2014.
//  Copyright (c) 2014 Malcolm Ramsay. All rights reserved.
//

#ifndef PARTICLE_H
#define PARTICLE_H

#include <iostream>
#include <vector>
#include <algorithm>
#include <fstream>
#include "functions.h"
#include "Vector2d.h"

struct particle_vars{
    double xpos = 0;
    double ypos = 0;
    int id = 0;
    int molid = 0;
    int type = 0;
    double radius = 1;
    double mass = 1;
};

// Particle class
class Particle {
public:
    Vector2d pos;
    int id;
    int molid;
    int type;
    double radius;
    double mass;
    std::vector<Particle *> my_neighbours;
    
    Particle();
    Particle(const Particle&);
    Particle(const particle_vars&);

    void append(Particle *);
    size_t numn() const;
    int set_pos(Vector2d);
    Vector2d pos_vect() const;
    int index() const;
    int mol_index() const;
    int order();
    void delete_neighbours();
    
    bool operator> (const Particle &b) const;
    bool operator>= (const Particle &b) const;
    
    bool operator< (const Particle &b) const;
    bool operator<= (const Particle &b) const;
    
    bool operator== (const Particle &b) const;
    bool operator!= (const Particle &b) const;
};



#endif

