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
#include <algorithm>
#include <fstream>
#include "functions.h"
#include "Vector.h"

using namespace LAlgebra;

struct particle_vars{
    double xpos = 0;
    double ypos = 0;
    int id = 0;
    int molid = 0;
    int type = 0;
    double radius = 1;
    double mass = 1;
};

// particle class
class Particle {
public:
    Vector<2> pos;
    int id;
    int molid;
    int type;
    double radius;
    double mass;
    std::vector<Particle *> my_neighbours;
    
    Particle();
    Particle(const Particle &);
    Particle(const particle_vars&);

    void append(Particle *);
    size_t numn() const;
    void set_pos(const Vector<2> &);
    Vector<2> pos_vect() const;
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

