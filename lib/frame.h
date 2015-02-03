//
//  frame.h
//  analysis
//
//  Created by Malcolm Ramsay on 7/12/2014.
//  Copyright (c) 2014 Malcolm Ramsay. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <string>
#include <exception>
#include <string.h>
#include <stdlib.h>
#include <cmath>

#include "particle.h"
#include "molecule.h"

#ifndef FRAME_H
#define FRAME_H

#define DELTA_T 0.005

class Frame {
    int atoms;
    bool coloured;
    bool neighbours;
    double xdim[2];
    double ydim[2];
    double zdim[2];
    double a;
    double b;
    double theta;
    
public:
    std::vector<particle> particles;
    std::vector<molecule> molecules;
    int timestep;
    int num_molecules;
    std::vector<std::mutex> m;
    
    Frame();
    
    void set_timestep(int);
    int get_timestep();
    void set_atoms(int);
    int num_atoms();
    void set_num_mol(int);
    int num_mol();
    void add_particle(particle *p);
    void add_link(int, int, bool);
    void add_link(int, int);
    double dist(vect,vect);
    vect direction(vect,vect);
    double get_area();
    double get_density();
    
    int set_crys(double, double, double);
    double get_a();
    double get_b();
    double get_theta();
    double get_height();
    double get_tilt();
    vect cartesian(vect v);
    vect fractional(vect v);
    
    double xmin();
    double xmax();
    double ymin();
    double ymax();
    double zmin();
    double zmax();
};

#endif /* defined(FRAME_H) */
