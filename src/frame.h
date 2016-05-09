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
#include <sstream>

#include "Particle.h"
#include "Molecule.h"

#ifndef FRAME_H
#define FRAME_H

#define DELTA_T 0.005

class Frame {
    int atoms;
    bool coloured;
    bool neighbours;
    double a;
    double b;
    double theta;
    
public:
    std::vector<Particle> particles;
    std::vector<Molecule> molecules;
    int timestep;
    int num_molecules;
    
    Frame();
    Frame(Frame const &);
    void set_timestep(int);
    int get_timestep();
    int get_time();
    void set_atoms(int);
    int num_atoms();
    void set_num_mol(int);
    int num_mol();
    void add_particle(Particle *p);
    void add_link(int, int, bool);
    void add_link(int, int);
    double dist(Vector2d,Vector2d);
    Vector2d direction(Vector2d,Vector2d);
    double get_area();
    double get_density();
    Molecule at(int);
    void update_links();
    
    int set_crys(double, double, double);
    double get_a();
    double get_b();
    double get_theta();
    double get_height();
    double get_tilt();
    Vector2d cartesian(Vector2d v);
    Vector2d fractional(Vector2d v);
};

#endif /* defined(FRAME_H) */
