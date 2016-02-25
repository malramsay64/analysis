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

/* TODO
 * Create intialiser stuct for Frame
 * Struct collated input data to create the frame
 * Generalise to 3D
 */

struct FrameVars {
    double a{0};
    double b{0};
    double theta{0};
    double step_size = 0.005;
    size_t num_atoms{0};
    size_t num_mols{0};
    int timestep{0};
};

class Frame {
    bool coloured;
    bool neighbours;
    double a;
    double b;
    double theta;
    double step_size = 0.005;
    
public:
    std::vector<Particle> particles;
    std::vector<Molecule> molecules;
    int timestep;
    // TODO particle/molecule binning

    Frame();
    Frame(const FrameVars &f);
    Frame(Frame const &);
    void set_timestep(int);
    int get_timestep() const;
    double get_step_size() const;
    int get_time() const;
    void set_atoms(int);
    size_t num_atoms() const;
    void set_num_mol(int);
    size_t num_mol() const;
    void add_particle(const Particle &p);
    void add_link(int, int, bool);
    void add_link(int, int);
    double get_area() const;
    double get_density() const;
    Molecule at(int) const;
    void update_links();
    void set_crys(double, double, double);
    double get_a() const;
    double get_b() const;
    double get_theta() const;
    double get_height() const;
    double get_tilt() const;
};

double dist(const Vector &, const Vector &, const Frame &);
Vector direction(const Vector &, const Vector &, const Frame&);
double angle(const Vector &, const Frame&);
Vector cartesian(const Vector &, const Frame&);
Vector fractional(const Vector &, const Frame&);

#endif /* defined(FRAME_H) */
