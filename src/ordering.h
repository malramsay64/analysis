//
//  mol_functions.h
//  analysis
//
//  Created by Malcolm Ramsay on 7/12/2014.
//  Copyright (c) 2014 Malcolm Ramsay. All rights reserved.
//

#ifndef MY_MOL_FUNCTIONS
#define MY_MOL_FUNCTIONS

#include <complex>
#include "Particle.h"
#include "Frame.h"
#include "neighbours.h"
#include "dyn_queue.h"

static double STRUCT_DIST = 0.3;

Vector<2> orientation(const Molecule &, const Frame &);
double angle(const Molecule &, const Frame &);
Vector<2> wrap_x(const Vector<2> &v, double a);
double com_colour(const Molecule &, const Frame &);
double mol_colour(const Molecule &, const Frame &);
double struct_relax(const Molecule &, const Frame &);
double hexatic(int n, const Molecule &, const Frame &);
double circle_ordering(const Molecule &);
int circle_colour(const Molecule &);
int neighbour_colour(const Molecule &, const Frame &);
double orient_ordering(const Molecule &);
int short_ordering(const Molecule &, const Frame &);

int tri_ordering(const Molecule &);
Molecule reorient(const Molecule &, const Frame &);

#endif /* defined(MY_MOL_FUNCTIONS) */
