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
#include "frame.h"
#include "my_mean.h"
#include "neighbours.h"
#include "dyn_queue.h"

static double STRUCT_DIST = 0.3;

Vector2d orientation(Molecule *, Frame *);
double angle(Molecule *, Frame *);
Vector2d wrap_x(Vector2d v, double a);
double com_colour(Molecule *, Frame *);
double mol_colour(Molecule * m, Frame * frame);
double struct_relax(Molecule * m, Frame * frame);
double hexatic(int n, Molecule* m1, Frame *frame);
double circle_ordering(Molecule *m);
int circle_colour(Molecule * m);
int neighbour_colour(Molecule * m, Frame *frame);
double orient_ordering(Molecule *m);
int short_ordering(Molecule *m, Frame * frame);

int tri_ordering(Molecule *mol);
Molecule reorient(Molecule *m, Frame *frame);

#endif /* defined(MY_MOL_FUNCTIONS) */
