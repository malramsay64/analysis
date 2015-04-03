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
#include "particle.h"
#include "frame.h"
#include "my_mean.h"
#include "neighbours.h"
#include "dyn_queue.h"

static double STRUCT_DIST = 0.3;

vect orientation(molecule *, Frame *);
double angle(molecule *, Frame *);
vect wrap_x(vect v, double a);
double com_colour(molecule *, Frame *);
double mol_colour(molecule * m, Frame * frame);
double struct_relax(molecule * m, Frame * frame);
double hexatic(int n, molecule* m1, Frame *frame);
double circle_ordering(molecule *m);
int circle_colour(molecule * m);
int neighbour_colour(molecule * m, Frame *frame);
double orient_ordering(molecule *m);
int short_ordering(molecule *m, Frame * frame);
molecule reorient(molecule *m, Frame *frame);

#endif /* defined(MY_MOL_FUNCTIONS) */
