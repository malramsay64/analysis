//
//  mol_functions.h
//  analysis
//
//  Created by Malcolm Ramsay on 7/12/2014.
//  Copyright (c) 2014 Malcolm Ramsay. All rights reserved.
//

#ifndef MY_MOL_FUNCTIONS
#define MY_MOL_FUNCTIONS

#include "particle.h"
#include "frame.h"
#include "my_mean.h"
#include "neighbours.h"
#include "dyn_queue.h"

static double STRUCT_DIST = 0.3;

vect orientation(molecule *, Frame *);
double angle(molecule *, Frame *);
vect wrap_x(vect v, double a);
int mol_colour(molecule * m, Frame * frame);
double struct_relax(molecule * m, Frame * frame);

#endif /* defined(MY_MOL_FUNCTIONS) */

