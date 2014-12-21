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
#include "angle_list.h"
#include "dyn_queue.h"

int get_colour(particle *p, Frame *frame);
int set_colour(Frame *frame);
int graph_colour(particle *p, Frame *frame);
double average_bonded_frac(particle *p, Frame * frame);
vect orientation(molecule *, Frame *);
double angle(molecule *, Frame *);
int colourise(Frame *frame);
angle_list like_me(Frame *frame);
double local_order(molecule * m, Frame * frame);
double global_order(molecule *m, Frame * frame);
double circle_order(particle * p, Frame * frame);
std::ostream& print_mol(std::ostream &, molecule *,Frame *);


#endif /* defined(MY_MOL_FUNCTIONS) */

