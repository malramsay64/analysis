//
//  neighbours.h
//  analysis
//
//  Created by Malcolm Ramsay on 7/12/2014.
//  Copyright (c) 2014 Malcolm Ramsay. All rights reserved.
//

#ifndef NEIGHBOURS_H
#define NEIGHBOURS_H

#include <cmath>
#include <queue>
#include <thread>
#include <iostream>

#include "functions.h"
#include "frame.h"
#include "mol_functions.h"
#include "dyn_queue.h"
#include "constants.h"

int short_range_order(Frame * frame);
int check_single(Frame *, std::list<particle *> *, dyn_queue<particle> *,int, std::ofstream *);
int order_type(molecule * m1, molecule * m2, Frame * frame);
int add_mol_neighbours(molecule * m1, molecule *m2);
int add_part_neighbours(particle *p1, particle *p2);

#endif /* endif(NEIGHBOURS_H) */

