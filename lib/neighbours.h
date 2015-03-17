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

//extern static int short_order_types;

int short_range_order(Frame * frame);
bool find_mol_neighbours(molecule * mol, Frame * frame, std::vector<std::vector<int>> *neigh_list);
int check_particles(molecule * mol1, molecule * mol2, Frame * frame);
int recompute_neighbours(molecule * mol, Frame * frame, std::vector<std::vector<int>> *neigh_list);
//int check_single(Frame *, std::list<particle *> *, dyn_queue *,int, std::ofstream *);
int order_type(molecule * m1, molecule * m2, Frame * frame);
int add_mol_neighbours(molecule * m1, molecule *m2);
int add_part_neighbours(particle *p1, particle *p2);
std::vector<int> short_neighbour_list(molecule * m, Frame * frame);


#endif /* endif(NEIGHBOURS_H) */

