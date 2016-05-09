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

//extern static int short_order_types;

int short_range_order(Frame * frame);
bool find_mol_neighbours(Molecule * mol, Frame * frame, std::vector<std::vector<int>> *neigh_list);
int check_particles(Molecule * mol1, Molecule * mol2, Frame * frame);
int recompute_neighbours(Molecule * mol, Frame * frame, std::vector<std::vector<int>> *neigh_list);
//int check_single(Frame *, std::list<Particle *> *, dyn_queue *,int, std::ofstream *);
int order_type(Molecule * m1, Molecule * m2, Frame * frame);
int add_mol_neighbours(Molecule * m1, Molecule *m2);
int add_part_neighbours(Particle *p1, Particle *p2);
std::vector<int> short_neighbour_list(Molecule * m, Frame * frame);
int find_neighbours(Frame *frame, std::vector<std::vector<int>> *neigh_list);


#endif /* endif(NEIGHBOURS_H) */

