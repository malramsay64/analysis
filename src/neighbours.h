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
#include <Stats.h>

#include "functions.h"
#include "Frame.h"
#include "ordering.h"
#include "dyn_queue.h"

//extern static int short_order_types;

int short_range_order(const Frame &frame);
bool find_mol_neighbours(Molecule &mol, Frame &frame, std::vector<std::vector<int>> &neigh_list);
int check_particles(Molecule & mol1, Molecule &mol2, Frame &frame);
int recompute_neighbours(Molecule &mol, Frame &frame, std::vector<std::vector<int>> &neigh_list);
int order_type(const Molecule &, const Molecule &, const Frame &);
void add_mol_neighbours(Molecule &m1, Molecule &m2);
void add_part_neighbours(Particle &p1, Particle &p2);
std::vector<int> short_neighbour_list(const Molecule &m, const Frame &frame);
int find_neighbours(Frame &frame, std::vector<std::vector<int>> &neigh_list);


#endif /* endif(NEIGHBOURS_H) */

