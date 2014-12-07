#include "functions.h"
#include "frame.h"
//#include "parallel.h"
#include "mol_functions.h"
#include <cmath>
#include <queue>
#include <thread>
#include <iostream>
//#include "output.h"

#ifndef NEIGHBOURS_H
#define NEIGHBOURS_H 

static int NUM_THREADS = 4;

int par_neigh(Frame *frame);
void *loop_neigh(Frame * frame, int begin, int end);
void find_neighbours(particle *a, Frame *frame);
int short_range_order(Frame * frame);
void randomise_orientation(Frame * frame);
int remove_neighbours(particle *, std::list<particle*>*); 
int check_single(Frame *, std::list<particle *> *, dyn_queue<particle> *,int, std::ofstream *);
int order_type(molecule * m1, molecule * m2, Frame * frame);
#endif
