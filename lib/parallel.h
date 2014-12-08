//
//  parallel.h
//  analysis
//
//  Created by Malcolm Ramsay on 7/12/2014.
//  Copyright (c) 2014 Malcolm Ramsay. All rights reserved.
//

#ifndef MY_PARALLEL
#define MY_PARALLEL

#include <thread>

#include "vect.h"
#include "particle.h"
#include "frame.h"
#include "mol_functions.h"
#include "neighbours.h"
#include "my_mean.h"

template <class type> void * mean_single_thread(std::vector<type> *l1, int begin, int end, double (*f)(type *, Frame *), Frame *, my_mean *res);
template <class type> void * mean_thread(std::vector<type> *l1, std::vector<type> *l2, int begin, int end, vect (*f)(type *, Frame *), double (*diff)(vect, vect, int), Frame *, int l, my_mean *res);

template <class type> my_mean par_mean(std::vector<type> *list, double (*f)(type *, Frame *), Frame *);

double v_legendre(vect v1, vect v2, int l);
double v_dist(vect v1, vect v2, int l);
vect COM(molecule *, Frame *);
double par_rotational_relaxation(int l, Frame *init, Frame *curr);
double par_average_moved(Frame *init, Frame *curr);
double par_MSD(Frame *init, Frame *curr);
double par_diffusion_constant(Frame *init, Frame *curr);
double bonded_frac(Frame *frame);
double par_average_rotation(Frame *init, Frame *curr);
int order_parameter(int reference, std::ofstream *file, Frame *frame);
double par_local_order(Frame * frame);
double par_global_order(Frame * frame);
double par_circle_order(Frame * frame);

#endif

