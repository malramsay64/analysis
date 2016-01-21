//
//  functions.h
//  analysis
//
//  This file defines the basic functions needed for analysis that are
//  not part of the standard libraries
//
//  Created by Malcolm Ramsay on 7/12/2014.
//  Copyright (c) 2014 Malcolm Ramsay. All rights reserved.
//

#include <cmath>
#include <iostream>
#include <string>
#include <vector>



#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#define EPS 5e-2

#ifndef PI
#define PI M_PI
#endif

double dot_product(double *v1, double *v2, int len);
double legendre(int l, double x);
double dist(double, double);
double my_mod(double, double);
double structure_factor(double q, std::vector<double> g, double rho, double dr);
double max_structure_factor(std::vector<double> g, double rho, double dr);
int pos_def_mod(int, int);
double map_angle(double x);

#endif /* defined(FUNCTIONS_H) */
