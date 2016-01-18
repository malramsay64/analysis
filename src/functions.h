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

#include <math.h>
#include <iostream>
#include <string>
#include <vector>

#ifndef XX
#define XX 0
#endif
#ifndef YY
#define YY 1
#endif
#ifndef PI
#define PI M_PI
#endif

#ifndef FUNCTIONS_H
#define FUNCTIONS_H

double dot_product(double *v1, double *v2, int len);
double legendre(int l, double x);
double dist(double, double);
double my_mod(double, double);
double structure_factor(double q, std::vector<double> g, double rho, double dr);
double max_structure_factor(std::vector<double> g, double rho, double dr);
int pos_def_mod(int, int);

#endif /* defined(FUNCTIONS_H) */
