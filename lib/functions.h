//
//  functions.h
//  analysis
//
//  Created by Malcolm Ramsay on 7/12/2014.
//  Copyright (c) 2014 Malcolm Ramsay. All rights reserved.
//

#include <math.h>
#include <iostream>
#include <string>

#ifndef XX
#define XX 0
#endif
#ifndef YY
#define YY 1
#endif
#ifndef PI
#define PI 3.14159265358979
#endif

#ifndef FUNCTIONS_H
#define FUNCTIONS_H

double dot_product(double *v1, double *v2, int len);
double legendre(int l, double x);
double dist(double, double);
double my_mod(double, double);

#endif /* defined(FUNCTIONS_H) */
