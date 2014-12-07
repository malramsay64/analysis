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

std::string split(std::string s, char delim);
std::string basename(std::string s);
double dot_product(double *v1, double *v2, int len);
double legendre(int l, double x);
double dist(double, double);

#endif
