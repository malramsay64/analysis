//
//  vect.h
//  analysis
//
//  Created by Malcolm Ramsay on 7/12/2014.
//  Copyright (c) 2014 Malcolm Ramsay. All rights reserved.
//

#include "functions.h"

#ifndef MY_VECT
#define MY_VECT

class vect {
public:
    double x;
    double y;
    
    vect();
    vect(double, double);
    vect(double *v);
    void normalise();
    double length();
    double orientation();
    void orthogonalise();
    
    vect operator+= (double);
    vect operator+= (vect);
    vect operator+= (int);
    vect operator-();
    vect operator/= (double i);
    vect operator+ (double i);
    vect operator+ (const vect &v);
    vect operator- (const vect &v);
    vect operator- (double i);
    vect operator+ (int i);
    vect operator/ (double i);
    vect operator* (double i);
    vect operator* (const vect v);
    vect operator/ (const vect v);
    friend std::ostream& operator<<(std::ostream &os, const vect &v){
        os << v.x << " " << v.y;
        return os;
    }
    
};


vect operator* (double i, vect v);

double dist(vect v1, vect v2);
vect direction(vect v1, vect v2);
double dot_product(vect v1, vect v2);
vect cos(vect v);
vect sin(vect v);
vect atan2(vect vx, vect vy);
double atan2(vect v);
double map_angle(double x);

#endif

