//
//  my_mean.h
//  analysis
//
//  Created by Malcolm Ramsay on 7/12/2014.
//  Copyright (c) 2014 Malcolm Ramsay. All rights reserved.
//

#ifndef MY_MEAN
#define MY_MEAN

#include "math.h"

class my_mean{
    int n;
    double mean;
    double M2;
    double weight;
public:
    my_mean();
    double add(double);
    double add(double, double);
    double add(my_mean);
    double get_mean();
    double get_stdev();
    double get_variance();
    double combine(my_mean);
    int get_count();
    
};

#endif /* defined(MY_MEAN) */
