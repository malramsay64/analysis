//
//  my_mean.cpp
//  analysis
//
//  Created by Malcolm Ramsay on 7/12/2014.
//  Copyright (c) 2014 Malcolm Ramsay. All rights reserved.
//

#include "my_mean.h"

my_mean::my_mean(){
    m = 0;
    stdev = 0;
    n = 0;
}

double my_mean::add(double val){
    double mprev = m;
    n += 1;
    m += (val - m)/(n);
    stdev += (val- m) * (val - mprev);
    return m;
}

double my_mean::add(my_mean val){
    m = (n*m + val.n*val.m)/(n+val.n);
    n += val.n;
    return m;
}

double my_mean::get_mean(){
    return m;
}

double my_mean::get_stdev(){
    return stdev;
}

double my_mean::combine(my_mean val){
    m = (n*m + val.n*val.m)/(n+val.n);
    n += val.n;
    return m;
}

int my_mean::get_count(){
    return n;
}