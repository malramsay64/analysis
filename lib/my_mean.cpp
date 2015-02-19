//
//  my_mean.cpp
//  analysis
//
//  Created by Malcolm Ramsay on 7/12/2014.
//  Copyright (c) 2014 Malcolm Ramsay. All rights reserved.
//

#include "my_mean.h"

my_mean::my_mean(){
    mean = 0;
    M2 = 0;
    n = 0;
}

double my_mean::add(double x){
    n += 1;
    double delta = x - mean;
    mean += delta/(n);
    M2 += delta * (x - mean);
    return mean;
}

double my_mean::add(my_mean val){
    mean = (n*mean + val.n*val.mean)/(n+val.n);
    n += val.n;
    return mean;
}

double my_mean::get_mean(){
    return mean;
}

double my_mean::get_stdev(){
    if (n < 2){
        return 0;
    }
    return sqrt(M2/(n));
}

double my_mean::get_variance(){
    if (n < 2){
        return 0;
    }
    return M2/(n);

}

double my_mean::combine(my_mean val){
    mean = (n*mean + val.n*val.mean)/(n+val.n);
    n += val.n;
    return mean;
}

int my_mean::get_count(){
    return n;
}