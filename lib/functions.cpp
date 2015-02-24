//
//  functions.cpp
//  analysis
//
//  Created by Malcolm Ramsay on 7/12/2014.
//  Copyright (c) 2014 Malcolm Ramsay. All rights reserved.
//

#include <assert.h>
#include "functions.h"



double dot_product(double *v1, double *v2, int len){
    double sum = 0;
    for (int i = 0; i < len; ++i){
        sum += v1[i]*v2[i];
    }
    return sum;
}

double dist(double x1, double x2){
    double x = x1 - x2;
    return acos(cos(x));
}

double legendre(int l, double x){
    switch (l){
        case 1:
            return x;
            break;
        case 2:
            return 2*pow(x,2) - 1;
            break;
        default:
            return 0;
            break;
        //case 3:
        //    return 0.5*(5*pow(x,3) - 3*x);
        //case 4:
        //    return (1./8)*(35*pow(x,4) - 30*pow(x,2) + 3);
        //default:
        //    return 0;
    }
}

double my_mod(double a, double b){
    double ret = fmod(a, b);
    if (ret < 0){
        ret += b;
    }
    return ret;
}


