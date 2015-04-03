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

int pos_def_mod(int a, int n){
    int res = a - (n * int(double(a)/n));
    while (res < 0) res+=n;
    return res;
}

double max_structure_factor(std::vector<double> g, double rho, double dr){
    double max_val = 0, max = 0, res;
    for (int i = 0; i < g.size(); i++){
        res = structure_factor(2*PI/(i*dr), g, rho, dr);
        if (res > max_val){
            max_val = res;
            max = 2*PI/(i*dr);
        }
    }
    
    return max_val;
}

double structure_factor(double q, std::vector<double> g, double rho, double dr){
    double s = 0;
    for (int i = 1; i < g.size(); i++){
        s+=dr*(i*dr)*(g.at(i)-1)*j0(q*(dr*i));
    }
    s *= 2*PI*rho;
    s += 1;
    return s;
}

