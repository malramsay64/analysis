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

std::string split(std::string s, char delim){
    int p;
    p = s.find_last_of(delim);
    return s.substr(0,p);
}

std::string basename(std::string s){
    size_t p;
    char delim = '/';
    p = s.find_last_of(delim); 
    return s.substr(p+1);
}

double legendre(int l, double x){
    switch (l){
        case 1:
            return x;
        case 2:
            return 0.5*(3*pow(x,2) - 1);
        case 3:
            return 0.5*(5*pow(x,3) - 3*x);
        case 4:
            return (1./8)*(35*pow(x,4) - 30*pow(x,2) + 3);
        default:
            return 0;
    }
}
