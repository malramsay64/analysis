//
//  main.cpp
//  radial
//
//  Created by Malcolm Ramsay on 5/02/2015.
//  Copyright (c) 2015 Malcolm Ramsay. All rights reserved.
//

#include <iostream>
#include <vector>
#include "vect.h"

using namespace std;

double rand_range(double max){
    return (rand()/double(RAND_MAX))*max;
}

int main(int argc, const char * argv[]) {

    int num_p = 1000;
    int points = 100;
    double x_range = 2*PI;
    double r_range = 1;
    double dr = r_range/points;
    
    vector<vect> particles = vector<vect>(num_p, vect());
    vector<int> r_dist = vector<int>(points, 0);
    
    // Generate random distribution of particles
    for (auto &p: particles){
        p = vect(rand_range(x_range), rand_range(x_range));
    }
    
    // Calculate radial distributions
    double d;
    for (auto p: particles){
        for (auto q: particles){
            d = dist(p,q);
            if (d < r_range){
                r_dist.at(int(d/dr))++;
            }
        }
    }
    
    // Output
    double area;
    double density = num_p/(x_range*x_range);
    for (int i = 0; i < points; i++){
        area = 2*PI*(i*dr)*dr;
        cout << i*dr << " " << r_dist.at(i) << " " << r_dist.at(i)/area << " " << r_dist.at(i)/(area*num_p)\
        << " " <<  r_dist.at(i)/(area*num_p*density) << endl;
    }
    
    
    return 0;
}
