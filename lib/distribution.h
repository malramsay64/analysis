//
//  distribution.h
//  analysis
//
//  Created by Malcolm Ramsay on 6/02/2015.
//  Copyright (c) 2015 Malcolm Ramsay. All rights reserved.
//

#ifndef MY_DISTRIBUTION
#define MY_DISTRIBUTION

#include <vector>
#include <fstream>
#include "constants.h"

class distribution{
    std::vector<int> dist;
    int sum;
    int elements;
    double delta_r;
  public:
    
    distribution();
    distribution(int);
    distribution(int, double);
    
    void add(int);
    void add(double);
    int add(distribution);
    double get_mean();
    int get_elements();
    int get_sum();
    int get_size();
    double get_delta_r();
    int at(int);
};

int print_distribution(distribution *d, std::string filename);
int print_radial_distribution(distribution *d, std::string filename, int nmol, double frame_area);
#endif /* defined(MY_DISTRIBUTION) */
