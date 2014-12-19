//
//  angle_list.h
//  analysis
//
//  Created by Malcolm Ramsay on 7/12/2014.
//  Copyright (c) 2014 Malcolm Ramsay. All rights reserved.
//

#ifndef ANGLE_LIST_H
#define ANGLE_LIST_H

#include <vector>
#include <iostream>
#include <cmath>


#include "my_mean.h"
#include "functions.h"
#include "constants.h"

static double deltaA = 5*PI/180;
static double deltaD = 0.15;

class angle_list{
    std::vector<my_mean> a;
    std::vector<my_mean> dist;
public:
    
    angle_list();
    
    int push(double angle);
    int push(double angle, double d);
    int print(std::ostream * file);
};


#endif /* defined(ANGLE_LIST_H) */
