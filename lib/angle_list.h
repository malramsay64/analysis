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
#include "constants.h"

static double deltaA = 5*PI/180;
static double deltaD = 0.15;

class angle_list{
    std::vector<my_mean> a;
    std::vector<my_mean> dist;
public:
    
    angle_list(){
       a = std::vector<my_mean>(0);
    };
    
    int push(double angle){
        return push(angle, 0);
    }
    
    int push(double angle, double d){
        //std::cout << angle << " " << d << std::endl;
        angle = angle;
        //std::cout << "Angle " << angle << std::endl;
        double dA;
        double dD;
        for (int i = 0; i < a.size(); i++){
            dA = asin(sin(angle-a.at(i).get_mean()));
            dD = my_mod(d - dist.at(i).get_mean(), dist.at(i).get_mean());
            if (fabs(dA) < deltaA && fabs(dD) < deltaD){
                a.at(i).add(a.at(i).get_mean()+dA);
                dist.at(i).add(dist.at(i).get_mean()+dD);
                                return i+1;
            }
        }
        if (d > 0){
            if (fabs(angle) < deltaA){
                std::cout << "Zero Angle"<< std::endl;
            }
            a.push_back(my_mean());
            dist.push_back(my_mean());
            a.back().add(angle);
            dist.back().add(d);
        }
        return a.size();
 
    }
    
    int print(std::ostream * file){
        std::cout << a.size() << std::endl;
        for (int i = 0; i < a.size(); i++){
            *file << (a.at(i).get_mean()) + PI/2 << " " << dist.at(i).get_mean() << " " << a.at(i).get_count() << std::endl;
        }
        return 0;
    }
};


#endif /* defined(ANGLE_LIST_H) */
