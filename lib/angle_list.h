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
    
    int push(double d){
        for (int i = 0; i < a.size(); i++){
            if (abs(d-a.at(i).get_mean()) < deltaA){
                a.at(i).add(d-a.at(i).get_mean());
                return i+1;
            }
        }
        a.push_back(my_mean());
        a.back().add(d);
        return a.size();
    }
    
    int push(double angle, double d){
        //std::cout << angle << " " << d << std::endl;
        angle = my_mod(angle - PI/2,PI);
        //std::cout << "Angle " << angle << std::endl;
        double delta;
        for (int i = 0; i < a.size(); i++){
            delta = acos(cos(2*(angle-a.at(i).get_mean())))/2;
            if (delta < deltaA && abs(d-dist.at(i).get_mean()) < deltaD){
                a.at(i).add(a.at(i).get_mean()+delta);
                return i+1;
            }
        }
        if (d > 0){
            a.push_back(my_mean());
            dist.push_back(my_mean());
            a.back().add(angle);
            dist.back().add(d);
        }
        return a.size();
 
    }
    
    int print(std::ostream * file){
        if (dist.size() == a.size()){
            for (int i = 0; i < a.size(); i++){
                *file << (a.at(i).get_mean())*180/PI << " " << dist.at(i).get_mean() << " " << a.at(i).get_count() << std::endl;
            }
        }
        else {
            for (int i = 0; i < a.size(); i++){
                *file << (a.at(i).get_mean())*180/PI << std::endl;
            }
        }
        return 0;
    }
};


#endif /* defined(ANGLE_LIST_H) */
