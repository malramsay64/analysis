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

#include "constants.h"

static double deltaA = 5*PI/180;
static double deltaD = 0.15;

class angle_list{
    std::vector<double> a;
    std::vector<int> count;
    std::vector<double> dist;
public:
    
    angle_list(){
        a = std::vector<double>(0,0);
    };
    
    int push(double d){
        for (int i = 0; i < a.size(); i++){
            if (abs(d-a.at(i)) < deltaA){
                a.at(i) += (d-a.at(i))/count.at(i);
                count.at(i)++;
                return i+1;
            }
        }
        a.push_back(d);
        count.push_back(1);
        return a.size();
    }
    
    int push(double angle, double d){
        for (int i = 0; i < a.size(); i++){
            if (abs(angle-a.at(i)) < deltaA && abs(d-dist.at(i)) < deltaD){
                a.at(i) += (angle-a.at(i))/count.at(i);
                dist.at(i) += (d-dist.at(i))/count.at(i);
                count.at(i)++;
                return i+1;
            }
        }
        a.push_back(angle);
        dist.push_back(d);
        count.push_back(1);
        return a.size();
 
    }
    
    int print(){
        if (dist.size() == a.size()){
            for (int i = 0; i < a.size(); i++){
                std::cout << (a.at(i)+PI)*180/PI << " " << dist.at(i) << std::endl;
            }
        }
        else {
            for (int i = 0; i < a.size(); i++){
                std::cout << (a.at(i)+PI)*180/PI << std::endl;
            }
        }
        return 0;
    }
};


#endif /* defined(ANGLE_LIST_H) */
