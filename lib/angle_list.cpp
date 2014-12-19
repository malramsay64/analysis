//
//  angle_list.cpp
//  analysis
//
//  Created by Malcolm Ramsay on 19/12/2014.
//  Copyright (c) 2014 Malcolm Ramsay. All rights reserved.
//

#include "angle_list.h"

using namespace std;

angle_list::angle_list(){
    a = std::vector<my_mean>(0);
};

int angle_list::push(double angle){
    return push(angle,0);
    double dA;
    vector<my_mean>::iterator i;
    for (i = a.begin(); i != a.end(); i++){
        dA = atan2(sin(angle - (*i).get_mean()), cos(angle-(*i).get_mean()));
        if (dA < 0.05){
            (*i).add((*i).get_mean()+dA);
            cout << "Colour " << (int) (i - a.begin() + 1)  << " Angle " << angle << endl;
            return (int) (i - a.begin() + 1);
        }
    }
    a.push_back(my_mean());
    a.back().add(angle);
    dist.push_back(my_mean());
    dist.back().add(0);
    cout << "Colour " << (int) (i - a.begin() + 1)  << " Angle " << angle << endl;
    return (int) a.size();
}

int angle_list::push(double angle, double d){
    double dA;
    double dD;
    for (int i = 0; i < a.size(); i++){
        dA = atan2(sin(angle - a.at(i).get_mean()), cos(angle-a.at(i).get_mean()));
        dD = my_mod(d - dist.at(i).get_mean(), dist.at(i).get_mean());
        if (fabs(dA) < deltaA && fabs(dD) < deltaD){
            a.at(i).add(a.at(i).get_mean()+dA);
            dist.at(i).add(dist.at(i).get_mean()+dD);
            return i+1;
        }
    }
    a.push_back(my_mean());
    dist.push_back(my_mean());
    a.back().add(angle);
    dist.back().add(d);
    return (int) a.size();
    
}

int angle_list::print(std::ostream * file){
    //std::cout << a.size() << std::endl;
    for (int i = 0; i < a.size(); i++){
        *file << (a.at(i).get_mean()) + PI/2 << " " << dist.at(i).get_mean() << " " << a.at(i).get_count() << std::endl;
    }
    return 0;
}
