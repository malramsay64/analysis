//
//  distribution.cpp
//  analysis
//
//  Created by Malcolm Ramsay on 6/02/2015.
//  Copyright (c) 2015 Malcolm Ramsay. All rights reserved.
//

#include "distribution.h"

using namespace std;

distribution::distribution(int size){
    dist = vector<int>(size, 0);
    delta_r = 1;
    sum = 0;
    elements = 0;
}

distribution::distribution(int size, double range){
    dist = vector<int>(size, 0);
    delta_r = range/size;
    sum = 0;
    elements = 0;
}

void distribution::add(int i){
    dist.at(i)++;
    sum += i;
    elements++;
}

void distribution::add(double d){
    if (d < dist.size()*delta_r){
        dist.at(int(d/delta_r))++;
        sum += int(d/delta_r);
        elements++;
    }
}

int distribution::add(distribution d){
    // Largest size
    if (d.get_size() > get_size()){
        dist.resize(d.get_size(),0);
    }
    // Same delta_r
    if (d.get_delta_r() != get_delta_r()){
        return 1;
    }
    for (int i = 0; i < d.get_size(); i++){
        dist.at(i) += d.at(i);
    }
    sum += d.get_sum();
    elements += d.get_elements();
    return 0;
}

double distribution::get_mean(){
    return double(sum)/elements;
}

int distribution::get_elements(){
    return elements;
}

int distribution::get_sum(){
    return sum;
}

int distribution::get_size(){
    return (int) dist.size();
}

double distribution::get_delta_r(){
    return delta_r;
}

int distribution::at(int i){
    return dist.at(i);
}

int print_distribution(distribution *d, string filename){
    ofstream file;
    file.open(filename.c_str());
    for (int i = 0; i < d->get_size(); i++){
        file << i << " " << d->at(i) << endl;
    }
    
    return 0;
}

int print_radial_distribution(distribution *d, string filename, int nmol, double frame_area){
    ofstream file;
    file.open(filename.c_str());
    double area;
    double density = nmol/frame_area;
    for (int i = 0; i < d->get_size(); i++){
        area = PI*(i*d->get_delta_r())*d->get_delta_r();
        file << i*d->get_delta_r() << " " << d->at(i)/(area*nmol*density) << endl;
    }
    return 0;
}