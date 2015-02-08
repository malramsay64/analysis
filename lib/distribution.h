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

template <class type>
class distribution{
    std::vector<type> dist;
    type sum;
    int elements;
    double delta_r;
  public:
    
    distribution();
    distribution(int);
    distribution(int, double);
    
    void add(int);
    void add(double);
    void add(distribution);
    //void add(distribution<int>);
    void add(std::vector<int>);
    double get_mean();
    int get_elements();
    int get_sum();
    int get_size();
    double get_delta_r();
    type at(int);
    
};

template <class type>
distribution<type>::distribution(int size){
    dist = std::vector<type>(size, 0);
    delta_r = 1;
    sum = 0;
    elements = 0;
}

template <class type>
distribution<type>::distribution(int size, double range){
    dist = std::vector<type>(size, 0);
    delta_r = range/size;
    sum = 0;
    elements = 0;
}

template <class type>
void distribution<type>::add(int i){
    dist.at(i)++;
    sum += i;
    elements++;
}

template <class type>
void distribution<type>::add(double d){
    if (d < dist.size()*delta_r){
        dist.at(int(d/delta_r))++;
        sum += int(d/delta_r);
        elements++;
    }
}

template <class type>
void distribution<type>::add(distribution<type> d){
    // Largest size
    if (d.get_size() > get_size()){
        dist.resize(d.get_size(),0);
    }
    // Same delta_r
    if (d.get_delta_r() != get_delta_r()){
        return;
    }
    for (int i = 0; i < d.get_size(); i++){
        dist.at(i) += d.at(i);
    }
    sum += d.get_sum();
    elements += d.get_elements();
}

/*
template <>
void distribution<double>::add(distribution<int> d) {
    for (auto &e: v){
        dist.at(e) += 1.0/v.size();
    }
}
*/

template <class type>
void distribution<type>::add(std::vector<int> v) {
    for (auto &e: v){
        dist.at(e) += 1.0/v.size();
    }
    elements++;
}

template <class type>
double distribution<type>::get_mean(){
    return double(sum)/elements;
}

template <class type>
int distribution<type>::get_elements(){
    return elements;
}

template <class type>
int distribution<type>::get_sum(){
    return sum;
}

template <class type>
int distribution<type>::get_size(){
    return (int) dist.size();
}

template <class type>
double distribution<type>::get_delta_r(){
    return delta_r;
}

template <class type>
type distribution<type>::at(int i){
    return dist.at(i);
}

template <class type>
int print_distribution(distribution<type> *d, std::string filename){
    std::ofstream file;
    file.open(filename.c_str());
    for (int i = 0; i < d->get_size(); i++){
        file << i << " " << d->at(i) << std::
        endl;
    }
    
    return 0;
}

int print_radial_distribution(distribution<int> *d, std::string filename, int nmol, double frame_area);

#endif /* defined(MY_DISTRIBUTION) */
