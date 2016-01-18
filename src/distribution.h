//
//  distribution.h
//  analysis
//
//  Created by Malcolm Ramsay on 6/02/2015.
//  Copyright (c) 2015 Malcolm Ramsay. All rights reserved.
//

#include <vector>
#include <fstream>
#include "constants.h"
#include "my_mean.h"
#include "output.h"

#ifndef MY_DISTRIBUTION
#define MY_DISTRIBUTION


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
    void add(int, int);
    void add(distribution);
    //void add(distribution<int>);
    void add(std::vector<int>);
    double get_mean();
    int get_elements();
    int get_sum();
    int get_size();
    double get_delta_r();
    type at(int);
    double fraction_at(int);
    std::vector<type> get_dist(double scale);
    
};

template <class type>
distribution<type>::distribution(){
    dist = std::vector<type>(0);
    elements = 0;
}

template <class type>
distribution<type>::distribution(int size){
    dist = std::vector<type>(size, 0);
    delta_r = 1;
    sum = 0;
    elements = 0;
}

template <> inline
distribution<my_mean>::distribution(int size){
    dist = std::vector<my_mean>(size);
    delta_r = 1;
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
        dist.at(pos_def_mod(int(d/delta_r),get_size()))++;
        sum += int(d/delta_r);
        elements++;
    }
}

template <class type>
void distribution<type>::add(int i, int j){
    dist.at(i) += j;
    sum += i*j;
    elements += j;
}

template <> inline
void distribution<my_mean>::add(int i, int j){
    dist.at(i).add(j);
    //sum = sum + j;
    elements++;
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
double distribution<type>::fraction_at(int i){
    return dist.at(i)/double(elements);
}

template <class type>
std::vector<type> distribution<type>::get_dist(double scale){
    std::vector<type> t = std::vector<type>(dist);
    for (auto &i: t){
        i *= scale;
    }
    return t;
}

template <class type>
int print_distribution(distribution<type> *d, std::string filename){
    std::ofstream file;
    file.open(filename.c_str());
    for (int i = 0; i < d->get_size(); i++){
        file << i << " " << d->at(i) << std::endl;
    }
    
    return 0;
}

template <class type> inline
int print_frac_distribution(distribution<type> *d, std::string filename){
    std::ofstream file;
    file.open(filename.c_str());
    for (int i = 0; i < d->get_size(); i++){
        file << i << " " << d->fraction_at(i) << std::endl;
    }
    
    return 0;
}


template <> inline
int print_distribution(distribution<my_mean> *d, std::string filename){
    std::ofstream file;
    file.open(filename.c_str());
    for (int i = 0; i < d->get_size(); i++){
        file << i << " " << d->at(i).get_mean() << std::endl;
    }
    
    return 0;
}

template <class type>
int print_time_distribution(distribution<type> *d, int time, std::ofstream *file){
    *file << time;
    for (int i = 0; i < d->get_size(); i++){
        *file << "," << d->at(i);
    }
    *file << std::endl;
    return 0;
}


#endif /* defined(MY_DISTRIBUTION) */
