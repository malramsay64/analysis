//
//  particle.cpp
//  analysis
//
//  Created by Malcolm Ramsay on 7/12/2014.
//  Copyright (c) 2014 Malcolm Ramsay. All rights reserved.
//

/*
 * particle.cpp
 * Copyright (C) 2014 malcolm <malcolm@macbook.local>
 *
 * Distributed under terms of the MIT license.
 */

#include "particle.h"

using namespace std;

particle::particle (){
    my_neighbours.reserve(MAX_ATOM_CONTACTS);
    mass = 1;
    pos = vect(0,0);
}

int particle::index(){
    return id-1;
}

int particle::m_i(){
    return molid-1;
}

int particle::numn(){
    return (int) my_neighbours.size();
}

void particle::append(particle *p){
    my_neighbours.push_back(p);
}

int particle::set_pos(vect v){
    pos = vect(v);
    return 0;
}

vect particle::pos_vect(){
    return pos;
}

int particle::order(){
    int small = 0, large = 0;
    for (auto &p: my_neighbours){
        if (p->type == 1){
            large++;
        }
        else{
            small++;
        }
    }
    // Ordering for d=1.637556 dimers
    if (type == 1 && large == 3 && small == 4 ){
        return 1;
    }
    else if (type == 2 && large == 4 && small == 1 ){
        return 1;
    }
    return 0;
}

bool particle::operator> (const particle &b){
    return id > b.id;
}

bool particle::operator>= (const particle &b){
    return id >= b.id;
}

bool particle::operator< (const particle &b){
    return id < b.id;
}

bool particle::operator<= (const particle &b){
    return id <= b.id;
}

bool particle::operator== (const particle &b){
    return id == b.id;
}

bool particle::operator!= (const particle &b){
    return id != b.id;
}

bool operator < (const particle &a, const particle &b){
    return a.id < b.id;
}

