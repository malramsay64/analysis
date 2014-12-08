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
    traversed = false;
}

int particle::index(){
    return id-1;
}

int particle::m_i(){
    return molid-1;
}

int particle::numn(){
    return my_neighbours.size();
}

int particle::n_large(){
    vector<particle *>::iterator p;
    int count = 0;
    for (p = my_neighbours.begin(); p!= my_neighbours.end();p++){
        if ((*p)->type == 1){
            count++;
        }
    }
    return count;
}

int particle::n_small(){
    vector<particle *>::iterator p;
    int count = 0;
    for (p = my_neighbours.begin(); p!= my_neighbours.end();p++){
        if ((*p)->type == 2){
            count++;
        }
    }
    return count;
}

void particle::append(particle *p){
    my_neighbours.push_back(p);
}

double particle::xpos(){
    return pos[XX];
}

int particle::set_xpos(double x){
    pos[XX] = x;
    return 0;
}

double particle::ypos(){
    return pos[YY];
}

int particle::set_ypos(double y){
    pos[YY] = y;
    return 0;
}

int particle::set_pos(vect v){
    set_xpos(v.x);
    set_ypos(v.y);
    return 0;
}

vect particle::pos_vect(){
    return vect(pos[XX], pos[YY]);
}

bool particle::get_traversed(){
    return traversed;
}

void particle::traverse(){
    traversed = true;
}

void particle::reset_traverse(){
    traversed = false;
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

