//
//  particle.cpp
//  analysis
//
//  Created by Malcolm Ramsay on 7/12/2014.
//  Copyright (c) 2014 Malcolm Ramsay. All rights reserved.
//

#include "Particle.h"

using namespace LAlgebra;

Particle::Particle(){
    mass = 1;
    molid = 0;
    type = 0;
    id = 0;
    radius = 0;
    my_neighbours = std::vector<Particle *>{};
    pos = Vector<2>{};
}

Particle::Particle(const Particle &p){
    mass = p.mass;
    pos = p.pos;
    molid = p.molid;
    type = p.type;
    id = p.id;
    radius = p.radius;
    my_neighbours = p.my_neighbours;
}

Particle::Particle(const particle_vars &p) {
    mass = p.mass;
    pos.r = {p.xpos, p.ypos};
    molid = p.molid;
    type = p.type;
    id = p.id;
    radius = p.radius;
}

void Particle::delete_neighbours(){
    my_neighbours = std::vector<Particle *>(0, nullptr);
}

int Particle::index() const{
    return id-1;
}

int Particle::mol_index() const{
    return molid-1;
}

size_t Particle::numn() const {
    return my_neighbours.size();
}

void Particle::append(Particle *p){
    my_neighbours.push_back(p);
}

void Particle::set_pos(const Vector<2> &v){
    pos = v;
}

Vector<2> Particle::pos_vect() const{
    return pos;
}

int Particle::order(){
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
    // Not including bonded particle
    if (type == 1 && large == 3 && small == 3 ){
        return 1;
    }
    else if (type == 2 && large == 3 && small == 1 ){
        return 1;
    }
    return 0;
}

bool Particle::operator> (const Particle &b) const {
    return id > b.id;
}

bool Particle::operator>= (const Particle &b)const {
    return id >= b.id;
}

bool Particle::operator< (const Particle &b) const{
    return id < b.id;
}

bool Particle::operator<= (const Particle &b) const{
    return id <= b.id;
}

bool Particle::operator== (const Particle &b) const{
    return id == b.id;
}

bool Particle::operator!= (const Particle &b) const {
    return id != b.id;
}
