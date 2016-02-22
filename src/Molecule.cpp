//
//  molecule.cpp
//  analysis
//
//  Created by Malcolm Ramsay on 7/12/2014.
//  Copyright (c) 2014 Malcolm Ramsay. All rights reserved.
//

#include "Molecule.h"

using namespace LAlgebra;

void Molecule::add_neighbour(Molecule *m){
    my_neighbours[m]++;
    contacts++;
}

void Molecule::delete_neighbours(){
    my_neighbours = std::map<Molecule *,int>{};
    contacts = 0;
    for (auto &p: atoms){
        p->delete_neighbours();
    }
}

void Molecule::delete_mol_neighbours(){
    my_neighbours = std::map<Molecule *, int>{};
    contacts = 0;   
}

double Molecule::get_mass() const {
    double mass = 0;
    for (auto const &a: atoms){
        mass += a->mass;
    }
    return mass;
}

Vector Molecule::get_COM() const {
    return pos_def_mod(com, 2*PI);
}

Vector Molecule::get_COM(){
    if (com != Vector{}){
        return pos_def_mod(com,2*PI);
    }
    else{
        calc_COM();
        return com;
    }
}

Vector Molecule::moved_COM() const{
    return com;
}

Vector Molecule::calc_COM(){
    double total = 0;
    Vector theta, xi, zeta;
    for (auto p: atoms){
        theta = p->pos_vect();
        xi += p->mass*cos(theta);
        zeta += p->mass*sin(theta);
        total += p->mass;
    }
    xi /= total;
    zeta /= total;
    com = atan2(-zeta,-xi) + PI;
    return com;
}

Vector Molecule::update_COM(){
    Vector old{com};
    com = Vector{old + direction(old, calc_COM())};
    return com;
}

Particle * Molecule::get_large() const{
    for (auto p: atoms){
        if (p->type == 1){
            return p;
        }
    }
    return nullptr;
}

double Molecule::update_orientation(double angle){
    double delta;
    delta = (angle-PI) - orientation;
    delta = atan2(sin(delta),cos(delta));
    orientation += delta;
    rotation += delta;
    return delta;
}

int Molecule::set_orientation(double angle){
    orientation = angle - PI;
    rotation = 0;
    return 0;
}

double Molecule::get_orientation() const{
    return orientation;
}

Vector Molecule::get_orient_vect() const{
    return Vector(std::valarray<double>{std::sin(get_orientation()), std::cos(get_orientation())});
}

double Molecule::get_rotation() const{
    return rotation;
}

Vector Molecule::atom_pos(int i) const{
    return atoms.at(i)->pos_vect();
}

size_t Molecule::num_particles() const{
    return atoms.size();
}

size_t Molecule::uniqe_contacts() const{
    return my_neighbours.size();
}

int Molecule::num_contacts() const{
    return contacts;
}

int Molecule::num_neighbours() const{
    return uniqe_contacts();
}

std::map <int, int> Molecule::pairing() const {
    std::map<int, int> d{};
    for (auto &n: my_neighbours){
        d[n.second]++;
    }
    return d;
}

int Molecule::max_pairing() const{
    int max = 0;
    for (auto &n: my_neighbours){
        if (n.second > max){
            max = n.second;
        }
    }
    return max;
}

int Molecule::same_period(){
    Vector com = get_COM();
    Vector d;
    std::vector<Particle *>::iterator i;
    for (auto &i: atoms){
        d = direction(get_COM(), i->pos_vect());
        i->set_pos(com + d);
    }
    return 0;
}

int Molecule::index() const{
    return id-1;
}

std::vector<Molecule *> Molecule::get_neighbours() {
    std::vector<Molecule *> ret{};
    for (auto i: my_neighbours){
        ret.push_back(i.first);
    }
    return ret;
}