//
//  molecule.cpp
//  analysis
//
//  Created by Malcolm Ramsay on 7/12/2014.
//  Copyright (c) 2014 Malcolm Ramsay. All rights reserved.
//

#include "Molecule.h"

using namespace LAlgebra;

Molecule::Molecule(){
}

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


Vector2d Molecule::get_COM() const {
    return pos_def_mod(com, 2*PI);
}

Vector2d Molecule::get_COM(){
    if (com != Vector2d{}){
        return pos_def_mod(com,2*PI);
    }
    else{
        calc_COM();
        return com;
    }
}

Vector2d Molecule::moved_COM() const{
    return com;
}

Vector2d Molecule::calc_COM(){
    double total = 0;
    Vector2d theta, xi, zeta;
    for (auto p: atoms){
        theta = p->pos_vect();
        xi += p->mass*cos(theta);
        zeta += p->mass*sin(theta);
        total += p->mass;
    }
    xi /= total;
    zeta /= total;
    theta = atan2(-zeta,-xi) + PI;
    com = theta;
    return com;
}

Vector2d Molecule::update_COM(){
    Vector2d old = Vector2d(com);
    com = Vector2d(old + direction(old, calc_COM()));
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

Vector2d Molecule::get_orient_vect() const{
    return Vector2d(sin(get_orientation()), cos(get_orientation()));
}

double Molecule::get_rotation() const{
    return rotation;
}

Vector2d Molecule::atom_pos(int i) const{
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
    Vector2d com = get_COM();
    Vector2d d;
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

bool Molecule::operator!=(const Molecule &b) const {
    return id != b.id;
}

bool Molecule::operator==(const Molecule &b) const {
    return id == b.id;
}

bool Molecule::operator<(const Molecule &b) const {
    return id < b.id;
}

bool Molecule::operator<=(const Molecule &b) const {
    return id <= b.id;
}

bool Molecule::operator>(const Molecule &b) const {
    return id > b.id;
}

bool Molecule::operator>=(const Molecule &b) const {
    return id >= b.id;
}