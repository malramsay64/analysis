//
//  molecule.cpp
//  analysis
//
//  Created by Malcolm Ramsay on 7/12/2014.
//  Copyright (c) 2014 Malcolm Ramsay. All rights reserved.
//

#include "molecule.h"

using namespace std;

molecule::molecule (){
    atoms.reserve(MOL_SIZE);
    my_neighbours.reserve(MAX_MOL_CONTACTS);
    nint = std::vector<int>(MAX_MOL_CONTACTS,0);
    contacts = 0;
    colour = 0;
    type = 0;
    com = vect();
}

void molecule::add_neighbour(molecule *m){
    int pos;
    pos = (int) (find(my_neighbours.begin(), my_neighbours.end(), m) - my_neighbours.begin());
    // Molecule already added
    if (pos != my_neighbours.size()){
        nint[pos]++;
    }
    else {
        my_neighbours.push_back(m);
        nint[pos]++;
    }
    contacts++;
}

void molecule::delete_neighbours(){
    my_neighbours = vector<molecule*>(0,0);
    nint = vector<int>(MAX_MOL_CONTACTS,0);
    contacts = 0;
}

double molecule::mass(){
    double mass = 0;
    std::vector<particle *>::iterator i;
    for ( i = atoms.begin(); i != atoms.end(); ++i){
        mass += (*i)->mass;
    }
    return mass;
}

vect molecule::COM(){
    if (com != vect()){
        return com;
    }
    else{
        com = calc_COM();
        return com;
    }
}

vect molecule::calc_COM(){
    double total = 0;
    vect theta, xi, zeta, out;
    for ( int i = 0; i < nump(); ++i){
        theta = atoms[i]->pos_vect();
        xi += atoms[i]->mass*cos(theta);
        zeta += atoms[i]->mass*sin(theta);
        total += atoms[i]->mass;
    }
    xi /= total;
    zeta /= total;
    theta = atan2(-zeta,-xi) + PI;
    out = theta;
    return out;
}

vect molecule::atom_pos(int i){
    return atoms[i]->pos_vect();
}

int molecule::nump(){
    return (int) atoms.size();
}

int molecule::uniqc(){
    return (int) my_neighbours.size();
}

int molecule::num_contacts(){
    return contacts;
}

int molecule::num_neighbours(){
    return uniqc();
}

distribution<int> molecule::pairing(){
    distribution<int> d = distribution<int>(MAX_MOL_CONTACTS);
    for (auto &n: nint){
        if ( n != 0){
            d.add(n);
        }
    }
    return d;
}

int molecule::max_pairing(){
    int max = 0;
    for (auto n: nint){
        if (n > max){
            max = n;
        }
    }
    return max;
}

int molecule::get_colour(){
    return colour;
}

int molecule::set_colour(int i){
    colour = i;
    return 0;
}

int molecule::same_period(){
    vect com = COM();
    vect d;
    vector<particle *>::iterator i;
    for (i = atoms.begin(); i != atoms.end(); i++){
        d = direction(COM(), (*i)->pos_vect());
        (*i)->set_pos(com + d);
    }
    return 0;
}

int molecule::index(){
    return id-1;
}
