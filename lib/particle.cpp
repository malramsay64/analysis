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

molecule::molecule (){
    atoms.reserve(MOL_SIZE);
    my_neighbours.reserve(MAX_MOL_CONTACTS);
    nint = std::vector<int>(MAX_MOL_CONTACTS,0);
    contacts = 0;
    colour = 0;
    traversed = false;
    type = 0;
}
/*
molecule::~molecule(){
    delete &atoms;
    delete &my_neighbours;
    delete &nint;
}
*/
void molecule::add_neighbour(molecule *m){
    unsigned int pos;
    pos = find(my_neighbours.begin(), my_neighbours.end(), m) - my_neighbours.begin();
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

double molecule::mass(){
    double mass = 0;
    std::vector<particle *>::iterator i;
    for ( i = atoms.begin(); i != atoms.end(); ++i){
        mass += (*i)->mass;
    }
    return mass;
}

vect molecule::COM(){
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
    return atoms.size();
}

int molecule::uniqc(){
    return my_neighbours.size();
}

int molecule::num_contacts(){
    return contacts;
}

int molecule::num_neighbours(){
    return uniqc();
}

int molecule::pairing(molecule *b){
    int pos;
    pos = std::find(my_neighbours.begin(), my_neighbours.end(), b) - my_neighbours.begin();
    return nint[pos];
}

int molecule::get_colour(){
    return colour;
}

int molecule::set_colour(int i){
    colour = i;
    return 0;
}

int molecule::graph_colour(){
    colour = 0;
    
    //int sum = 0;
    vector<int> val (6,0);
    vector<molecule*>::iterator i;
    //cout << "Sizeof val " << val.size() << endl; 
    for (i = my_neighbours.begin(); i != my_neighbours.end(); i++){
        //cout << "Colour: " << (*i)->get_colour() << endl;
        if ((*i)->get_colour()){
            val.at((*i)->get_colour()-1) = 1;
        }
        //sum += (int) pow(2,(*i)->get_colour()-1);
    }
    int j;
    for (j = 0; j < val.size(); j++){    
        if (val.at(j) == 0){
            j++;
            break;
        }
    }
    colour = j;
    //colour = (int) floor(log2(15-sum)+1);
    //cout << colour << endl;
    //colour = 0;
    return colour;
}

bool molecule::get_traversed(){
    return traversed;
}

void molecule::traverse(){
    traversed = true;
}

void molecule::reset_traverse(){
    traversed = false;
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

