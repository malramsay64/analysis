#include "frame.h"
#include <sstream>

using namespace std;

Frame::Frame(){
    num_molecules = 0;
    molecules = vector<molecule>();
    particles = vector<particle>();
    coloured = false;
    neighbours = false;
}

double Frame::dist(vect v1, vect v2){
    return direction(v1,v2).length();
}

vect Frame::direction(vect v1, vect v2){
    return cartesian(::direction(v1,v2));
}

void Frame::set_atoms(int a){
    atoms = a;
    particles.reserve(atoms);
}

int Frame::num_atoms(){
    return atoms;
}

void Frame::set_num_mol(int n){
    num_molecules = n;
    molecules.resize(n);
    m = vector<mutex>(n);
}

int Frame::num_mol(){
    return (int) molecules.size();
}

void Frame::add_particle(particle *p){
    particles.push_back(particle(*p));
}

void Frame::add_link(int molpos, int ppos){
    add_link(molpos, ppos, true);
}

void Frame::add_link(int molpos, int ppos, bool b_sort=true){
    molecule *m = &molecules[molpos];
    m->atoms.push_back(&particles[ppos]);
    if (b_sort){
        sort(m->atoms.begin(), m->atoms.end());
    }
}

void Frame::set_timestep(int t){
    timestep = t;
}

int Frame::get_timestep(){
    return timestep;
}


// Crystal Coordinates
int Frame::set_crys(double a, double b, double theta){
    this->a = a;
    this->b = b;
    this->theta = theta;
    return 0;
}

double Frame::get_a(){
    return a;
}

double Frame::get_b(){
    return b;
}

double Frame::get_theta(){
    return theta;
}

double Frame::get_height(){
    return b*sin(theta);
}

double Frame::get_tilt(){
    return b*cos(theta);
}

double Frame::get_area(){
    return get_a()*get_height();
}

double Frame::get_density(){
    return num_atoms()/get_area();
}

vect Frame::cartesian(vect v){
    v = v/(2*PI);
    v.x = v.x*a + v.y*b*cos(theta);
    v.y = v.y*b*sin(theta);
    return v;
}

vect Frame::fractional(vect v){
    v.x = v.x*(1/a) + v.y*(-cos(theta)/(a*sin(theta)));
    v.y = v.y/(b*sin(theta));
    v = v*2*PI;
    return v;
    
}

