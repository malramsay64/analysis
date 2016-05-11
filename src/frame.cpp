#include "frame.h"


using namespace std;

Frame::Frame(){
    num_molecules = 0;
    atoms = 0;
    molecules = vector<Molecule>();
    particles = vector<Particle>();
    coloured = false;
    neighbours = false;
}

Frame::Frame(Frame const &frame){
    num_molecules = frame.num_molecules;
    atoms = frame.atoms;
    molecules = frame.molecules;
    particles = frame.particles;
    coloured = frame.coloured;
    neighbours = frame.neighbours;
    a = frame.a;
    b = frame.b;
    theta = frame.theta;
    timestep = frame.timestep;
    update_links();
}

double Frame::dist(Vector2d v1, Vector2d v2){
    return direction(v1,v2).length();
}

Vector2d Frame::direction(Vector2d v1, Vector2d v2){
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
}

int Frame::num_mol(){
    return (int) molecules.size();
}

void Frame::add_particle(Particle *p){
    particles.push_back(Particle(*p));
}

void Frame::add_link(int molpos, int ppos){
    add_link(molpos, ppos, true);
}

void Frame::add_link(int molpos, int ppos, bool b_sort=true){
    Molecule *m = &molecules[molpos];
    m->atoms.push_back(&particles[ppos]);
    if (b_sort){
        sort(m->atoms.begin(), m->atoms.end());
    }
}

void Frame::update_links(){
    for (auto &m: molecules){
        // Links to particles
        for (int i = 0; i < m.num_particles(); i++){
            m.atoms.at(i) = &particles.at(m.atoms.at(i)->index());
        }
        // Links to neighbours
    }
}

Molecule Frame::at(int i){
    return molecules.at(i);
}

void Frame::set_timestep(int t){
    timestep = t;
}

int Frame::get_timestep(){
    return timestep;
}

double Frame::get_time(){
    return get_timestep()*STEP_SIZE;
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

Vector2d Frame::cartesian(Vector2d v){
    v = v/(2*PI);
    v.x = v.x*a + v.y*b*cos(theta);
    v.y = v.y*b*sin(theta);
    return v;
}

Vector2d Frame::fractional(Vector2d v){
    v.x = v.x*(1/a) + v.y*(-cos(theta)/(a*sin(theta)));
    v.y = v.y/(b*sin(theta));
    v = v*2*PI;
    return v;
}

