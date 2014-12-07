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

/*
Frame::~Frame(){
    delete &molecules;
    delete &particles;
    delete &m;
}
*/

double Frame::dist(vect v1, vect v2){
    return direction(v1,v2).length();
}

vect Frame::direction(vect v1, vect v2){
    return ::direction(v1,v2) * size();
}

vect Frame::size(){
    return vect(xlength(), ylength());
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
    return molecules.size();
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

// Querying structure

double Frame::num_contacts(vector<int> *dist){
    int n, total = 0;
    vector<molecule>::iterator m;
    for (m=molecules.begin(); m!= molecules.end(); ++m){
        n = (*m).num_contacts();
        if (n > (*dist).size()){
            (*dist).resize(n, 0);
        }
        ++(*dist)[n];
        total += n;
    }
    return total/num_molecules;
}

double Frame::num_neighbours(vector<int> *dist){
    int n, total = 0;
    for ( int i = 0; i < num_mol(); ++i){
        n = molecules[i].num_neighbours();
        if (n > (*dist).size()){
            (*dist).resize(n, 0);
        }
        ++(*dist)[n];
        total += n;
    }
    return total/num_molecules;
}

double Frame::pairing(vector<int> *dist){
    std::vector<molecule>::iterator mols;
    int n = 0; 
    int total = 0;
    for ( mols = molecules.begin(); mols != molecules.end(); ++mols){
        for (unsigned int i = 0; i < (*mols).my_neighbours.size(); i++){
            (*dist)[(*mols).nint[i]]++;
            total += (*mols).nint[i];
            n++;
        }
    }
    return (double) total/n;
}

int Frame::reset_traverse(){
    vector<particle>::iterator i;
    for (i = particles.begin(); i != particles.end(); i++){
        (*i).reset_traverse();
    }
    return 0;
}

int Frame::set_coloured(){
    coloured = true;
    return 0;
}

int Frame::set_neighbours(){
    neighbours = true;
    return 0;
}

int Frame::reset_neighbours(){
    neighbours = false;
    return 0;
}

bool Frame::get_neighbours(){
    return neighbours;
}

bool Frame::get_coloured(){
    return coloured;
}

// Setting
void Frame::setx(double min, double max){
    xdim[0] = min;
    xdim[1] = max;
}
void Frame::sety(double min, double max){
    ydim[0] = min;
    ydim[1] = max;
}
void Frame::setz(double min, double max){
    zdim[0] = min;
    zdim[1] = max;
}

// Length
double Frame::xlen(){
    return 2*PI;//xdim[1] - xdim[0];
}
double Frame::ylen(){
    return 2*PI;//ydim[1] - ydim[0];
}
double Frame::zlen(){
    return zdim[1] - zdim[0];
}

double Frame::length(){
    double l = (xlength()+ylength())/2; 
    if (abs(l - xlength()) > EPS){
        cerr << "Box not square " << xlength() << " " << ylength() << endl;
    }
    return l;
}

void Frame::set_timestep(int t){
    timestep = t;
}

int Frame::get_timestep(){
    return timestep;
}

double Frame::xlength(){
    return (xmax() - xmin())/xlen();
}

double Frame::ylength(){
    return (ymax() - ymin())/ylen();
}

double Frame::zlength(){
    return (zmax() - zmin())/zlen();
}

// Min
double Frame::xmin(){
    return xdim[0];
}
double Frame::ymin(){
    return ydim[0];
}
double Frame::zmin(){
    return zdim[0];
}

// Max
double Frame::xmax(){
    return xdim[1];
}
double Frame::ymax(){
    return ydim[1];
}
double Frame::zmax(){
    return zdim[1];
}

