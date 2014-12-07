#include <iostream>
#include <fstream>
#include <string>
#include <exception>
#include <string.h>
#include <stdlib.h>
#include <cmath>
#include "particle.h"
#include <assert.h>
#include <mutex>
//#include "parallel.h"

#ifndef FRAME_H
#define FRAME_H

#define DELTA_T 0.005

#define cilk_for for

class Frame {
    int atoms;
    bool coloured;
    bool neighbours;
    double xdim[2];
    double ydim[2];
    double zdim[2];

  public:
    std::vector<particle> particles;
    std::vector<molecule> molecules;
    int timestep;
    int num_molecules;
    std::vector<std::mutex> m;
    

    Frame(); 
    //Frame(std::ifstream *);
    //~Frame();

    double par_rotational_relaxation(int l, Frame *init);
    double rotational_relaxation(int, Frame * init);
    double MSD(int i, Frame *init);
    double MSD(Frame *init);
    double diffusion_constant(Frame *init);
    double average_moved(Frame *);
    double par_average_moved(Frame *init);
    double average_bonded_frac();
    double average_rotation(Frame *);
    void set_timestep(int);
    int get_timestep();
    void set_atoms(int);
    int num_atoms();
    void set_num_mol(int);
    int num_mol();
    void add_particle(particle *p);
    void add_link(int, int, bool);
    void add_link(int, int);
    double dist(vect,vect);
    vect direction(vect,vect);
    vect size();
    int reset_traverse();
    int set_coloured();
    int set_neighbours();
    int reset_neighbours();
    bool get_neighbours();
    bool get_coloured();

    double num_contacts(std::vector<int> *);
    double num_neighbours(std::vector<int> *);
    double pairing(std::vector<int> *);
    

    void setx(double, double);
    void sety(double, double);
    void setz(double, double);
    
    double xlen();
    double ylen();
    double zlen();

    double length();
    double xlength();
    double ylength();
    double zlength();
    
    double xmin();
    double xmax();
    double ymin();
    double ymax();
    double zmin();
    double zmax();
};

#endif
