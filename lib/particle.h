#include <iostream>
#include <vector>
#include <list>
#include <algorithm>
#include <math.h>
#include <fstream>
#include "functions.h"
#include <pthread.h>
#include "vect.h"

#ifndef PI
#define PI 3.1415265358978323
#endif

#define MAX_MOL_CONTACTS 200
#define MAX_ATOM_CONTACTS 100
#define MOL_SIZE 3

#define XX 0
#define YY 1

#define EPS 5e-2

#define ORDER_LEN 15

#ifndef PARTICLE_H
#define PARTICLE_H



// Particle class
class particle {
    bool traversed;
  public:
    int id;
    int molid;
    int type;
    double radius;
    double mass;
    double pos[2];
    std::vector<particle *> my_neighbours;

    particle();
    //~particle();
    void append(particle *);
    int numn();
    int n_large();
    int n_small();
    double xpos();
    double ypos();
    int set_xpos(double);
    int set_ypos(double);
    int set_pos(vect);
    vect pos_vect();
    int index();
    int m_i();
    bool get_traversed();
    void traverse();
    void reset_traverse();

    bool operator> (const particle &b);
    bool operator>= (const particle &b);
    
    bool operator< (const particle &b);
    bool operator<= (const particle &b);
 
    bool operator== (const particle &b);
    bool operator!= (const particle &b);
} ;

bool operator < (const particle &a, const particle &b);

class molecule {
    int colour;
    bool traversed;
  public:
    int type;
    int id;
    std::vector<particle *> atoms;    
    std::vector<molecule *> my_neighbours;    
    int contacts;
    std::vector<int> nint;

    molecule();
    //~molecule();
    void add_neighbour(molecule *);
    int uniqc();
    int nump();
    double mass();
    vect COM();
    vect atom_pos(int);
    int num_contacts();
    int num_neighbours();
    int pairing(molecule *);
    int get_colour();
    int set_colour(int);
    int graph_colour();
    bool get_traversed();
    void traverse();
    void reset_traverse();
    int same_period();
} ;

#endif
