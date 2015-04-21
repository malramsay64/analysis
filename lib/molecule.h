//
//  molecule.h
//  analysis
//
//  Created by Malcolm Ramsay on 7/12/2014.
//  Copyright (c) 2014 Malcolm Ramsay. All rights reserved.
//

#include "particle.h"
#include "distribution.h"

#ifndef MY_MOLECULE
#define MY_MOLECULE


class molecule {
    double rotation;
    double orientation;
    int colour;
    vect com;
public:
    int type;
    int id;
    std::vector<particle *> atoms;
    std::vector<molecule *> my_neighbours;
    int contacts;
    std::vector<int> nint;
    
    molecule();
    void add_neighbour(molecule *);
    void delete_neighbours();
    int uniqc();
    int nump();
    double mass();
    vect COM();
    vect calc_COM();
    vect update_COM();
    vect atom_pos(int);
    int num_contacts();
    int num_neighbours();
    distribution<int> pairing();
    int max_pairing();
    int get_colour();
    int set_colour(int);
    int graph_colour();
    int same_period();
    int index();
    double update_orientation(double);
    double get_orientation();
    int set_orientation(double);
    double get_rotation();
    vect get_orient_vect();
    
} ;


#endif /* defined(MY_MOLECULE) */
