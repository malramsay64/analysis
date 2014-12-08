//
//  molecule.h
//  analysis
//
//  Created by Malcolm Ramsay on 7/12/2014.
//  Copyright (c) 2014 Malcolm Ramsay. All rights reserved.
//

#include "particle.h"

#ifndef MY_MOLECULE
#define MY_MOLECULE


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


#endif /* defined(MY_MOLECULE) */
