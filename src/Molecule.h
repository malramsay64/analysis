//
//  molecule.h
//  analysis
//
//  Created by Malcolm Ramsay on 7/12/2014.
//  Copyright (c) 2014 Malcolm Ramsay. All rights reserved.
//

#ifndef MY_MOLECULE
#define MY_MOLECULE

#include "Particle.h"
#include <Stats.h>
#include <map>

class Molecule {
    double rotation{0};
    double orientation{0};
    int colour{0};
    Vector2d com{};
public:
    int type{0};
    int id{0};
    std::vector<Particle *> atoms{};
    std::map<Molecule *, int> my_neighbours{};
    int contacts{0};

    Molecule();
    void add_neighbour(Molecule *);
    void delete_neighbours();
    size_t uniqe_contacts() const;
    size_t num_particles() const ;
    double get_mass() const;
    Particle * get_large();
    Vector2d get_COM();
    Vector2d moved_COM() const ;
    Vector2d calc_COM();
    Vector2d update_COM();
    Vector2d atom_pos(int) const ;
    int num_contacts() const ;
    int num_neighbours() const ;
    std::map<int, int> pairing() const ;
    int max_pairing() const ;
    int get_colour() const;
    int set_colour(int);
    int graph_colour();
    int same_period();
    int index() const ;
    double update_orientation(double);
    double get_orientation() const ;
    int set_orientation(double);
    double get_rotation() const ;
    Vector2d get_orient_vect() const;
    void delete_mol_neighbours();
} ;


#endif /* defined(MY_MOLECULE) */
