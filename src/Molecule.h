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

using namespace LAlgebra;

class Molecule {
    double rotation{0};
    double orientation{0};
    int colour{0};
    Vector com{};
public:
    int type{0};
    int id{0};
    std::vector<Particle *> atoms{};
    std::map<Molecule *, int> my_neighbours{};
    int contacts{0};
    
    Molecule() {};
    void add_neighbour(Molecule *);
    void delete_neighbours();
    void delete_mol_neighbours();
    size_t uniqe_contacts() const;
    size_t num_particles() const ;
    int num_contacts() const ;
    int num_neighbours() const ;
    double get_mass() const;
    Particle * get_large() const;
    Vector get_COM();
    Vector get_COM() const;
    Vector moved_COM() const ;
    Vector calc_COM();
    Vector update_COM();
    Vector atom_pos(int) const ;
    std::map<int, int> pairing() const ;
    int max_pairing() const ;
    int same_period();
    int index() const ;
    double update_orientation(double);
    double get_orientation() const ;
    int set_orientation(double);
    double get_rotation() const ;
    Vector get_orient_vect() const;
    std::vector<Molecule *> get_neighbours();
} ;

inline bool operator> (const Molecule &lhs, const Molecule &rhs) { return lhs.id > rhs.id; };
inline bool operator>= (const Molecule &lhs, const Molecule &rhs) { return lhs.id >= rhs.id; };
inline bool operator< (const Molecule &lhs, const Molecule &rhs) { return lhs.id < rhs.id; };
inline bool operator<= (const Molecule &lhs, const Molecule &rhs) { return lhs.id <= rhs.id; };
inline bool operator== (const Molecule &lhs, const Molecule &rhs) { return lhs.id == rhs.id; }
inline bool operator!= (const Molecule &lhs, const Molecule &rhs) { return lhs.id != rhs.id; };

#endif /* defined(MY_MOLECULE) */
