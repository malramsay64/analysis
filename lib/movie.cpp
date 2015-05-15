//
//  movie.cpp
//  analysis
//
//  Created by Malcolm Ramsay on 9/02/2015.
//  Copyright (c) 2015 Malcolm Ramsay. All rights reserved.
//

#include "movie.h"

using namespace std;

double movie_colour(molecule * mol, Frame * frame){
    double mol_colour = 0;
    return mol->num_contacts();
    vector<int> order = short_neighbour_list(mol, frame);
    // Has Antiparallel ordering
    for (auto j: order){
        if (j == 4){
            mol_colour = 1;
        }
    }
    return mol_colour;
}

int print_movie(ofstream * file, Frame * frame, molecule * mol){
    if (mol == 0){
        *file << "ITEM: TIMESTEP" << endl << frame->timestep << endl;
        *file << "ITEM: NUMBER OF ATOMS" << endl << frame->num_atoms() << endl;
        *file << "ITEM: BOX BOUNDS xy xz yz pp pp pp" << endl;
        *file << 0 << " " << frame->get_a() << " " << frame->get_tilt() << endl;
        *file << 0 << " " << frame->get_height() << " " << 0 << endl;
        *file << "-0.5 0.5 0" << endl;
        *file << "ITEM: ATOMS id mol type x y z vx vy vz" << endl;
    }
    else {
        vect cart_com, d;
        cart_com = frame->cartesian(mol->COM());
        cart_com = wrap_x(cart_com, frame->get_a());
        double mol_colour = movie_colour(mol, frame);
        std::vector<particle *>::iterator i;
        for (i = mol->atoms.begin(); i != mol->atoms.end(); i++){
            d = frame->cartesian(direction(mol->COM(), (*i)->pos_vect()));
            *file << (*i)->id << " " << mol->id << " " << (*i)->type << " " << cart_com + d << " 0 " << mol_colour << " 0 0" << endl;
        }
    }
    return 0;
}
