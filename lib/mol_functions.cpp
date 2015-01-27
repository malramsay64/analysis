//
//  mol_functions.cpp
//  analysis
//
//  Created by Malcolm Ramsay on 7/12/2014.
//  Copyright (c) 2014 Malcolm Ramsay. All rights reserved.
//

#include "mol_functions.h"

using namespace std;

vect orientation(molecule *m, Frame *frame){
    vect com = m->COM();
    vect v = frame->direction(com, m->atom_pos(0));
    if (m->nump() == 2 || (v.length() < EPS)){
        v = frame->direction(com, m->atoms.at(1)->pos_vect());
        v.orthogonalise();
    }
    v.normalise();
    return v;
}

double angle(molecule *m, Frame * frame){
    return atan2(orientation(m, frame)) + PI;
}

std::ostream& print_mol(std::ostream &os, molecule *mol, Frame *frame){
    vect d, com;
    
    com = frame->cartesian(mol->COM());
    com = wrap(com, frame->get_a(), frame->get_height());
    std::vector<particle *>::iterator i;

    for (i = mol->atoms.begin(); i != mol->atoms.end(); i++){
        d = frame->cartesian(direction(mol->COM(), (*i)->pos_vect()));
        os << com + d << " " << (*i)->radius << " " << mol->get_colour() << " " << mol->id << std::endl;
    }
    
    return os;
}


