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

ostream& print_mol(ostream &os, molecule *mol, Frame *frame){
    vect d, com;
    
    com = frame->cartesian(mol->COM());
    com = wrap_x(com, frame->get_a());
    
    vector<particle *>::iterator i;

    for (i = mol->atoms.begin(); i != mol->atoms.end(); i++){
        d = frame->cartesian(direction(mol->COM(), (*i)->pos_vect()));
        os << com + d << " " << (*i)->radius << " " << mol->get_colour() << " " << mol->id << endl;
    }
    
    return os;
}

vect wrap_x(vect v, double a){
    double x = v.x/a;
    x = x*2*PI;
    x = atan2(sin(x),cos(x))+PI;
    x = x/(2*PI);
    x = x*a;
    return vect(x,v.y);
}

vect wrap(vect v){
    return atan2(sin(v),cos(v))+PI;
}



