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

int mol_colour(molecule * m, Frame * frame){
    vector<int> list;
    list = short_neighbour_list(m, frame);
    int max = 0;
    for (auto i: list){
        if ( i > max){
            max = i;
        }
    }
    return max;
}


double struct_relax(molecule * m, Frame * init){
    my_mean relax;
    molecule * m2 = &init->molecules.at(m->index());
    for (int i = 0; i < m->nump(); i++){
        if (init->dist(m->atom_pos(i), m2->atom_pos(i)) > STRUCT_DIST){
            relax.add(0);
        }
        else
            relax.add(1);
    }
    return relax.get_mean();
}

double hexatic(int n, molecule* m1, Frame *frame){
    my_mean mean;
    complex<double> i = sqrt(-1);
    double theta;
    for (auto &m2: m1->my_neighbours){
        theta = frame->direction(m2->COM(), m1->COM()).angle();
        mean.add(exp(n*theta*i).real());
    }
    return mean.get_mean();
}

