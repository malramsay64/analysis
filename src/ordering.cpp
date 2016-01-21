//
//  mol_functions.cpp
//  analysis
//
//  Created by Malcolm Ramsay on 7/12/2014.
//  Copyright (c) 2014 Malcolm Ramsay. All rights reserved.
//

#include "ordering.h"

using namespace std;

Vector2d orientation(molecule *m, Frame *frame){
    Vector2d com = m->COM();
    if (m->nump() == 1){
        return Vector2d(0, 0);
    }
    Vector2d v = frame->direction(com, m->atom_pos(0));
    if (m->nump() == 2 || (v.length() < EPS)){
        v = frame->direction(com, m->get_large()->pos_vect());
        v.orthogonalise();
    }
    v.normalise();
    return v;
}

double angle(molecule *m, Frame * frame){
    return atan2(orientation(m, frame)) + PI;
}

Vector2d wrap_x(Vector2d v, double a){
    double x = v.x/a;
    x = x*2*PI;
    x = atan2(sin(x),cos(x))+PI;
    x = x/(2*PI);
    x = x*a;
    return Vector2d(x, v.y);
}

Vector2d wrap(Vector2d v){
    return atan2(sin(v),cos(v))+PI;
}

double mol_colour(molecule * m, Frame * frame){
    //return circle_colour(m);
    return my_mod(angle(m, frame),PI)/PI;
}


double com_colour(molecule * m, Frame * frame){
    return hexatic(6, m, frame);
}

int circle_colour(molecule * m){
    int colour = 0;
    for (auto &p: m->atoms){
        colour += p->order();
    }
    return colour;
}

int neighbour_colour(molecule * m, Frame *frame){
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
        else {
            relax.add(1);
        }
    }
    return relax.get_mean();
}

double hexatic(int n, molecule* m1, Frame *frame){
    double theta;
    complex<double> sum = complex<double>(0,0);
    complex<double> i = complex<double>(0,1);
    for (auto &m2: m1->my_neighbours){
        theta = frame->direction(m2->COM(), m1->COM()).angle();
        sum += (1./m1->num_neighbours())*exp(6*theta*i);
    }
    return abs(sum);
}

double circle_ordering(molecule *m){
    my_mean order;
    for (auto &p: m->atoms){
        order.add(p->order());
    }
    return order.get_mean();
}

double orient_ordering(molecule *m){
    my_mean mean;
    for (auto m2: m->my_neighbours){
        mean.add(pow(dot_product(m->get_orient_vect(), m2->get_orient_vect()),2));
    }
    return mean.get_mean();
}

int short_ordering(molecule *mol, Frame * frame){
    int type = 0;
    for (auto mol2: mol->my_neighbours){
        type = order_type(mol, mol2, frame);
        if (type == 2 || type == 3 || type == 4 || type == 5){
            return type;
        }
    }
    return type;
}

int tri_ordering(molecule *mol){
    if (mol->max_pairing() == 4){
        for (int i = 0; i < mol->num_neighbours(); i++){
            if (mol->nint.at(i) == 4){
                if (fabs(fmod(mol->my_neighbours.at(i)->get_orientation() - mol->get_orientation(), PI)) < 0.1 ){
                    return 1;
                }
            }
        }
    }
    return 0;
}

molecule reorient(molecule *m, Frame* frame){
    Vector2d d;
    
    double delta_t = m->get_orientation();
    molecule n = molecule(*m);
    for (auto p: n.atoms){
        d = frame->direction(p->pos_vect(), n.COM());
        p->set_pos(Vector2d(d.length() * sin(d.angle() - delta_t), d.length() * cos(d.angle() - delta_t)));
    }
    return n;
}
