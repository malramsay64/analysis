//
//  mol_functions.cpp
//  analysis
//
//  Created by Malcolm Ramsay on 7/12/2014.
//  Copyright (c) 2014 Malcolm Ramsay. All rights reserved.
//

#include "ordering.h"

using namespace std;

Vector2d orientation(const Molecule &m, const Frame &frame){
    Vector2d com = m.get_COM();
    if (m.num_particles() == 1){
        return Vector2d{0, 0};
    }
    Vector2d v = direction(com, m.atom_pos(0), frame);
    if (m.num_particles() == 2 || (v.length() < EPS)){
        v = direction(com, m.get_large()->pos_vect(), frame);
        v.orthogonalise();
    }
    v.normalise();
    return v;
}

double angle(const Molecule &m, const Frame &frame){
    return atan2(orientation(m, frame)) + PI;
}

Vector2d wrap_x(const Vector2d &v, double a){
    double x = v.x/a;
    x = x*2*PI;
    x = atan2(sin(x),cos(x))+PI;
    x = x/(2*PI);
    x = x*a;
    return Vector2d(x, v.y);
}

Vector2d wrap(const Vector2d &v){
    return atan2(sin(v),cos(v))+PI;
}

double mol_colour(const Molecule &m, const Frame &frame){
    return my_mod(angle(m, frame),PI)/PI;
}


double com_colour(const Molecule &m, const Frame &frame){
    return hexatic(6, m, frame);
}

int circle_colour(const Molecule &m){
    int colour = 0;
    for (auto &p: m.atoms){
        colour += p->order();
    }
    return colour;
}

int neighbour_colour(const Molecule &m, const Frame &frame){
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

double struct_relax(const Molecule &m, Frame &init){
    Stats::Stats relax;
    Molecule &m2 = init.molecules.at(m.index());
    for (int i = 0; i < m.num_particles(); i++){
        if (dist(m.atom_pos(i), m2.atom_pos(i), init) > STRUCT_DIST){
            relax.push(0);
        }
        else {
            relax.push(1);
        }
    }
    return relax.getMean();
}

double hexatic(int n, const Molecule &m1, const Frame &frame){
    double theta;
    complex<double> sum = complex<double>(0,0);
    complex<double> i = complex<double>(0,1);
    for (auto &m2: m1.my_neighbours){
        theta = direction(m2.first->get_COM(), m1.get_COM(), frame).angle();
        sum += (1./m1.num_neighbours())*exp(6*theta*i);
    }
    return abs(sum);
}

double circle_ordering(const Molecule &m){
    Stats::Stats order;
    for (auto &p: m.atoms){
        order.push(p->order());
    }
    return order.getMean();
}

double orient_ordering(const Molecule &m){
    Stats::Stats mean;
    for (auto m2: m.my_neighbours){
        mean.push(pow(dot_product(m.get_orient_vect(), m2.first->get_orient_vect()),2));
    }
    return mean.getMean();
}

int short_ordering(const Molecule &mol, const Frame &frame){
    int type = 0;
    for (auto mol2: mol.my_neighbours){
        type = order_type(mol, *mol2.first, frame);
        if (type == 2 || type == 3 || type == 4 || type == 5){
            return type;
        }
    }
    return type;
}

int tri_ordering(const Molecule &mol){
    if (mol.max_pairing() == 4){
        for (auto neigh: mol.my_neighbours) {
            if(neigh.second == 4){
                if (fabs(fmod(neigh.first->get_orientation() - mol.get_orientation(), PI)) < 0.1 ){
                    return 1;
                }
            }
        }
    }
    return 0;
}

Molecule reorient(const Molecule &m, const Frame &frame){
    Vector2d d{};
    double delta_t = m.get_orientation();
    Molecule n{m};
    for (auto p: n.atoms){
        d = direction(p->pos_vect(), n.get_COM(), frame);
        p->set_pos(Vector2d(d.length() * sin(d.angle() - delta_t), d.length() * cos(d.angle() - delta_t)));
    }
    return n;
}
