//
//  mol_functions.cpp
//  analysis
//
//  Created by Malcolm Ramsay on 7/12/2014.
//  Copyright (c) 2014 Malcolm Ramsay. All rights reserved.
//

#include "mol_functions.h"

using namespace std;


double average_bonded_frac(particle *p, Frame * frame){
    double mean = 0, stdev = 0, mprev, d, frac;
    vect v;
    vector<particle *>::iterator i;
    int num = 1;
    for (i = p->my_neighbours.begin(); i != p->my_neighbours.end(); i++){
        v = frame->direction(p->pos_vect(), (*i)->pos_vect());
        d = v.length();
        frac = d/(p->radius + (*i)->radius);
        mprev = mean;
        mean += (frac - mean)/(num);
        stdev += (frac-mean) * (frac - mprev);
        num++;
    }
    if (mean == 0){
        return 1;
    }
    return mean;
}

vect orientation(molecule *m, Frame *frame){
    vect com = m->COM();
    vect v = frame->direction(com, m->atoms.front()->pos_vect());
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

int get_colour(particle *p, Frame *frame){
    return frame->molecules.at(p->m_i()).get_colour();
}

int graph_colour(particle *p, Frame *frame){
    return frame->molecules.at(p->m_i()).graph_colour();
}


angle_list like_me(Frame *frame){
    vector<molecule>::iterator mol;
    angle_list angles;
    vect com1, com2, d;
    for (mol = frame->molecules.begin(); mol != frame->molecules.end(); mol++){
        dyn_queue<molecule> queue;
        molecule * m;
        queue.push(&(*mol), 0);
        queue.pop();
        m = queue.pop();
        com1 = (*mol).COM();
        int c = 0;
        while (queue.get_depth() < 4 && m){
            if (m->get_colour() == (*mol).get_colour()){
                com2 = m->COM();
                d = direction(com1,com2);
                angles.push(atan2(d), (d*frame->size()).length());
                c++;
            }
            m = queue.pop();
        }
    }
    return angles;
}


int colourise(Frame * frame){
    if (frame->get_coloured()){
        return 0;
    }
    angle_list a;
    vector<molecule>::iterator m;
    for (m = frame->molecules.begin(); m != frame->molecules.end(); m++){
        (*m).set_colour(a.push(angle(&(*m),frame)));
    }
    frame->set_coloured();
    return 0;
}

int graph_colourise(Frame *frame){
    // Frame already coloured
    if (frame->get_coloured()){
        return 0;
    }
    dyn_queue<particle> q = dyn_queue<particle>(&(frame->particles.front()));
    particle *p = q.pop();
    while (p){
        graph_colour(p, frame);
        p = q.pop();
    }
    frame->set_coloured();
    return 0;
}

double local_order(molecule * m, Frame * frame){
    vector<molecule *>::iterator i;
    double sum = 0;
    for (i = m->my_neighbours.begin(); i != m->my_neighbours.end(); i++){
        sum += pow(dot_product(orientation(m,frame),orientation((*i), frame)),2);
    }
    return sum/m->num_neighbours();
}

double global_order(molecule * m, Frame * frame){
    vector<molecule>::iterator i;
    my_mean mean;
    for (i = frame->molecules.begin(); i != frame->molecules.end(); i++){
        mean.add(pow(dot_product(orientation(m,frame),orientation(&(*i), frame)),2));
    }
    return mean.get_mean();
}

double circle_order(particle * p, Frame * frame){
    int large = p->n_large();
    int small = p->n_small();
    if (p->type == 1){
        if (large == 3 && small == 4){
            return 1;
        }
    }
    else{
        if (large == 4 && small == 1){
            return 1;
        }
    }
    return 0;
}

