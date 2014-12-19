//
//  parallel.cpp
//  analysis
//
//  Created by Malcolm Ramsay on 7/12/2014.
//  Copyright (c) 2014 Malcolm Ramsay. All rights reserved.
//

#include "parallel.h"


using namespace std;


template <class type> my_mean par_mean(vector<type> *list, double (*f)(type *, Frame *), Frame *frame){
    int size = (int) (list->size()/NUM_THREADS);
    vector<thread> threads;
    int begin, end;
    vector<my_mean>means(NUM_THREADS);
    for (int i = 0; i < NUM_THREADS; i++){
        begin = i*size;
        end = (i+1)*size;
        threads.push_back(thread(&mean_single_thread<type>, list, begin, end, f, frame, &means[i]));
    }
    my_mean mean;
    for (int j = 0; j < NUM_THREADS; j++){
        threads.at(j).join();
        mean.combine(means.at(j));
    }
    return mean;
}


template <class type> void * mean_single_thread(vector<type> *l1, int begin, int end, double (*f)(type *, Frame *), Frame *frame, my_mean *res){
    double val;
    typename vector<type>::iterator i1;
    for (i1 = l1->begin()+begin; i1 != l1->begin()+end && i1 != l1->end(); i1++){
        val = (*f)(&(*i1), frame);
        res->add(val);
    }
    return 0;
}


template <class type> my_mean par_mean(vector<type> *list1, vector<type> *list2, vect (*f)(type *, Frame *), double (*diff)(vect, vect, int), Frame *frame, int l){
    int size = ceil(list1->size()/NUM_THREADS);
    vector<thread> threads;
    threads.reserve(NUM_THREADS);
    int begin, end;
    vector<my_mean> means(NUM_THREADS);
    for (int i = 0; i < NUM_THREADS; i++){
        begin = i*size;
        end = (i+1)*size;
        threads.push_back(thread(&mean_thread<type>, list1, list2, begin, end, f, diff, frame, l, &means[i]));
    }
    my_mean mean;
    for (int j = 0; j < NUM_THREADS; j++){
        threads.at(j).join();
        // Last thread
        mean.combine(means.at(j));
    }
    return mean;
}

template <class type> void * mean_thread(vector<type> *l1, vector<type> *l2, int begin, int end, vect (*f)(type*, Frame *), double (*diff)(vect, vect, int), Frame *frame, int l, my_mean *res){
    vect x1,x2;
    double val;
    typename vector<type>::iterator i1, i2;
    for (i1 = l1->begin()+begin, i2 = l2->begin()+begin;\
         i1 != l1->begin()+end && i1 != l1->end() &&\
         i2 != l2->begin()+end && i1 != l2->end();\
         i1++, i2++){
        
        x1 = (*f)(&(*i1), frame);
        x2 = (*f)(&(*i2), frame);
        val = (diff)(x1,x2,l);
        res->add(val);
    }
    return 0;
}

double v_angle(vect v1, vect v2, int l){
    return dist((atan2(v1)),(atan2(v2)));
}

double v_legendre(vect v1, vect v2, int l){
    return legendre(l, dot_product(v1,v2));
}

double v_dist(vect v1, vect v2, int l){
    return dist(v1,v2);
}

double dist_sq(vect v1, vect v2, int l){
    return pow(dist(v1,v2),2);
}

vect COM(molecule *m, Frame * frame){
    return m->COM()*frame->size();
}

double par_rotational_relaxation(int l, Frame *init, Frame *curr){
    return par_mean<molecule>(&init->molecules, &curr->molecules, &orientation, &v_legendre, init, l).get_mean();
}

double par_average_moved(Frame *init, Frame *curr){
    return par_mean<molecule>(&init->molecules, &curr->molecules, &COM, &v_dist, init, 0).get_mean();
}

double par_MSD(Frame *init, Frame *curr){
    return par_mean<molecule>(&init->molecules, &curr->molecules, &COM, &dist_sq, init, 0).get_mean();
}

double par_diffusion_constant(Frame *init, Frame *curr){
    int steps = curr->timestep - init->timestep;
    return par_MSD(init, curr)/(4*steps*DELTA_T);
}

double bonded_frac(Frame *frame){
    return par_mean<particle>(&frame->particles, &average_bonded_frac, frame).get_mean();
}

double par_average_rotation(Frame *init, Frame *curr){
    return par_mean<molecule>(&init->molecules, &curr->molecules, &orientation, &v_angle, init, 0).get_mean();
}

double par_local_order(Frame *frame){
    return par_mean<molecule>(&frame->molecules, &local_order, frame).get_mean();
}

double par_global_order(Frame * frame){
    return par_mean<molecule>(&frame->molecules, &global_order, frame).get_mean();
}

double par_circle_order(Frame * frame){
    return par_mean<particle>(&frame->particles, &circle_order, frame).get_mean();
}

int order_parameter(int reference, ofstream *file, Frame *frame){
    par_neigh(frame);
    vector<molecule>::iterator mol;
    for (mol = frame->molecules.begin(); mol != frame->molecules.end();  mol++) {
        vect com = mol->COM();
        vect v,d;
        double rot = angle(&(*mol), frame);
        double theta;
        // For each particle
        vector<particle>::iterator atom;
        for (atom = frame->particles.begin(); atom != frame->particles.end(); atom++) {
            v = atom->pos_vect();
            d = frame->direction(v,com);
            theta = atan2(d) + PI - rot + PI/2;
            if (d.length() < ORDER_LEN){
                *file << theta << "," << d.length() << "," << atom->type << endl;
            }
        }
        // Center of Mass
        vector<molecule>::iterator m;
        for (m = frame->molecules.begin(); m != frame->molecules.end(); m++){
            v = m->COM();
            d = frame->direction(v,com);
            theta = atan2(d) + PI - rot + PI/2;
            //orientation = angle(&(*m), frame) - angle(&(*mol),frame);
            if (d.length() < ORDER_LEN){
                *file << theta << "," << d.length() << "," << 3 << ","  << endl;
            }
        }
    }
    return 0;
}




