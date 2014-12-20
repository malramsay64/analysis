//
//  neighbours.cpp
//  analysis
//
//  Created by Malcolm Ramsay on 7/12/2014.
//  Copyright (c) 2014 Malcolm Ramsay. All rights reserved.
//

#include "neighbours.h"

using namespace std;

static float R_FACTOR = pow(2,1./6);

int par_neigh(Frame *frame){
    // Neighbours already calculated
    if (frame->get_neighbours()){
        return 0;
    }
    
    int per_thread = ceil(frame->num_atoms()/NUM_THREADS);
    vector<thread> threads;
    threads.reserve(NUM_THREADS);
    
    // Divide into threads
    int begin, end;
    for (int i = 0; i < NUM_THREADS; i++){
        begin = i*per_thread;
        if (i == NUM_THREADS-1){
            end = frame->num_atoms();
        }
        else{
            end = (i+1)*per_thread;
        }
        threads.push_back(thread(&loop_neigh, frame, begin, end));
    }
    
    // Join threads
    for (int j=0; j < NUM_THREADS; j++){
        threads.at(j).join();
    }
    frame->set_neighbours();
    return 0;
}

void *loop_neigh(Frame * frame, int begin, int end){
    // Loop through all particles in thread
    for (int i = begin; i < end && i < frame->num_atoms(); i++){
        find_neighbours(&frame->particles.at(i), frame);
    }
    return 0;
}


// This function looks for all nearest neighbours it updates the nearest neighbour list of both the particle that it found
// and itself.
// We are concerned about the position of both large and small particles, but add them to different lists for identification
void find_neighbours(particle *a, Frame *frame){
    double d;
    // starts looking at the atom after itself up to the end
    vector<particle>::iterator i;
    for (i = frame->particles.begin()+a->index()+1; i != frame->particles.end(); i++){
        particle *b = &(*i);
        // Check not in same molecule
        d = frame->dist(a->pos_vect(), b->pos_vect());
        // Bonded if the distance between them is less than
        // the sum of the radii multiplied by a factor
        if (d < R_FACTOR*(a->radius + b->radius)){
            // If different molecules
            if (a->molid != b->molid){
                //lock(frame->m.at(a->m_i()),frame->m.at(b->m_i()));
            }
            // In same molecule, only single lock needed
            else{
                //frame->m.at(a->m_i()).lock();
            }
            
            // Add particles as neighbours
            a->append(b);
            b->append(a);
            
            // Add Molecules as neighbours
            if (a->molid != b->molid){
                molecule *amol = &frame->molecules.at(a->m_i());
                molecule *bmol = &frame->molecules[b->molid-1];
                
                amol->add_neighbour(bmol);
                bmol->add_neighbour(amol);
                //frame->m.at(a->m_i()).unlock();
                //frame->m.at(b->m_i()).unlock();
            }
            else {
                //frame->m.at(a->m_i()).unlock();
            }
        }
    }
}

int short_range_order(Frame * frame){
    ofstream file;
    file.open("short_order.csv", ios::out);
    int colour = 0;
    vector<int> count(6,0);
    int total = 0;
    // Print m1
    vector<molecule>::iterator m1;
    for (m1 = frame->molecules.begin(); m1 != frame->molecules.end(); m1++){
        vect com1 = (*m1).COM();
        vect d;
        double theta;
        double rot = angle(&(*m1), frame);
        //file << (*m1).nump() << endl;
        for (int k = 0; k < (*m1).nump(); k++){
            //cout << k << endl;
            d = frame->direction((*m1).atom_pos(k),com1);
            theta = atan2(d) + PI - rot + PI/2;
            // In format for gnuplot: theta, dist, circle radius, colour
            file << theta << "," << d.length() << "," << (*m1).atoms[k]->radius << "," << 1 << endl;
        }
        
        for (int j = 0; j < (*m1).num_neighbours(); j++){
            // Multiple interactions between particles
            if ((*m1).nint[j] > 1){
                colour = 0;
                molecule *m2 = (*m1).my_neighbours[j];
                // Theta
                colour = order_type(&(*m1), m2, frame);
                count[colour]++;
                total++;
                if (colour){
                    // Order parameter data
                    vect com1 = (*m1).COM();
                    vect d;
                    double theta;
                    double rot = angle(&(*m1), frame);
                    for (int k = 0; k < m2->nump(); k++){
                        d = frame->direction(m2->atom_pos(k),com1);
                        theta = atan2(d) + PI - rot + PI/2;
                        file << theta << "," << d.length() << "," << 0.04 << "," << colour << endl;
                    }
                }
            }
        }
    }
    file.close();
    ofstream dist;
    dist.open("short_order_dist.csv", ios::out);
    for (int i = 0; i < count.size(); i++){
        dist << i << "," << count.at(i)/(double) total  << endl;
    }
    dist.close();
    return 0;
}

int order_type(molecule * m1, molecule * m2, Frame * frame){
    double Rb = (*m1).atoms[0]->radius;
    double Rh = (*m1).atoms[1]->radius;
    // Theta
    double theta = angle(m1, frame) - angle(m2, frame);
    // Body-Body Dist (body = atoms[0])
    double Dbb = frame->dist(m1->atom_pos(0), m2->atom_pos(0));
    // Body-Head Dist (head = atoms[1])
    double Dbh = min( frame->dist(m1->atom_pos(1), m2->atom_pos(0)), \
                     frame->dist(m1->atom_pos(0), m2->atom_pos(1)) );
    // Head-Head Dist
    double Dhh = frame->dist(m1->atom_pos(1), m2->atom_pos(1));
    
    double deltaD = 0.15;
    double deltaT = 5*PI/180;
    if (fabs( Dbh - (Rb + Rh) ) < deltaD){
        if (fabs( Dhh - (Rh + Rh) ) < deltaD){
            if (fabs( theta - PI )  < deltaT){
                // Anti Parallel 1
                return 3;
            }
            else if (fabs( Dbb - (Rb + Rb)) < deltaD){
                // Chiral
                return 5;
            }
        }
        else if (fabs( Dbb - (Rb + Rb)) < deltaD){
            if (fabs( theta - 0 ) < deltaT){
                // Parallel
                return 2;
            }
            else if (fabs( theta - PI ) < deltaT){
                // Anti Parallel 2
                return 4;
            }
        }
    }
    return 0;
}


/*
 void randomise_orientation(Frame * frame){
 ofstream file;
 file.open("trj/out.lammpstrj", ios::out);
 //print(frame, &file);
 // Initialise atoms
 int mol = 1;
 int num_mols = frame->num_mol();
 list<particle *> single;
 dyn_queue<particle> large;
 vector<particle *> small;
 vector<particle>::iterator i;
 // Set all particles to unallocated
 frame->reset_traverse();
 // Pick random first molecule
 int r = rand() % num_mols;
 particle *p = frame->molecules.at(r).atoms.front();
 large = dyn_queue<particle>(p);
 // While molecules unallocated
 while (p){
 // Take large particle from the front
 p = large.pop();
 cout << "Pop : " <<  p->id << endl;
 small = vector<particle *>(0);
 vector<particle *>::iterator it;
 // Find all neighbours
 for (it = p->my_neighbours.begin(); it != p->my_neighbours.end(); it++){
 if ((*it)->type == 2) {
 small.push_back(*it);
 }
 }
 r = rand() % small.size();
 p->molid = mol;
 small.at(r)->molid = mol;
 small.at(r)->traverse();
 frame->add_link(mol-1, p->index(), false);
 frame->add_link(mol-1, small.at(r)->index(), false);
 mol++;
 remove_neighbours(p,&single);
 remove_neighbours(small.at(r),&single);
 single.remove(small.at(r));
 //single.remove(p);
 //print(frame, &file);
 mol = check_single(frame, &single, &large, mol, &file);
 //cout << "Done " << mol << " of " << num_mols << endl;
 }
 vector<particle>::iterator z;
 for (z = frame->particles.begin(); z != frame->particles.end(); z++){
 //cout << "ID: " << (*z).id << " MolID: " << (*z).molid << endl;
 }
 //cout << "Neighbours" << endl;
 file.close();
 frame->reset_neighbours();
 par_neigh(frame);
 //cout << "Exit" << endl;
 }
 
 int check_single(Frame * frame, list<particle *> *single, dyn_queue<particle> *large, int mol, ofstream *file){
 particle *p;
 //cout << "Check " << single->size() << endl;
 while (!single->empty()){
 vector<particle *> neigh(0);
 p = single->front();
 single->pop_front();
 if (p->molid == 0 && p->type == 2){
 vector<particle *>::iterator it;
 for (it = p->my_neighbours.begin(); it != p->my_neighbours.end(); it++){
 if ((*it)->type == 1 && (*it)->molid == 0){
 neigh.push_back(*it);
 }
 }
 if (neigh.size()){
 p->molid = mol;
 neigh.front()->molid = mol;
 p->traverse();
 neigh.front()->traverse();
 large->remove(neigh.front());
 cout << "Allocate Single " <<  p->id << " " << neigh.front()->id << endl;
 frame->add_link(mol-1, neigh.front()->index(), false);
 frame->add_link(mol-1, p->index(), false);
 mol++;
 remove_neighbours(p,single);
 remove_neighbours(neigh.front(),single);
 
 }
 }
 if (p->molid == 0 && p->type == 1){
 vector<particle *>::iterator it;
 for (it = p->my_neighbours.begin(); it != p->my_neighbours.end(); it++){
 if ((*it)->type == 2 && (*it)->molid == 0){
 neigh.push_back(*it);
 }
 else if ((*it)->type == 1 && (*it)->molid == 0){
 cout << "Push: " << (*it)->id << endl;
 large->push(*it);
 }
 }
 if (neigh.size()){
 p->molid = mol;
 neigh.front()->molid = mol;
 cout << "Allocate large Single " <<  p->id << " " << neigh.front()->id << endl;
 frame->add_link(mol-1, neigh.front()->index(), false);
 frame->add_link(mol-1, p->index(), false);
 mol++;
 remove_neighbours(neigh.front(),single);
 remove_neighbours(p,single);
 
 }
 }
 print(frame, file);
 }
 return mol;
 }
 
 int remove_neighbours(particle *p, list<particle *> *single){
 vector<particle *>::iterator it1,it2;
 //cout << "Neighbours: " << p->numn() << endl;
 for (it1 = p->my_neighbours.begin(); it1 != p->my_neighbours.end(); it1++){
 for (it2 = (*it1)->my_neighbours.begin(); it2 != (*it1)->my_neighbours.end(); it2++){
 if ((*it2)->id == p->id){
 it2--;
 (*it1)->my_neighbours.erase(it2+1);
 }
 }
 if ((*it1)->molid == 0){
 if ((*it1)->n_large() == 1 && (*it1)->type == 2){
 cout << "Add Single" << (*it1)->id << endl;
 single->push_back(*it1);
 }
 else if ((*it1)->n_small() == 1 && (*it1)->type == 1){
 //cout << "Add Single" << (*it1)->id << endl;
 //single->push_back(*it1);
 }
 }
 }
 p->my_neighbours.clear();
 return 0;
 }
 */

