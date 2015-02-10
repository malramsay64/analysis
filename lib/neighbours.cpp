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
double deltaD = 0.30;
double deltaT = 10*PI/180;
static int short_order_types = 7;

static double neigh_cutoff = 5.0;
static double neigh_alert = 7.0;

int check_particles(molecule * mol1, molecule * mol2, Frame * frame){
    double d_atom;
    for (auto &p2: mol2->atoms){
        for (auto &p1: mol1->atoms){
            d_atom = frame->dist(p1->pos_vect(), p2->pos_vect());
            d_atom /= p1->radius + p2->radius;
            if (d_atom < R_FACTOR){
                // Add Neighbours
                add_mol_neighbours(mol1, mol2);
                add_part_neighbours(p1,p2);
            }
        }
    }
    return 0;
}

bool find_mol_neighbours(molecule * mol, Frame * frame, vector<vector<int>> *neigh_list){
    vect com;
    vector<int> * neighbours = &neigh_list->at(mol->index());
    // mol properties
    com = mol->COM();
    bool recompute = false;
    
    // Check no neighbours in list
    if (neighbours->size() == 0){
        // Iterate through all molecules with id > mol
        vector<molecule>::iterator mol2;
        for (mol2 = frame->molecules.begin()+mol->id; mol2 != frame->molecules.end(); mol2++){
            double d_com = frame->dist(com, (*mol2).COM());
            if (d_com < neigh_cutoff){
                neighbours->push_back(mol2->index());
                // Check interparticle distnces
                check_particles(mol, &(*mol2), frame);
            }
        }
    }
    // Iterate through neighbour list
    else {
        for (auto &m: *neighbours){
            molecule *mol2 = &(frame->molecules.at(m));
            
            double d_com = frame->dist(com, mol2->COM());
            if (d_com < neigh_cutoff){
                check_particles(mol, mol2, frame);
            }
            else if (d_com > neigh_alert){
                // Flag local neighbour list for recomputation
                recompute = true;
            }
        }
    }
    return recompute;
}

int recompute_neighbours(molecule * mol, Frame * frame, vector<vector<int>> *neigh_list){
    // Recompute local neighbour list
    // Reset particle neighbour list
    neigh_list->at(mol->index()) = vector<int>();
    vector<int> * neighbours;
    neighbours = &neigh_list->at(mol->index());
    vect com = mol->COM();
    molecule *mol2;
    dyn_queue queue = dyn_queue(mol);
    mol2 = queue.pop();
    while (queue.get_depth() < 2 && mol2){
        mol2 = queue.pop();
        if (mol2->id > mol->id){
            double d_com = frame->dist(com, mol2->COM());
            if (d_com < neigh_cutoff){
                neighbours->push_back(mol2->index());
                check_particles(mol, mol2, frame);
            }
        }
    }
    // All neighbours
    for (auto &m: *neighbours){
        molecule *m1 = &(frame->molecules.at(m-1));
        neigh_list->at(m1->id-1) = vector<int>();
        neighbours = &neigh_list->at(m1->id-1);
        dyn_queue queue = dyn_queue(m1);
        mol2 = queue.pop();
        while (queue.get_depth() < 2 && mol2){
            mol2 = queue.pop();
            if (mol2 > mol){
                double d_com = frame->dist(m1->COM(), mol2->COM());
                if (d_com < neigh_cutoff){
                    neighbours->push_back(mol2->id);
                    check_particles(m1, mol2, frame);
                }
            }
        }
    }
    return 0;
}


bool check_mol_neighbours(molecule *m1, molecule * m2, Frame * frame){
    vector<particle *>::iterator p1, p2;
    double d;
    for (p1 = m1->atoms.begin(); p1 != m1->atoms.end(); p1++){
        for (p2 = m2->atoms.begin(); p2 != m2->atoms.end(); p2++){
            d = frame->dist((*p1)->pos_vect(), (*p2)->pos_vect());
            if (d < R_FACTOR*((*p1)->radius + (*p2)->radius)){
                return true;
            }
        }
    }
    return false;
}

int add_mol_neighbours(molecule * m1, molecule *m2){
    m1->add_neighbour(m2);
    m2->add_neighbour(m1);
    return 0;
}

int add_part_neighbours(particle *p1, particle *p2){
    p1->append(p2);
    p2->append(p1);
    return 0;
}

vector<int> short_neighbour_list(molecule * m, Frame * frame){
    vector<int> dist = vector<int>();
    int type;
    for (auto &neigh: m->my_neighbours){
        type = order_type(m, neigh, frame);
        if (type != 0) {
            dist.push_back(type);
        }
    }
    if (dist.size() == 0){
        dist.push_back(0);
    }
    return dist;
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
        for (int k = 0; k < (*m1).nump(); k++){
            d = frame->direction((*m1).atom_pos(k),com1);
            theta = atan2(d) + PI - rot + PI/2;
            // In format for gnuplot: theta, dist, circle radius, colour
            file << theta << "," << d.length() << "," << (*m1).atoms[k]->radius << "," << 1 << endl;
        }
        
        for (int j = 0; j < (*m1).num_neighbours(); j++){
            // Multiple interactions between particles
            if ((*m1).nint[j] > 1){
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
    if (fabs( Dbh - (Rb + Rh) ) < deltaD){
        if (fabs( Dhh - (Rh + Rh) ) < deltaD){
            if (dist( theta - PI, 2*PI )  < deltaT){
                // Anti Parallel 1
                return 3;
            }
            else if (fabs( Dbb - (Rb + Rb)) < deltaD){
                // Chiral
                return 5;
            }
        }
        else if (fabs( Dbb - (Rb + Rb)) < deltaD){
            if (dist (theta, 2*PI) < deltaT){
                // Parallel
                return 2;
            }
            else if (dist( theta - PI, 2*PI ) < deltaT){
                // Anti Parallel 2
                return 4;
            }
        }
    }
    if (my_mod(theta-PI/2, PI) < deltaT){
        // Perpendicular
        return 6;
    }
    return 0;
}
