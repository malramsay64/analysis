//
//  modular.cpp
//  analysis
//
//  Created by Malcolm Ramsay on 6/02/2015.
//  Copyright (c) 2015 Malcolm Ramsay. All rights reserved.
//

#include "modular.h"

using namespace std;

vector<vector<int>> mod_neigh_list;

static double radial_plot = 15;
static double radial_cutoff = 20;
static int short_order_types = 7;


int mod_analyse(Frame * frame, std::vector<Frame *> key_frames, int print, int movie, int dist){
    
    if (key_frames.size() == 0){
        mod_neigh_list = vector<vector<int>>(frame->num_mol(), vector<int>());
    }
    
    distribution num_neigh = distribution(MAX_MOL_CONTACTS);
    distribution num_contact = distribution(MAX_MOL_CONTACTS);
    distribution pairing = distribution(MAX_MOL_CONTACTS);
    distribution radial = distribution(1000, radial_plot);
    
    /*
     * Find Neighbours
     */
    for (auto &mol: frame->molecules){
        find_mol_neighbours(&mol, frame, &mod_neigh_list);
    }
    
    /*
     * Analysis
     */
    for (auto &mol: frame->molecules){
        
        num_neigh.add(mol.num_neighbours());
        num_contact.add(mol.num_contacts());
        pairing.add(mol.pairing());
        
        dyn_queue r = dyn_queue(&mol);
        molecule *mol2 = r.pop();
        double distance = 0;
        while ( mol2 && distance < radial_cutoff ){
            distance = frame->dist(mol.COM(), mol2->COM());
            radial.add(distance);
            mol2 = r.pop();
        }

        
    }
    if (print){
        // Print num neighbours
        cout << "Num Neighbours: " << num_neigh.get_mean() << endl;
        print_distribution(&num_neigh, "stats/neighbours.dat");
        
        // Print num contacts
        cout << "Num Contacts: " << num_contact.get_mean() << endl;
        print_distribution(&num_contact, "stats/contacts.dat");
        
        // Print num contacts
        cout << "Pairing: " << pairing.get_mean() << endl;
        print_distribution(&pairing, "stats/pairing.dat");
        
        // Print radial distribution
        print_radial_distribution(&radial, "radial_dist.dat", frame->num_mol(), frame->get_area());
        
    }
    return 0;
}
