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
    
    distribution<int> num_neigh = distribution<int>(MAX_MOL_CONTACTS);
    distribution<int> num_contact = distribution<int>(MAX_MOL_CONTACTS);
    distribution<int> pairing = distribution<int>(MAX_MOL_CONTACTS);
    distribution<int> radial = distribution<int>(1000, radial_plot);
    
    distribution<double> short_order = distribution<double>(short_order_types);
    
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
        
        // Short range order
        short_order.add(short_neighbour_list(&mol, frame));
        
        dyn_queue r = dyn_queue(&mol);
        molecule *mol2 = r.pop();
        double distance = 0;
        
        // All molecules within cutoff
        while ( mol2 && distance < radial_cutoff ){
            distance = frame->dist(mol.COM(), mol2->COM());
            radial.add(distance);

            
            mol2 = r.pop();
        }

        
    }
    if (print){
        // Print num neighbours
        cout << "Num Neighbours: " << num_neigh.get_mean() << endl;
        print_distribution<int>(&num_neigh, "stats/neighbours.dat");
        
        // Print num contacts
        cout << "Num Contacts: " << num_contact.get_mean() << endl;
        print_distribution<int>(&num_contact, "stats/contacts.dat");
        
        // Print num contacts
        cout << "Pairing: " << pairing.get_mean() << endl;
        print_distribution<int>(&pairing, "stats/pairing.dat");
        
        // Print radial distribution
        print_radial_distribution(&radial, "radial_dist.dat", frame->num_mol(), frame->get_area());
        
        // Print short range order
        print_distribution<double>(&short_order, "short_order_dist.dat");
        
    }
    return 0;
}
