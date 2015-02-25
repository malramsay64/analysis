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

map<int, my_mean> collate_MSD, collate_MFD, collate_c1, collate_c2, collate_struct;

ofstream MSD_file, rotations_file, movie_file, short_order_file, struct_file;

static double radial_plot = 15;
static double radial_cutoff = 20;
static int short_order_types = 7;


int mod_analyse(Frame * frame, std::vector<Frame *> key_frames, int print, int movie, int dist){
    
    if (key_frames.size() == 0){
        mod_neigh_list = vector<vector<int>>(frame->num_mol(), vector<int>());
        movie_file.open("trj/movie.lammpstrj");
        print_movie(&movie_file, frame);
        
        // Dynamics
        MSD_file.open("msd.csv");
        rotations_file.open("rotation.csv");
        struct_file.open("struct.csv");
        
        // Short Order
        short_order_file.open("short_order_hist.csv");
        short_order_file << "Timestep,No Colour,None,Parallel,Anti Parallel 1,Anti Parallel 2,Chiral,Perpendicular" << endl;

    }
    
    distribution<int> num_neigh = distribution<int>(MAX_MOL_CONTACTS);
    distribution<int> num_contact = distribution<int>(MAX_MOL_CONTACTS);
    distribution<int> pairing = distribution<int>(MAX_MOL_CONTACTS);
    distribution<int> radial = distribution<int>(1000, radial_plot);
    
    distribution<my_mean> pair_contact = distribution<my_mean>(MAX_MOL_CONTACTS);
    distribution<my_mean> pair_neigh = distribution<my_mean>(MAX_MOL_CONTACTS);
    
    my_mean neigh_frac;
    //my_mean struct_func;
    
    
    int num_key_frames = (int) key_frames.size();
    vector<my_mean> MSD(num_key_frames), MFD(num_key_frames), c1(num_key_frames),\
    c2(num_key_frames), struct_func(num_key_frames);
    
    // Short order histogram
    distribution<double> short_order = distribution<double>(short_order_types);
    
    
    ofstream gnuplot;
    ofstream short_range;
    
    /*
     * Find Neighbours
     */
    bool recompute = false;
    vector<molecule *> recompute_list;
    for (auto &mol: frame->molecules){
        recompute = find_mol_neighbours(&mol, frame, &mod_neigh_list);
        if (recompute) {
            break;
        }
    }
    if (recompute){
        // Delete any neighbours
        mod_neigh_list = vector<vector<int>>(frame->num_mol(), vector<int>());
        for (auto &mol: frame->molecules){
            mol.delete_neighbours();
        }
        for (auto &mol: frame->molecules){
            find_mol_neighbours(&mol, frame, &mod_neigh_list);
        }
    }

    /*
     * Analysis
     */
    for (auto &mol: frame->molecules){
        num_neigh.add(mol.num_neighbours());
        num_contact.add(mol.num_contacts());
        pairing.add(mol.pairing());
        
        pair_contact.add(mol.max_pairing(), mol.num_contacts());
        pair_neigh.add(mol.max_pairing(), mol.num_neighbours());
        
        // Short range order
        short_order.add(short_neighbour_list(&mol, frame));
        
        
        // MSD / Rotations
        double phi;
        int k = 0;
        double displacement;
        for (auto key: key_frames) {
            molecule * mol2 = &key->molecules.at(mol.index());
            displacement = frame->dist(mol.COM(), mol2->COM());
            // Displacement
            MSD.at(k).add(pow(displacement,2));
            MFD.at(k).add(pow(displacement,4));
            // Rotations
            phi = dot_product(orientation(&mol, frame), orientation(mol2, frame));
            c1.at(k).add(legendre(1,phi));
            c2.at(k).add(legendre(2,phi));
            // Structural Relaxation
            //cout << k << " " <<  struct_relax(&mol, key) << endl;
            struct_func.at(k).add(struct_relax(&mol, key));
            
            k++;
        }
        
        dyn_queue r = dyn_queue(&mol);
        molecule *mol2 = r.pop();
        double distance = 0;
        
        // All molecules within cutoff
        while ( mol2 && distance < radial_cutoff ){
            distance = frame->dist(mol.COM(), mol2->COM());
            radial.add(distance);

            mol2 = r.pop();
        }
        if (movie){
            print_movie(&movie_file, frame, &mol);
        }
        if (print){
            print_short_order(&short_range, &mol, frame);
        }
    }
    
    for (int k = 0; k < key_frames.size(); k++){
        int steps = frame->timestep - key_frames.at(k)->timestep;
        collate_MSD[steps].add(MSD.at(k).get_mean());
        collate_MFD[steps].add(MFD.at(k).get_mean());
        collate_c1[steps].add(c1.at(k).get_mean());
        collate_c2[steps].add(c2.at(k).get_mean());
        collate_struct[steps].add(struct_func.at(k).get_mean());
        //cout << steps << " " << struct_func.at(k).get_mean() << endl;
    }

    print_time_distribution(&short_order, frame->timestep, &short_order_file);
    
    
    if (print && key_frames.size() > 0){
        // Print num neighbours
        //cout << "Num Neighbours: " << num_neigh.get_mean() << endl;
        print_distribution<int>(&num_neigh, "stats/neighbours.dat");
        
        // Print num contacts
        //cout << "Num Contacts: " << num_contact.get_mean() << endl;
        print_distribution<int>(&num_contact, "stats/contacts.dat");
        
        // Print num contacts
        //cout << "Pairing: " << pairing.get_mean() << endl;
        print_distribution<int>(&pairing, "stats/pairing.dat");
        
        // Print radial distribution
        print_radial_distribution(&radial, "radial_dist.dat", frame->num_mol(), frame->get_area());
        
        // Print short range order
        //print_distribution<double>(&short_order, "short_order_dist.dat");
        
        // Displacement
        double alpha;
        MSD_file << "Timestep,MSD,MFD,\\alpha" << endl;
        for (auto d: collate_MSD){
            alpha = collate_MFD.at(d.first).get_mean() / (2*pow(d.second.get_mean(),2)) - 1;
            MSD_file << d.first << "," << d.second.get_mean() << "," << collate_MFD.at(d.first).get_mean() << "," << alpha << endl;
        }
        
        // Rotations
        int t1 = 0, t2 = 0;
        rotations_file << "Timestep,P1(t),P2(t)" << endl;
        for (auto c: collate_c1){
            rotations_file << c.first << "," << c.second.get_mean() << "," << \
            collate_c2.at(c.first).get_mean() << "," << endl;
            
            if (c.second.get_mean() < 1/CONST_E && t1 == 0){
                t1 = c.first;
            }
            else if (collate_c2.at(c.first).get_mean() < 1/CONST_E && t2 == 0){
                t2 = c.first;
            }
        }
        print_relax_time("t1: ", t1);
        print_relax_time("t2: ", t2);
        
        // Structure
        int ts = 0;
        struct_file << "Timestep,F(t),\\xi(t)"  << endl;
        for (auto c: collate_struct){
            struct_file << c.first << "," << c.second.get_mean() \
            << "," << c.second.get_variance() << endl;
            if (c.second.get_mean() < 1/CONST_E && ts == 0){
                ts = c.first;
            }
        }
        print_relax_time("Struct-relax: ", ts);
        cout << "Diffusion-constant: " << \
        collate_MSD.at(frame->timestep).get_mean()/(4*frame->timestep*STEP_SIZE)  \
        << endl;
        
        print_distribution(&pair_contact, "stats/pair_contact.dat");
        print_distribution(&pair_neigh, "stats/pair_neigh.dat");
        
        print_moved(key_frames.front(), frame);
    }
    if  (print){
        print_frame(frame);
        
    }
    return 0;
}
