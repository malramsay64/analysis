//
//  analyse.cpp
//  analysis
//
//  Created by Malcolm Ramsay on 14/01/2015.
//  Copyright (c) 2015 Malcolm Ramsay. All rights reserved.
//

#include "analyse.h"

using namespace std;

static float R_FACTOR = pow(2,1./6);

int analyse(Frame *frame, vector<Frame *> key_frames, int print=0, int movie=0){
    /* 
     * Variables
     */
    unsigned long num_key_frames = key_frames.size();
    angle_list frame_colour;
    
    // Dynamic Quantities
    vector<my_mean> MSD(num_key_frames), c1(num_key_frames), c2(num_key_frames), c3(num_key_frames), c4(num_key_frames);
    
    double phi;
    vector<molecule>::iterator mol, mol2;
    vect com1, com2, orient1, orient2;
    
    // Bond fraction
    my_mean neigh_frac;
    
    // Short range order
    vector<int> short_order_count(6,0);
    int short_order_colour;
    int total_short_order;
    
    // Number of contacts
    vector<int> num_contact_dist;
    int total_contact = 0;
    
    // Number of neighbours
    vector<int> num_neigh_dist;
    int total_neigh = 0;
    
    // Pairing of Neighbours
    // TODO
    
    /*
     * Output files
     */
    ofstream gnuplot;
    ofstream short_order;
    ofstream short_order_hist;
    ofstream hist_neigh;
    ofstream hist_contact;
    char fname[30];
    
    // Short range order histogram - append to existing
    short_order_hist.open("short_order_hist.csv", ofstream::app);
    short_order_hist << frame->get_timestep() << endl;
    
    if (print){
        // Short range order
        short_order.open("short_order.csv");
        //Print central molecule
        // TODO
        
        // Gnuplot
        snprintf(fname,30, "trj_contact/%10i.dat", frame->timestep);
        gnuplot.open(fname);
        gnuplot << frame->get_a() << " " << frame->get_height() << endl << endl;
        
        // Neighbours Histogram
        snprintf(fname, 30, "stats/num_neighbours.dat");
        hist_neigh.open(fname);
        
        // Contacts Histogram
        snprintf(fname, 30, "stats/num_contacts.dat");
        hist_contact.open(fname);
    
    }
    
    /*
     * Parallelise TODO
     */
    for (mol = frame->molecules.begin(); mol != frame->molecules.end(); mol++){
        /*
         * Properties of mol
         */
        // Center of mass
        com1 = (*mol).COM();
        // Orientation of molecule
        orient1 = orientation(&(*mol), frame);
        // Colour of molecule
        (*mol).set_colour(frame_colour.push(atan2(orient1)+PI));
        
        /*
         * Find Neighbours
         * Regions TODO
         */
        for (mol2 = mol+1; mol2 != frame->molecules.end(); mol2++){
            vector<particle *>::iterator p1, p2;
            double d = frame->dist(com1, com2);
            if (d < 4){
                for (p1 = (*mol).atoms.begin(); p1 != (*mol).atoms.end(); p1++){
                    for (p2 = (*mol2).atoms.begin(); p2 != (*mol2).atoms.end(); p2++){
                        d = frame->dist((*p1)->pos_vect(), (*p2)->pos_vect());
                        d /= (*p1)->radius + (*p2)->radius;
                        if (d < R_FACTOR){
                            // Average Neighbour Bond fraction
                            neigh_frac.add(d);
                            // Find Neighbours
                            add_mol_neighbours(&*mol, &*mol2);
                            // Short range order
                            short_order_colour = order_type(&(*mol), &(*mol2), frame);
                            short_order_count.at(short_order_colour)++;
                            total_short_order++;
                            if (print){
                                // Print short order
                                // TODO
                                
                                // Print order
                                // TODO
                            }
                        }
                    }
                }
            }
            else if (d < 15){
                if (print){
                    // print order
                    // TODO
                }
            }
        }

        /*
         * Every frame
         */
        
        vector<Frame *>::iterator key;
        int frame_no = 0;
        for (key = key_frames.begin(); key != key_frames.end(); key++){
            
            // MSD -> from each key frame
            com2 = (*key)->molecules.at((*mol).id).COM();
            MSD.at(frame_no).add(frame->direction(com1,com2).length());
            
            // Rotation -> from each key frame
            orient2 = orientation(&((*key)->molecules.at((*mol).id)), *key);
            phi = dot_product(orient1, orient2);
            c1.at(frame_no).add(legendre(1,phi));
            c2.at(frame_no).add(legendre(2,phi));
            c3.at(frame_no).add(legendre(3,phi));
            c4.at(frame_no).add(legendre(4,phi));
            
            
            frame_no++;
        }
        
        // Composition
        
        // Crystal Front
        
    
        /*
         * Some frames
         */
        if (print){
            int n;
            // Print angles
            
            // Hist num contact
            n = (*mol).num_contacts();
            if (n > num_contact_dist.size()){
                num_contact_dist.resize(n, 0);
            }
            ++num_contact_dist.at(n);
            total_contact += n;

            
            // Hist num neighbours
            n = (*mol).num_neighbours();
            if (n > num_neigh_dist.size()){
                num_neigh_dist.resize(n, 0);
            }
            ++num_neigh_dist.at(n);
            total_neigh += n;
            
            // Hist Pairing
            // TODO
            
            // Print frame
            print_mol(gnuplot, &*mol, frame);
            gnuplot << endl;
        }
        /*
         * Video output
         */
        if (movie){
            // Print lammpstrj frame (every frame)
        }
    }
    /*
     * Output
     */
    
    vector<int>::iterator i;
    // Collate key frames
    
    // Relaxation times
    
    // Print data
    
    // Percentage in each short range ordering
    int count = 0;
    for (i = short_order_count.begin(); i != short_order_count.end(); i++){
        short_order_hist << count << "," << *i << endl;
    }
    short_order_hist << endl;
    short_order_hist.close();
    
    
    if (print){
        vector<int>::iterator i;
        
        // Print num contact
        cout << "Num Contacts: " << double(total_contact)/frame->num_mol() << endl;
        for (i = num_contact_dist.begin(); i != num_contact_dist.end(); i++){
            hist_contact << i - num_contact_dist.begin() << " " << *i << endl;
        }
        hist_contact.close();

        
        // Print num neighbours
        cout << "Num Neighbours: " << double(total_neigh)/frame->num_mol() << endl;
        for (i = num_neigh_dist.begin(); i != num_neigh_dist.end(); i++){
            hist_neigh << i - num_neigh_dist.begin() << " " << *i << endl;
        }
        hist_neigh.close();
        
        // Print Pairing
        // TODO
    }
    return 0;
}