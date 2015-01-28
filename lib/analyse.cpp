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
static float neighbour_size = 15;
static int n_angle_bins = 360;
static double angle_range = 2*PI;
static int com_colour = 3;
static int short_order_types = 7;

/*
 * Output files
 */

ofstream short_order_hist;
ofstream MSD_file;
ofstream rotations_file;
ofstream movie_file;

// Neighbour List
vector<vector<int>> neigh_list;

map<int, my_mean> collate_MSD, collate_c1, collate_c2, collate_c3, collate_c4;

int analyse(Frame *frame, vector<Frame *> key_frames, int print, int movie){
    
    /*
     * Variables
     */
    unsigned long num_key_frames = 1; //key_frames.size();
    angle_list frame_colour;
    
    // Dynamic Quantities
    vector<my_mean> MSD(num_key_frames), c1(num_key_frames), c2(num_key_frames), c3(num_key_frames), c4(num_key_frames);
    
    double phi;
    vector<molecule>::iterator mol, mol2;
    vect com1, com2, orient1, orient2;
    
    // Bond fraction
    my_mean neigh_frac;
    
    // Short range order
    vector<double> short_order_count(short_order_types,0);
    int short_order_colour = 0;
    int total_short_order = 0;
    
    // Local Order
    my_mean mean_local_order;
    
    // Global Order
    my_mean mean_global_order;
    
    // Circle Order
    my_mean mean_circle_order;
    
    // Angles
    vector<int> absolute_angles(n_angle_bins, 0);
    vector<int> relative_angles(n_angle_bins, 0);
    
    // Number of contacts
    vector<int> num_contact_dist;
    int total_contact = 0;
    
    // Number of neighbours
    vector<int> num_neigh_dist;
    int total_neigh = 0;
    
    // Pairing of Neighbours
    // TODO
    
    
    char fname[40];
    
    
    // Open new file if first;
    if (key_frames.size() == 0){
        
        // MSD file
        MSD_file.open("msd.csv");

        // Rotations file
        rotations_file.open("rotation.csv");
        rotations_file << "Steps, C1, C2, C3, C4" << endl;
        
        // Neighbour List
        neigh_list.resize(frame->num_mol(), vector<int>());
        
        // Movie File
        movie_file.open("trj/movie.lammpstrj");

    }
    
    
    ofstream short_order;
    ofstream order_parameter_file;
    ofstream gnuplot;
    ofstream angles;
    if (print){
        // Short range order
        short_order.open("short_order.csv");
        
        // Order Parameter
        snprintf(fname,40, "order/%010i.csv", frame->timestep);
        order_parameter_file.open(fname);
        
        // Gnuplot
        snprintf(fname,40, "trj_contact/%010i.dat", frame->timestep);
        gnuplot.open(fname);
        gnuplot << frame->get_a() << " " << frame->get_height() << endl << endl;
        
        // Angles
        snprintf(fname,40, "trj_contact/%010i-angles.csv", frame->timestep);
        angles.open(fname);
    }
    if (movie){
        // Lammpstrj frame info
        movie_file << "ITEM: TIMESTEP" << endl << frame->timestep << endl;
        movie_file << "ITEM: NUMBER OF ATOMS" << endl << frame->num_atoms() << endl;
        movie_file << "ITEM: BOX BOUNDS xy xz yz pp pp pp" << endl;
        movie_file << 0 << " " << frame->get_a() << " " << frame->get_tilt() << endl;
        movie_file << 0 << " " << frame->get_height() << " " << 0 << endl;
        movie_file << "-0.5 0.5 0" << endl;
        movie_file << "ITEM: ATOMS id mol type x y z vx vy vz" << endl;
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
        // Absolute Angle
        int angle_bin = (int) (angle(&(*mol), frame) / (angle_range/(n_angle_bins)));
        absolute_angles.at(angle_bin)++;

        my_mean my_local_order = my_mean();
        if (print){
            vect d;
            for (auto & x: (*mol).atoms){
                d = frame->direction(x->pos_vect(), com1);
                double phi = atan2(d) + 5*PI/2 - atan2(orient1);
                
                // Print mol for short range order
                short_order << phi << "," << d.length() << "," << x->radius << "," << 1 << endl;
                
                // Print mol for order
                order_parameter_file << phi << "," << d.length() << "," << x->type << endl;
            }
            d = vect();
            phi = 0;
            order_parameter_file << phi << "," << d.length() << "," << com_colour << endl;
        }
        
        /*
         * Find Neighbours
         */
        if (key_frames.size() == 0){
            for (mol2 = mol+1; mol2 != frame->molecules.end(); mol2++){
                vector<particle *>::iterator p1, p2;
                /*
                 * Properties of mol2
                 */
                com2 = (*mol2).COM();
                orient2 = orientation(&(*mol2), frame);
                // Global Order
                mean_global_order.add(pow(dot_product(orient1, orient2),2));
                // Relative Angles
                angle_bin = ((int) (angle(&(*mol), frame) - angle(&(*mol2), frame) / (angle_range/(n_angle_bins))) + n_angle_bins) % n_angle_bins;
                relative_angles.at(angle_bin)++;

                double d_com = frame->dist(com1, com2);
                double d_atom;
                if (d_com < neighbour_size){
                    neigh_list.at((*mol).id-1).push_back((*mol2).id);
                    if (print){
                        // print order atomic
                        for (auto & x: (*mol2).atoms){
                            vect d = frame->direction(x->pos_vect(), com1);
                            double phi = atan2(d) + 5*PI/2 - atan2(orient1);
                            order_parameter_file << phi << "," << d.length() << "," << x->type << endl;
                        }
                        // Print order COM
                        vect d = frame->direction(com2, com1);
                        double phi = atan2(d) + 5*PI/2 - atan2(orient1);
                        order_parameter_file << phi << "," << d.length() << "," << com_colour << endl;
                    }
                    for (auto &p2: (*mol2).atoms){
                        for (auto &p1: (*mol).atoms){
                            d_atom = frame->dist(p1->pos_vect(), p2->pos_vect());
                            d_atom /= p1->radius + p2->radius;
                            if (d_atom < R_FACTOR){
                                // Average Neighbour Bond fraction
                                neigh_frac.add(d_atom);
                            
                                // Find Neighbours
                                add_mol_neighbours(&(*mol), &(*mol2));
                                add_part_neighbours(p1,p2);
                            
                                // Local Order
                                mean_local_order.add(pow(dot_product(orient1, orient2),2));
                                my_local_order.add(pow(dot_product(orient1, orient2),2));
                            }
                        }
                    }
                }
            }
        }
        else {
            vector<int> neighbours = neigh_list.at((*mol).id-1);
            molecule * moln;
            for (auto &m: neighbours){
                moln = &frame->molecules.at(m-1);
                vector<particle *>::iterator p1, p2;
                /*
                 * Properties of mol2
                 */
                com2 = moln->COM();
                orient2 = orientation(moln, frame);
                // Global Order
                mean_global_order.add(pow(dot_product(orient1, orient2),2));
                // Relative Angles
                angle_bin = ((int) (angle(&(*mol), frame) - angle(moln, frame) / (angle_range/(n_angle_bins))) + n_angle_bins) % n_angle_bins;
                relative_angles.at(angle_bin)++;

                double d_atom;
                
                if (print){
                    // print order atomic
                    for (auto & x: moln->atoms){
                        vect d = frame->direction(x->pos_vect(), com1);
                        double phi = atan2(d) + 5*PI/2 - atan2(orient1);
                        order_parameter_file << phi << "," << d.length() << "," << x->type << endl;
                    }
                    // Print order COM
                    vect d = frame->direction(com2, com1);
                    double phi = atan2(d) + 5*PI/2 - atan2(orient1);
                    order_parameter_file << phi << "," << d.length() << "," << com_colour << endl;
                }
                
                for (p2 = moln->atoms.begin(); p2 != moln->atoms.end(); p2++){
                    for (p1 = (*mol).atoms.begin(); p1 != (*mol).atoms.end(); p1++){
                        d_atom = frame->dist((*p1)->pos_vect(), (*p2)->pos_vect());
                        d_atom /= (*p1)->radius + (*p2)->radius;
                        if (d_atom < R_FACTOR){
                            // Average Neighbour Bond fraction
                            neigh_frac.add(d_atom);
                            
                            // Find Neighbours
                            add_mol_neighbours(&(*mol), moln);
                            add_part_neighbours(*p1,*p2);
                            
                            // Local Order
                            mean_local_order.add(pow(dot_product(orient1, orient2),2));
                            my_local_order.add(pow(dot_product(orient1, orient2),2));
                        }
                    }
                }
                
            }
        }
        /*
         * Every frame
         */
        
        // Short Order
        vector<int> this_short_order;
        for (int j = 0; j < (*mol).num_neighbours(); j++){
            molecule *m2 = (*mol).my_neighbours[j];
            short_order_colour = order_type(&(*mol), m2, frame);
            if (short_order_colour){
                this_short_order.push_back(short_order_colour);
                if (print){
                    vect d;
                    double theta;
                    for (int k = 0; k < m2->nump(); k++){
                        d = frame->direction(m2->atom_pos(k),com1);
                        theta = atan2(d) + 5*PI/2 - atan2(orient1);
                        short_order << theta << "," << d.length() << "," << 0.04 << "," << short_order_colour << endl;
                    }
                }
            }
                    }
        if (this_short_order.size() == 0){
            short_order_count.at(0) += 1;
            
        }
        else {
            for (auto &i: this_short_order){
                short_order_count.at(i) += 1.0/this_short_order.size();
            }
        }
        total_short_order++;
        
        int frame_no = 0;
        if (key_frames.size() > 0){
            Frame * key = key_frames.front();
        //for (auto key: key_frames){
            
            // MSD -> from each key frame
            com2 = key->molecules.at((*mol).id-1).COM();
            MSD.at(frame_no).add(frame->direction(com1,com2).length());
            
            // Rotation -> from each key frame
            orient2 = orientation(&(key->molecules.at((*mol).id-1)), key);
            phi = dot_product(orient1, orient2);
            c1.at(frame_no).add(legendre(1,phi));
            c2.at(frame_no).add(legendre(2,phi));
            c3.at(frame_no).add(legendre(3,phi));
            c4.at(frame_no).add(legendre(4,phi));
            
            frame_no++;
        }
        
        // Composition
        
        // Crystal Front
        
        // Circle Order
        // TODO
        
        /*
         * Some frames
         */
        if (print){
            int n;
            // Print angles
            // TODO
            
            // Hist num contact
            n = (*mol).num_contacts();
            if (n+1 > num_contact_dist.size()){
                num_contact_dist.resize(n+1, 0);
            }
            ++num_contact_dist.at(n);
            total_contact += n;

            
            // Hist num neighbours
            n = (*mol).num_neighbours();
            if (n+1 > num_neigh_dist.size()){
                num_neigh_dist.resize(n+1, 0);
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
            // id mol type x y z vx vy vz
            vect cart_com, d;
            double mol_colour = 0;
            
            cart_com = frame->cartesian(mol->COM());
            cart_com = wrap_x(cart_com, frame->get_a());
            
            if ((*mol).num_contacts() > 8){
                mol_colour = 1;
            }
            std::vector<particle *>::iterator i;
            for (i = mol->atoms.begin(); i != mol->atoms.end(); i++){
                d = frame->cartesian(direction(mol->COM(), (*i)->pos_vect()));
                movie_file << (*i)->id << " " << mol->id << " " << (*i)->type << " " << cart_com + d << " 0 " << mol_colour << " 0 0" << endl;
            }
        }
    }
    /*
     * Output
     */
    
    
    // Collate key frames
    for (int k = 0; k < key_frames.size(); k++){
        int steps = frame->timestep - key_frames.at(k)->timestep;
        collate_MSD[steps].add(MSD.at(k).get_mean());
        collate_c1[steps].add(c1.at(k).get_mean());
        collate_c2[steps].add(c2.at(k).get_mean());
        collate_c3[steps].add(c3.at(k).get_mean());
        collate_c4[steps].add(c4.at(k).get_mean());
    }

    // Percentage in each short range ordering
    if (key_frames.size() == 0){
        short_order_hist.open("short_order_hist.csv");
        short_order_hist << "Timestep,No Colour, None, Parallel, Anti Parallel 1, Anti Parallel 2, Chiral, Perpendicular" << endl;
    }
    else{
        short_order_hist.open("short_order_hist.csv", ofstream::app);
    }
    short_order_hist << frame->get_timestep();
    for (auto &i: short_order_count){
        short_order_hist << "," << double(i)/total_short_order;
    }
    short_order_hist << endl;
    short_order_hist.close();
    
    if (key_frames.size() > 0){
        // MSD plot
        MSD_file  << frame->timestep << "," << MSD.front().get_mean() << endl;
    
        // Rotations Plot
        rotations_file << frame->timestep << "," << c1.front().get_mean() << "," << c2.front().get_mean() << \
        "," << c3.front().get_mean() << "," << c4.front().get_mean() << endl;
    }
    
    // Only want to print this the last time
    if (print && key_frames.size() > 0){
        
        // Relaxation times
        int t1 = 0, t2 = 0, t3 = 0, t4 = 0;
        for (auto &x: collate_c1){
            if (x.second.get_mean() < 0){
                t1 = x.first;
                break;
            }
        }
        for (auto &x: collate_c2){
            if (x.second.get_mean() < 0.25){
                t2 = x.first;
                break;
            }
        }
        for (auto &x: collate_c3){
            if (x.second.get_mean() < 0){
                t3 = x.first;
                break;
            }
        }
        for (auto &x: collate_c4){
            if (x.second.get_mean() < 0.1406){
                t4 = x.first;
                break;
            }
        }
        cout << "T1: " << t1 << endl;
        cout << "T2: " << t2 << endl;
        cout << "T3: " << t3 << endl;
        cout << "T4: " << t4 << endl;
        
        
        // Print Diffusion Constant
        my_mean diffusion;
        for (auto &x: collate_MSD){
            diffusion.add(x.second.get_mean()/(4*x.first*DELTA_T));
        }
        cout << "Diffusion Constant: " << diffusion.get_mean() << endl;
        
        // Print Average bond fraction
        cout << "Average Bond Fraction: " << neigh_frac.get_mean() << endl;
        
        // Print Local Order
        cout << "Local Order: " << mean_local_order.get_mean() << endl;
        
        // Print Global Order
        cout << "Global Order: " << mean_global_order.get_mean() << endl;
        
        // Print Circle Order
        // TODO
        
        // Print num contact
        ofstream hist_contact;
        snprintf(fname, 40, "stats/contacts.dat");
        hist_contact.open(fname);
        cout << "Num Contacts: " << double(total_contact)/frame->num_mol() << endl;
        for (int i = 0; i < num_contact_dist.size(); i++){
            hist_contact << i << " " << num_contact_dist.at(i) << endl;
        }
        hist_contact.close();

        // Print num neighbours
        ofstream hist_neigh;
        snprintf(fname, 40, "stats/neighbours.dat");
        hist_neigh.open(fname);
        cout << "Num Neighbours: " << double(total_neigh)/frame->num_mol() << endl;
        for (int i = 0; i < num_neigh_dist.size(); i++){
            hist_neigh << i << " " << num_neigh_dist.at(i) << endl;
        }
        hist_neigh.close();
        
        // Print Pairing
        // TODO
    }
    if (print){
        
        // Print angle
        for (int j = 0; j < n_angle_bins; j++){
            angles << angle_range/(n_angle_bins) * (j) * 180/PI << ","\
            << absolute_angles.at(j)/(float) frame->num_mol() << ","\
            << relative_angles.at(j)/(((float) frame->num_mol())/2 * frame->num_mol()) << endl;
        }
        // Close Files
        short_order.close();
        gnuplot.close();
        angles.close();
        order_parameter_file.close();
    }
    return 0;
}