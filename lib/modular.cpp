//
//  modular.cpp
//  analysis
//
//  Created by Malcolm Ramsay on 6/02/2015.
//  Copyright (c) 2015 Malcolm Ramsay. All rights reserved.
//

#include "modular.h"

using namespace std;
// Neighbour List
vector<vector<int>> mod_neigh_list;

map<int, my_mean> collate_MSD, collate_MFD, collate_c1, collate_c2, collate_struct;
vector<map<int,my_mean>> collate_regio_c1, collate_regio_c2, collate_regio_MSD;

vector<my_mean> regio_orientation, regio_circle, velocity;

ofstream MSD_file, rotations_file, movie_file, short_order_file, struct_file, regio_file, order_file;



static double radial_plot = 15;
static double radial_cutoff = 20;
static int short_order_types = 7;
static int regio_res = 50;
static int radial_res = 360;
static int theta_res = 360;
double dtheta = 2*PI/theta_res;

double rot_cut = 0.5;
double dif_cut = 0.3;

int mod_analyse(Frame * frame, std::vector<Frame *> key_frames, int print, int dist){
    
    // Initialise Variables
    int num_key_frames = (int) key_frames.size();
    
    vector<my_mean> MSD(num_key_frames), MFD(num_key_frames), c1(num_key_frames),\
    c2(num_key_frames), struct_func(num_key_frames);
    
    vector<vector<my_mean>> regio_c1, regio_c2, regio_MSD;
    //vector<my_mean> regio_orientation, regio_circle;
    distribution<int> num_neigh, num_contact, pairing, radial, radial_part, num_neigh_mob[4], num_contact_mob[4];
    distribution<my_mean> pair_contact, pair_neigh;
    distribution<double> short_order;
    vector<distribution<int>> radial2d_rel, radial2d_abs, radial2d_large, radial2d_small, radial2d_part;
    
    my_mean neigh_frac;
    my_mean hexatic_order;
    my_mean circle_order;
    my_mean orientational_order;
    
    ofstream gnuplot;
    ofstream short_range;
    
    if (key_frames.size() == 0){
        mod_neigh_list = vector<vector<int>>(frame->num_mol(), vector<int>());
        
        // Dynamics
        MSD_file.open("msd.csv");
        rotations_file.open("rotation.csv");
        struct_file.open("struct.csv");
        
        if (time_structure){
            // Short Order
            short_order_file.open("short_order_hist.csv");
            short_order_file << "Timestep,No Colour,None,Parallel,Anti Parallel 1,Anti Parallel 2,Chiral,Perpendicular" << endl;
        
            // Hexatic Ordering
            order_file.open("order.csv");
            order_file << "Timestep,Hexatic,Circle,Orientation" << endl;
        }
        
        if (regio){
            // Regio relaxations
            collate_regio_c1 = vector<map<int,my_mean>>(regio_res);
            collate_regio_c2 = vector<map<int,my_mean>>(regio_res);
            collate_regio_MSD = vector<map<int,my_mean>>(regio_res);
            regio_orientation = vector<my_mean>(regio_res);
            regio_circle = vector<my_mean>(regio_res);
        }
        
        if (movie){
            movie_file.open("trj/movie.lammpstrj");
            print_movie(&movie_file, frame);
        }
        
        if (m_orient){
            velocity = vector<my_mean>(radial_res);
        }
    }
    
    if (time_structure || print){
        num_neigh = distribution<int>(MAX_MOL_CONTACTS);
        num_contact = distribution<int>(MAX_MOL_CONTACTS);
        pairing = distribution<int>(MAX_MOL_CONTACTS);
        for (int i = 0; i < 4; i ++){
            num_contact_mob[i] = distribution<int>(MAX_MOL_CONTACTS);
            num_neigh_mob[i] = distribution<int>(MAX_MOL_CONTACTS);
        }
        radial = distribution<int>(radial_res, radial_plot);
        radial_part = distribution<int>(radial_res, radial_plot);
        radial2d_rel = vector<distribution<int>>(theta_res, distribution<int>(theta_res, radial_plot));
        radial2d_abs = vector<distribution<int>>(theta_res, distribution<int>(theta_res, radial_plot));
        radial2d_large = vector<distribution<int>>(theta_res, distribution<int>(theta_res, radial_plot));
        radial2d_small = vector<distribution<int>>(theta_res, distribution<int>(theta_res, radial_plot));
        radial2d_part = vector<distribution<int>>(theta_res, distribution<int>(theta_res, radial_plot));
    
        pair_contact = distribution<my_mean>(MAX_MOL_CONTACTS);
        pair_neigh = distribution<my_mean>(MAX_MOL_CONTACTS);
        
    }
    
    if (regio){
        regio_c1 = vector<vector<my_mean>>(regio_res, vector<my_mean>(num_key_frames));
        regio_c2 = vector<vector<my_mean>>(regio_res, vector<my_mean>(num_key_frames));
        regio_MSD = vector<vector<my_mean>>(regio_res, vector<my_mean>(num_key_frames));
    }
    
    if (time_structure || print){
        short_order = distribution<double>(short_order_types);
        
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
    }

    
    /*
     * Analysis
     */
    for (auto &mol: frame->molecules){
        if (time_structure || print){
            num_neigh.add(mol.num_neighbours());
            num_contact.add(mol.num_contacts());
            pairing.add(mol.pairing());
            
            
            pair_contact.add(mol.max_pairing(), mol.num_contacts());
            pair_neigh.add(mol.max_pairing(), mol.num_neighbours());
        
            // Short range order
            short_order.add(short_neighbour_list(&mol, frame));
        
            // Order
            hexatic_order.add(fabs(hexatic(6, &mol, frame)));
            circle_order.add(circle_ordering(&mol));
            orientational_order.add(orient_ordering(&mol));
            
            // Regio
            if (regio){
                int index;
                // x divide
                if (frame->get_a() > frame->get_height()){
                    index = pos_def_mod(int (mol.COM().x/(2*PI/regio_res)),regio_res);
                }
                // y divide
                else {
                    index = pos_def_mod(int (mol.COM().y/regio_res),regio_res);
                }
                regio_circle.at(index).add(circle_ordering(&mol));
                regio_orientation.at(index).add(orient_ordering(&mol));
            }
        }
            
        // MSD / Rotations
        double phi;
        int k = 0;
        double displacement;
        for (auto key: key_frames) {
            molecule * mol2 = &key->molecules.at(mol.index());
            displacement = frame->cartesian(mol2->moved_COM()-mol.moved_COM()).length();
            // Displacement
            MSD.at(k).add(pow(displacement,2));
            MFD.at(k).add(pow(displacement,4));
            // Rotations
            phi = dot_product(orientation(&mol, frame), orientation(mol2, frame));
            c1.at(k).add(legendre(1,phi));
            c2.at(k).add(legendre(2,phi));
            // Structural Relaxation
            struct_func.at(k).add(struct_relax(&mol, key));
            
            
            if (time_structure || print){
                // Contact distributions
                if (displacement > dif_cut){
                    num_neigh_mob[1].add(mol.num_neighbours());
                    num_contact_mob[1].add(mol.num_contacts());
                }
                else {
                    num_neigh_mob[0].add(mol.num_neighbours());
                    num_contact_mob[0].add(mol.num_contacts());
                    
                }
                
                if (fabs(angle(mol2, frame) - angle(&mol, frame)) > rot_cut){
                    num_neigh_mob[2].add(mol.num_neighbours());
                    num_contact_mob[2].add(mol.num_contacts());
                }
                else {
                    num_neigh_mob[3].add(mol.num_neighbours());
                    num_contact_mob[3].add(mol.num_contacts());
                }
            }
            
            // Regio
            if (regio){
                int index;
                // x divide
                if (key->get_a() > key->get_height()){
                    index = pos_def_mod(int (mol2->COM().x/(2*PI/regio_res)), regio_res);
                }
                // y divide
                else {
                    index = pos_def_mod(int (mol2->COM().y/regio_res), regio_res);
                }
                regio_c1.at(index).at(k).add(legendre(1,phi));
                regio_c2.at(index).at(k).add(legendre(2,phi));
                regio_MSD.at(index).at(k).add(pow(displacement,2));
            }
            k++;
        }
        
        
        if (time_structure || print){
            dyn_queue r = dyn_queue(&mol);
            molecule *mol2 = r.pop();
            vect direction;
        
            // All molecules within cutoff
            while ( mol2 && direction.length() < radial_cutoff){
                direction = frame->direction(mol.COM(), mol2->COM());
                if (direction.length() > 0){
                    radial.add(direction.length());
                    radial2d_rel.at(pos_def_mod(int((direction.angle()-mol.get_orientation())/dtheta), theta_res)).add(direction.length());
                    radial2d_abs.at(pos_def_mod(int((direction.angle())/dtheta), theta_res)).add(direction.length());
                    for (auto p: mol2->atoms){
                        direction = frame->direction(mol.COM(), p->pos_vect());
                        if (p->type == 1){
                            radial2d_large.at(pos_def_mod(int((direction.angle()-mol.get_orientation())/dtheta), theta_res)).add(direction.length());
                        }
                        else if (p->type == 2){
                            radial2d_small.at(pos_def_mod(int((direction.angle()-mol.get_orientation())/dtheta), theta_res)).add(direction.length());
                        }
                        for (auto p1: mol.atoms){
                            direction = frame->direction(p1->pos_vect(), p->pos_vect());
                            radial_part.add(direction.length());
                            radial2d_part.at(pos_def_mod(int(direction.angle()/dtheta), theta_res)).add(direction.length());
                            radial_part.add(direction.length());
                        }
                    }
                }
                mol2 = r.pop();
            }
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
        if (regio){
            for (int i = 0; i < regio_res; i++){
                collate_regio_c1.at(i)[steps].add(regio_c1.at(i).at(k).get_mean());
                collate_regio_c2.at(i)[steps].add(regio_c2.at(i).at(k).get_mean());
                collate_regio_MSD.at(i)[steps].add(regio_MSD.at(i).at(k).get_mean());
            }
        }
    }
    
    if (time_structure){
        print_time_distribution(&short_order, frame->timestep, &short_order_file);
        order_file << frame->timestep << "," << hexatic_order.get_mean()\
                << "," << circle_order.get_mean() \
                << "," << orientational_order.get_mean() << endl;
    }
    
    double max_alpha = 0;
    int max_alpha_time = 0;
    if (print && key_frames.size() > 0){
        // num neighbours
        print_distribution<int>(&num_neigh, "stats/neighbours.dat");
        
        // Print num contacts
        print_distribution<int>(&num_contact, "stats/contacts.dat");
        
        // Print num contacts
        print_distribution<int>(&pairing, "stats/pairing.dat");
        
        // Print radial distribution
        print_radial_distribution(&radial, "radial_dist.dat", frame->num_mol(), frame->get_area());
        print_radial_distribution(&radial_part, "radial_part.dat", frame->num_atoms(), frame->get_area());
        print_radial2d_distributions("radial2d.dat", frame, &radial2d_abs, {&radial2d_rel, &radial2d_small, &radial2d_large, &radial2d_part});
        
        if (m_orient){
            ofstream orient;
            orient.open("velocity.dat");
            int i = 0;
            for (auto v: velocity){
                orient << dtheta*i << " " << v.get_mean() << endl;
                i++;
            }
            orient.close();
        }
        
        // Order Parameters
        
        ofstream ordering;
        ordering.open("order.log");
        ordering << "Hexatic: " << hexatic_order.get_mean() << endl;
        ordering << "Circle-order: " << circle_order.get_mean() << endl;
        ordering << "Frac-6-fold: " << num_neigh.fraction_at(6) << endl;
        ordering << "Orientational: " << orientational_order.get_mean() << endl;
        ordering << "Structure-factor: " << max_structure_factor(get_radial_distribution(&radial, frame->num_mol(), frame->get_area()), frame->get_density(), 0.015) << endl;
        ordering.close();

        // Print short range order
        //print_distribution<double>(&short_order, "short_order_dist.dat");
        
        // Displacement
        double alpha;
        MSD_file << "Timestep,MSD,MFD,\\alpha" << endl;
        for (auto d: collate_MSD){
            alpha = collate_MFD.at(d.first).get_mean() / (2*pow(d.second.get_mean(),2)) - 1;
            if (alpha > max_alpha){
                max_alpha = alpha;
                max_alpha_time = d.first;
            }
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
        struct_file << "Timestep,F,chi"  << endl;
        for (auto c: collate_struct){
            struct_file << c.first << "," << c.second.get_mean() \
            << "," << relax_time(c.second.get_variance()) << endl;
            //cout << c.second.get_mean() << " " << c.second.get_variance() << " " << c.second.get_count() << endl;
            if (c.second.get_mean() < 1/CONST_E && ts == 0){
                ts = c.first;
            }
        }
        print_relax_time("Struct-relax: ", ts);
        // Calculate diffusion constant
        my_mean diff_const;
        for (auto m1: collate_MSD){
            for (auto m2: collate_MSD){
                if (m1.first > relax_time(ts) && m1.first < m2.first){
                    diff_const.add((m2.second.get_mean() - m1.second.get_mean())/(4*(m2.first -m1.first)*STEP_SIZE), m2.first-m1.first);
                }
            }
        }
        cout << "Diffusion-constant: " << print_relax_time(diff_const.get_mean()) << endl;
        if (regio){
            // Regio relaxations
            regio_file.open("regio.csv");
            regio_file << "Timestep,Position,MSD,R1,R2,circle,orientation" << endl;
            
            ofstream regio_relax;
            regio_relax.open("regio-relax.csv");
            regio_relax << "Pos,t1,t2,circle,orientation" << endl;
            
            int regio_t1, regio_t2;
            double delta;
            int first = 0;
            if (frame->get_a() > frame->get_height()){
                delta = frame->get_a()/regio_res;
            }
            else {
                delta = frame->get_height()/regio_res;
            }
            for (int i = 0; i < regio_res; i++){
                regio_t1 = 0;
                regio_t2 = 0;
                for (auto &c: collate_regio_MSD.at(i)){
                    if (first == 0){
                        first  = c.first;
                    }
                    regio_file << c.first << "," << i*delta << "," << c.second.get_mean() \
                    << "," << collate_regio_c1.at(i).at(c.first).get_mean() \
                    << "," << collate_regio_c2.at(i).at(c.first).get_mean() \
                    << "," << regio_circle.at(i).get_mean() \
                    << "," << regio_orientation.at(i).get_mean() << endl;
                    
                    if (collate_regio_c1.at(i).at(c.first).get_mean() < 1/CONST_E && regio_t1 == 0){
                        regio_t1 = c.first;
                    }
                    else if (collate_regio_c2.at(i).at(c.first).get_mean() < 1/CONST_E && regio_t2 == 0){
                        regio_t2 = c.first;
                    }
                }
                regio_file << endl;
                regio_relax << i*delta << "," << print_relax_time(regio_t1) \
                << "," << print_relax_time(regio_t2) \
                << "," << regio_circle.at(i).get_mean()
                << "," << regio_orientation.at(i).get_mean()
                << endl;
            }
        }
        print_rot_diff(key_frames, frame);
        
        print_distribution(&pair_contact, "stats/pair_contact.dat");
        print_distribution(&pair_neigh, "stats/pair_neigh.dat");
        print_distributions("stats/diff_contact.dat", MAX_MOL_CONTACTS, {&num_contact_mob[0],&num_contact_mob[1],&num_contact_mob[2],&num_contact_mob[3]});
        print_distributions("stats/diff_neigh.dat", MAX_MOL_CONTACTS, {&num_neigh_mob[0],&num_neigh_mob[1],&num_neigh_mob[2],&num_neigh_mob[3]});
        
        int key = 0;
        for (auto f: key_frames){
            key  = f->timestep;
            if (f->timestep > max_alpha_time){
                break;
            }
        }
        print_moved(key_frames.front(), key_frames.at(key));
    }
    if  (print){
        print_frame(frame);
        
    }

    return 0;
}
