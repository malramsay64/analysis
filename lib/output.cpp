//
//  output.cpp
//  analysis
//
//  Created by Malcolm Ramsay on 9/02/2015.
//  Copyright (c) 2015 Malcolm Ramsay. All rights reserved.
//

#include "output.h"

using namespace std;

int print_map(std::map<int, my_mean> map, std::ofstream * file){
    for (auto m: map){
        *file << m.first << "," << m.second.get_mean() << endl;
    }
    return 0;
}

int print_mol(ostream *file, molecule *mol, Frame *frame){
    vect d, com;
    
    com = frame->cartesian(mol->COM());
    com = wrap_x(com, frame->get_a());
    
    particle * p;
    
    for (int i = 0; i < mol->atoms.size(); i++){
        p = mol->atoms.at((i+1) % mol->atoms.size());
        d = frame->cartesian(direction(mol->COM(), p->pos_vect()));
        *file << com + d << " " << p->radius << " " << mol->get_colour() << " " << mol->id << endl;
    }
    *file << endl;
    
    return 0;
}

int print_short_order(ofstream * file, molecule * mol,  Frame * frame){
    int type;
    vect d;
    double theta;
    for (auto mol2: mol->my_neighbours){
        type = order_type(mol, mol2, frame);
        for (int k = 0; k < mol2->nump(); k++){
            d = frame->direction(mol2->atom_pos(k),mol->COM());
            theta = atan2(d) + 5*PI/2 - atan2(orientation(mol, frame));
            *file << theta << "," << d.length() << "," << 0.04 << "," << type << endl;
        }
    }
    return 0;
}

int print_radial_distribution(distribution<int> *d, string filename, int nmol, double frame_area){
    ofstream file;
    file.open(filename.c_str());
    double area;
    double density = 2*nmol/frame_area;
    for (int i = 0; i < d->get_size(); i++){
        area = PI*(i*d->get_delta_r())*d->get_delta_r();
        file << i*d->get_delta_r() << " " << d->at(i)/(area*nmol*density) << endl;
    }
    return 0;
}

int print_relax_time(string s, int t){
    if (t){
        cout << s << t << endl;
    }
    else {
        cout << s << 1 << endl;
    }
    return 0;
}

int print_frame(Frame * frame){
    char fname[40];
    ofstream gnuplot;
    ofstream complot;
    snprintf(fname,40, "trj_contact/%010i.dat", frame->timestep);
    gnuplot.open(fname);
    complot.open("complot.dat");
    gnuplot << frame->get_a() << " " << frame->get_height() << endl << endl;
    complot << frame->get_a() << " " << frame->get_height() << endl << endl;
    
    for (auto &m: frame->molecules){
        vect d, com;
        com = frame->cartesian(m.COM());
        com = wrap_x(com, frame->get_a());
        
        particle * p;
        complot << com << endl;
        for (int i = 0; i < m.atoms.size(); i++){
            p = m.atoms.at((i+2) % m.atoms.size());
            d = frame->cartesian(direction(m.COM(), p->pos_vect()));
            gnuplot << com + d << " " << p->radius << " " << mol_colour(&m,frame) << " " << m.id << endl;
        }
        gnuplot << endl;
        //print_mol(&gnuplot, &m, frame);
    }
    
    return 0;
}