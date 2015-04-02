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

vector<double> get_radial_distribution(distribution<int> *d, int nmol, double frame_area){
    double area;
    double density = 2*nmol/frame_area;
    vector<double> radial;
    for (int i = 0; i < d->get_size(); i++){
        area = PI*(i*d->get_delta_r())*d->get_delta_r();
        radial.push_back(d->at(i)/(area*nmol*density));
    }
    return radial;
}


int print_relax_time(string s, int t){
    cout << s << print_relax_time(t) << endl;
    return 0;
}

int relax_time(int t){
    if (t == 0){
        return INFINITY;
    }
    return t;
}

double relax_time(double t){
    if (t == 0){
        return INFINITY;
    }
    return t;
}

string print_relax_time(int t){
    if (t){
        stringstream ss;
        ss << setprecision(5) << scientific << double (t);
        return ss.str();
    }
    return "NAN";
}

string print_relax_time(double t){
    if (t != 0){
        stringstream ss;
        ss << setprecision(5) << scientific << double (t);
        return ss.str();
    }
    return "NAN";
}

int print_rot_diff(vector<Frame *> key_frames, Frame * frame){
    ofstream rot_diff;
    rot_diff.open("rot_diff.csv");
    rot_diff << "Molecule,rotation,diffusion" << endl;
    double rot, diff, total_rot;
    int i;
    for (auto m: frame->molecules){
        total_rot = 0;
        i = m.index();
        for (auto key: key_frames){
            rot = key->at(i).get_rotation();
            diff = frame->dist(key->at(i).COM(),key_frames.front()->at(i).COM());
            rot_diff << m.id << "," << fabs(rot) << "," << diff << endl;
        }
        rot = m.get_rotation();
        diff = frame->dist(m.COM(), key_frames.front()->at(i).COM());
        rot_diff << m.id << "," << fabs(rot) << "," << diff << endl;
    }
    return 0;
}


int print_moved(Frame * init, Frame * final){
    ofstream file;
    file.open("moved.dat");
    file << init->get_a() << " " << init->get_height() << endl << endl;
    
    vect com1, com2;
    double rotation;
    for (auto &m: init->molecules){
        com1 = m.COM();
        com2 = final->molecules.at(m.index()).COM();
        rotation = final->molecules.at(m.index()).get_rotation();
        
        file << init->cartesian(com1) << " " << init->direction(com1, com2) << " " << rotation << endl;
    }
    file.close();
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
        
        particle * p;
        complot << com << " " << com_colour(&m, frame) << endl;
        for (int i = 0; i < m.atoms.size(); i++){
            p = m.atoms.at((i+2) % m.atoms.size());
            d = frame->cartesian(direction(m.COM(), p->pos_vect()));
            gnuplot << com + d << " " << p->radius << " " << \
            m.get_orientation() << " " << m.num_neighbours() << " " << circle_ordering(&m) \
            << " " << short_ordering(&m,frame) << " " << m.id << endl;
        }
        gnuplot << endl;
    }
    gnuplot.close();
    complot.close();
    return 0;
}

int print_motion(vector<Frame *> key_frames){
    Frame * last;
    vect com1,com2;
    vect orient1, orient2;
    for (auto key: key_frames){
        if (last){
            for (auto m1: key->molecules){
                com1 = m1.COM();
                com2 = last->molecules.at(m1.index()).COM();
                orient1 = m1.get_orient_vect();
                orient2 = last->molecules.at(m1.index()).get_orient_vect();
            }
        }
    }
    return 0;
}