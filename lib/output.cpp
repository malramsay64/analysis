//
//  output.cpp
//  analysis
//
//  Created by Malcolm Ramsay on 7/12/2014.
//  Copyright (c) 2014 Malcolm Ramsay. All rights reserved.
//

#include "output.h"

using namespace std;


int print_prop(vector<Frame *> frames, double (*func)(Frame *, Frame *), string fname){
    ofstream file;
    file.open(fname.c_str(), ios::out);
    file << "Timestep,function" << endl;
    vector<Frame *>::iterator f;
    for (f = frames.begin(); f != frames.end(); f++){
        file << (*f)->timestep << "," \
        << (*func)(frames.front(),(*f)) << endl;
    }
    file.close();
    return 0;

}

int MSD(vector<Frame *> frames){
    return print_prop(frames, &par_MSD, "msd.csv");
}

int diffusion(vector<Frame *> frames){
    return print_prop(frames, &par_diffusion_constant, "diffusion.csv");
}

int rotation(vector <Frame *> frames){
    ofstream file;
    file.open("rotation.csv", ios::out);
    file << "Timestep,C1,C2,C3,C4" << endl;
    double c1,c2,c3,c4, tend;
    int t1 = 0, t2 = 0, t3 = 0, t4 = 0;
    vector<Frame *>::iterator f;
    for (f = frames.begin()+1; f != frames.end() ; f++){
        c1 = par_rotational_relaxation(1,frames.front(),*f);
        c2 = par_rotational_relaxation(2,frames.front(),*f);
        c3 = par_rotational_relaxation(3,frames.front(),*f);
        c4 = par_rotational_relaxation(4,frames.front(),*f);
        file << (*f)->timestep << "," << c1 << "," << c2 << "," << c3 << "," << c4 << endl;
        
        // Relaxation Times
        if (c1 < 0 && t1 == 0){
            t1 = (*f)->timestep;
        }
        if (c2 < 0.25 && t2 == 0){
            t2 = (*f)->timestep;
        }
        if (c3 < 0 && t3 == 0){
            t3 = (*f)->timestep;
        }
        if (c4 < 0.1406 && t4 == 0){
            t4 = (*f)->timestep;
        }
    }
    tend = frames.back()->timestep;
    cout << "C1 Relaxation: " << t1/tend << endl;
    cout << "C2 Relaxation: " << t2/tend << endl;
    cout << "C3 Relaxation: " << t3/tend << endl;
    cout << "C4 Relaxation: " << t4/tend << endl;
    file.close();
    return 0;
}

// Outputs the distance and angle of each particle
int order_parameter(string o, int reference, Frame *frame){
    ofstream file;
    char fname[30];
    snprintf(fname,30, "%s/%010i.csv", o.c_str(), frame->timestep);
    file.open(fname, ios::out);
    order_parameter(reference, &file, frame);
    file.close();
    return 0;
}

void print(double (Frame::*f)(vector<int> *), Frame *frame, string str){
    char fname[30];
    snprintf(fname, 30, "stats/%s.csv", str.c_str());
    vector<int> dist(15,0);
    ofstream file;
    file.open(fname, ios::out);
    cout << str << ": " << (frame->*f)(&dist)/10 << endl;
    vector<int>::iterator i;
    for (i = dist.begin(); i != dist.end(); i++){
        file << i - dist.begin() << " " << *i << endl;
    }
    file.close();
}

int stats(Frame *frame){
    par_neigh(frame);
    print(&Frame::num_contacts, frame , "num_contacts");
    print(&Frame::num_neighbours, frame, "num_neighbours");
    print(&Frame::pairing, frame, "pairing");
    return 0;
}

int print(Frame *frame, ofstream *file){
    vector<particle>::iterator p;
    *file << "ITEM: TIMESTEP" << endl;
    *file << frame->timestep << endl;
    *file << "ITEM: NUMBER OF ATOMS" << endl;
    *file << frame->num_atoms() << endl;
    *file << "ITEM: BOX BOUNDS pp pp pp" << endl;
    *file << frame->xmin() << " " << frame->xmax() << endl;
    *file << frame->ymin() << " " << frame->ymax() << endl;
    *file << frame->zmin() << " " << frame->zmax() << endl;
    *file << "ITEM: ATOMS id mol type x y z vx vy vz" << endl;
    for (p = frame->particles.begin(); p != frame->particles.end(); p++){
        //*file << (*p).id << " " << (*p).molid << " " << (*p).type << " " << (*p).pos_vect() << " 0" << endl;
        int m = 0;
        if ((*p).get_traversed()){
            m = 1;///(*p).molid;
        }
        *file << (*p).id <<  " 1 " << (*p).type << " " << (*p).pos_vect()*frame->size() << " 0 0 " << "0" << " " << m << endl;
    }
    //file->close();
    return 0;
}

int print_angles(Frame * frame){
    ofstream file;
    char fname[40];
    snprintf(fname,40, "trj_contact/%010i-angles.csv", frame->timestep);
    file.open(fname, ios::out);
    double range = 2*PI;
    int n_bins = 720;
    int bin;
    vector<int> absolute(n_bins, 0);
    vector<int> relative(n_bins, 0);
    vector<molecule>::iterator mol;
    vector<molecule>::iterator mol2;
    for (mol = frame->molecules.begin(); mol != frame->molecules.end(); mol++){
        // Angle relative to frame
        bin = (int) (angle(&(*mol), frame) / (range/(n_bins)));
        absolute.at(bin)++;
        for (mol2 = mol+1; mol2 != frame->molecules.end(); mol2++){
            bin = ((int) (angle(&(*mol), frame) - angle(&(*mol2), frame) / (range/(n_bins))) + n_bins) % n_bins;
            relative.at(bin)++;
        }
    }
    for (int i = 0; i < n_bins; i++){
        file << range/(n_bins) * (i) * 180/PI << ","\
        << absolute.at(i)/(float) frame->num_mol() << ","\
        << relative.at(i)/(((float) frame->num_mol())/2 * frame->num_mol()) << endl;
    }
    return 0;
}

int print_gnuplot(Frame * frame){
    par_neigh(frame);
    colourise(frame);
    ofstream file;
    char fname[30];
    snprintf(fname,30, "trj_contact/%010i.dat", frame->timestep);
    file.open(fname);
    vector<particle>::iterator p;
    int molid = 1;
    for (p = frame->particles.begin(); p != frame->particles.end(); p++){
        if (molid != (*p).molid){
            file << endl;
            molid = (*p).molid;
        }
        file << (*p).pos_vect()*frame->size() << " " << (*p).radius << " " << get_colour(&(*p), frame) << endl;
        //file << (*p).pos_vect()*frame->size() << " " << (*p).radius << " " << (*p).type << endl;
    }
    file.close();
    return 0;
}


