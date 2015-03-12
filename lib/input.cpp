//
//  input.cpp
//  analysis
//
//  Created by Malcolm Ramsay on 7/12/2014.
//  Copyright (c) 2014 Malcolm Ramsay. All rights reserved.
//

#include "input.h"
#include <algorithm>
#include <sstream>

using namespace std;

int read_data(std::ifstream *myfile, Frame *frame){
    double x[3], y[3], z[3];
    double zp;
    int num_atoms, timestep;
    int mols = 0;
    string line;
    double a,b,theta;
    // TIMESTEP
    if (!getline(*myfile, line)){
        throw 20;
    }
    *myfile >> timestep;
    frame->set_timestep(timestep);
    // NUMBER OF ATOMS
    getline(*myfile, line);
    getline(*myfile, line);
    *myfile >> num_atoms;
    frame->set_atoms(num_atoms);
    // BOX BOUNDS
    getline(*myfile, line);
    getline(*myfile, line);
    
    // X
    getline(*myfile, line);
    stringstream sx(line);
    x[2] = 0;
    sx >> x[0] >> x[1] >> x[2];
    a = x[1]-fabs(x[2]) - x[0];
    
    // Y
    getline(*myfile, line);
    stringstream sy(line);
    y[2] = 0;
    sy >> y[0] >> y[1] >> y[2];
    b = y[1]-fabs(y[2]) - y[0];
    
    // Z
    getline(*myfile, line);
    stringstream sz(line);
    z[2] = 0;
    sz >> z[0] >> z[1] >> z[2];
    
    // Crystal parameters
    theta = atan(b/fabs(x[2]));
    if (x[2] < 0){
        theta = PI-theta;
    }
    b = b/sin(theta);
    frame->set_crys(a, b, theta);
    
    //ITEM: ATOMS
    getline(*myfile, line);
    particle *p;
    p = new particle();
    double pos[2] = {0,0};
    for (int i=0; i < frame->num_atoms(); i++){
        getline(*myfile, line);
        stringstream sp(line);
        sp >> p->id >> p->molid >> p->type >> p->radius >> pos[0] >> pos[1] >> zp;
        // Shift axes of box to 0
        p->set_pos(vect(pos[0], pos[1]) - vect(x[0], y[0]));
        
        // Convert to fractional coordinates
        p->set_pos(frame->fractional(p->pos_vect()));
        
        // Radius read in as diameter
        p->radius = p->radius/2;
        frame->add_particle(p);
        // Get num mols
        if (p->molid > mols){
            mols = p->molid;
        }
    }
    delete p;
    frame->set_num_mol(mols);
    sort(frame->particles.begin(), frame->particles.end());
    vector<molecule>::iterator j;
    for (j = frame->molecules.begin(); j != frame->molecules.end(); ++j){
        (*j).id = (int) (j-frame->molecules.begin()) + 1;
    }
    
    vector<particle>::iterator i;
    for (i = frame->particles.begin(); i != frame->particles.end(); ++i){
        frame->add_link(i->molid-1, i->id-1);
    }
    
    // Put all particles from molecule on same frame based on COM
    // Update angle of molecule
    for (auto &m: frame->molecules){
        m.same_period();
        m.set_orientation(angle(&m,frame));
    }
    return 0;
}

int update(std::ifstream *myfile, Frame *frame){
    if (frame->num_atoms() == 0){
        read_data(myfile, frame);
        return 0;
    }
    double x[3], y[3], z[3];
    double zp;
    int timestep;
    string line;
    double a,b,theta;
    // TIMESTEP
    if (!getline(*myfile, line)){
        throw 20;
    }
    *myfile >> timestep;
    frame->set_timestep(timestep);
    // NUMBER OF ATOMS
    getline(*myfile, line);
    getline(*myfile, line);
    //*myfile >> num_atoms;
    //frame->set_atoms(num_atoms);
    // BOX BOUNDS
    getline(*myfile, line);
    getline(*myfile, line);
    // X
    getline(*myfile, line);
    stringstream sx(line);
    x[2] = 0;
    sx >> x[0] >> x[1] >> x[2];
    a = x[1]-fabs(x[2]) - x[0];
    
    // Y
    getline(*myfile, line);
    stringstream sy(line);
    y[2] = 0;
    sy >> y[0] >> y[1] >> y[2];
    b = y[1]-fabs(y[2]) - y[0];
    
    // Z
    getline(*myfile, line);
    stringstream sz(line);
    z[2] = 0;
    sz >> z[0] >> z[1] >> z[2];
    
    // Crystal parameters
    theta = atan(b/fabs(x[2]));
    if (x[2] < 0){
        theta = PI-theta;
    }
    b = b/sin(theta);
    frame->set_crys(a, b, theta);
    
    //ITEM: ATOMS
    getline(*myfile, line);
    particle *p;
    p = new particle();
    double pos[2] = {0,0};
    for (int i=0; i < frame->num_atoms(); i++){
        getline(*myfile, line);
        stringstream sp(line);
        sp >> p->id >> p->molid >> p->type >> p->radius >> pos[0] >> pos[1] >> zp;
        // Update coordinates
        frame->particles.at(p->index()).set_pos(frame->fractional(vect(pos[0], pos[1]) - vect(x[0], y[0])));

    }
    delete p;
    for (auto &m: frame->molecules){
        m.calc_COM();
        m.delete_neighbours();
        m.same_period();
        m.update_orientation(angle(&m,frame));
    }
    return 0;
}


int skip_frame(std::ifstream *myfile){
    string line;
    int num_atoms;
    if (!getline(*myfile, line)){
        throw 20;
    }
    // NUMBER OF ATOMS
    getline(*myfile, line);
    getline(*myfile, line);
    *myfile >> num_atoms;
    // BOX BOUNDS
    getline(*myfile, line);
    getline(*myfile, line);
    // X
    getline(*myfile, line);
    // Y
    getline(*myfile, line);
    // Z
    getline(*myfile, line);
    // Crystal parameters
    //ITEM: ATOMS
    getline(*myfile, line);
    for (int i=0; i < num_atoms; i++){
        getline(*myfile, line);
    }
    return 0;
}



