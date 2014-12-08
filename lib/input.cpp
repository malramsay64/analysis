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
    int a, timestep;
    int mols = 0;
    string line;
    // TIMESTEP
    if (!getline(*myfile, line)){
        //cout << line << endl;
        throw 20;
    }
    *myfile >> timestep;
    frame->set_timestep(timestep);
    // NUMBER OF ATOMS
    getline(*myfile, line);
    getline(*myfile, line);
    *myfile >> a;
    frame->set_atoms(a);
    // BOX BOUNDS
    getline(*myfile, line);
    getline(*myfile, line);
    
    
    // X
    getline(*myfile, line);
    stringstream sx(line);
    x[2] = 0;
    sx >> x[0] >> x[1] >> x[2];
    frame->setx(x[0], x[1]-abs(x[2]));
    
    // Y
    getline(*myfile, line);
    stringstream sy(line);
    y[2] = 0;
    sy >> y[0] >> y[1] >> y[2];
    frame->sety(y[0], y[1] - abs(y[2]));
    
    // Z
    getline(*myfile, line);
    stringstream sz(line);
    z[2] = 0;
    sz >> z[0] >> z[1] >> z[2];
    frame->setz(z[0], z[1]-abs(z[2]));
    
    //ITEM: ATOMS
    getline(*myfile, line);
    particle *p;
    p = new particle();
    for (int i=0; i < a; i++){
        *myfile >> p->id >> p->molid >> p->type >> p->radius >> p->pos[XX] >> p->pos[YY] >> zp;
        getline(*myfile, line);
        // Shift axes of box to 0
        p->set_xpos(p->xpos() - frame->xmin());
        p->set_ypos(p->ypos() - frame->ymin());
        
        // Scale box length to 2*PI
        p->set_xpos(p->xpos()/(frame->xlength()));
        p->set_ypos(p->ypos()/(frame->ylength()));
        // Radius read in as diameter
        p->radius = p->radius/2;
        frame->add_particle(p);
        // Get num mols
        if (p->molid > mols){
            mols = p->molid;
        }
    }
    delete p;
    //cout << mols << endl;
    frame->set_num_mol(mols);
    sort(frame->particles.begin(), frame->particles.end());
    vector<molecule>::iterator j;
    for (j = frame->molecules.begin(); j != frame->molecules.end(); ++j){
        (*j).id = (j-frame->molecules.begin()) + 1;
    }
    
    vector<particle>::iterator i;
    for (i = frame->particles.begin(); i != frame->particles.end(); ++i){
        frame->add_link(i->molid-1, i->id-1);
        //cout << i->molid << endl;
    }
    
    for (j = frame->molecules.begin(); j != frame->molecules.end(); ++j){
        (*j).same_period();
        // Put all particles on same frame based on COM
    }
    //cout << molecules.size()<< endl;
    return 0;
}


