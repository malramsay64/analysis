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
#include <exception>


void useage(){
    std::cout <<\
                  "useage: analysis -i <filename> [-qfmdrto] [-s <steps>] [-k <steps]"\
                  ""\
                  "-i <filename>"\
                  "      filename to read"\
    << std::endl;
    /* Read arguments
* Valid arguments are:
* -i <file>  input filename
* -q look at quenched configuration
* -s <n> look at every n steps
* -f fast, only first and last
* -m movie, output movie
* -d distance moved in quench
* -k <n> key frame rate
* -r regio data
* -t time depedent structure
* -o orientation of motion
*/
}

static struct option long_options[] = {
    {}
};


struct options {
    std::string inFname;
    bool quenched{false};
    bool fast{false};
    bool movie{false};
};

Frame read_frame(std::istream &is) {
    FrameVars fv;
    std::string s;
    // TIMESTEP
    is >> s;
    is >> fv.timestep;
    // NUMBER OF ATOMS
    is >> s;
    is >> fv.num_atoms;
    // BOX BOUNDS
    is >> s;
    Vector x, y, z, tilt;
    is >> x;
    is >> y;
    is >> z;
    is >> tilt;
    fv.a = x[1]-std::fabs(x[2]) - x[0];
    fv.b = y[1]-std::fabs(y[2]) - y[0];
    fv.theta = std::atan(fv.b/std::fabs(x[2]));
    if (x[2] < 0){
        fv.theta = PI-fv.theta;
    }
    fv.b = fv.b/std::sin(fv.theta);

    Frame frame{fv};
    int mols;
    // ITEM: ATOMS
    is >> s;
    for (auto &p: frame.particles){
        is >> p;
        p.molid > mols ? mols = p.molid: mols;
    }

    frame.set_num_mol(mols);
    std::sort(frame.particles.begin(), frame.particles.end());
    std::vector<Molecule>::iterator j;
    for (int i=0; i<frame.molecules.size(); i++){
        frame.molecules[i].id = i+1;
    }

    for (auto &p: frame.particles){
        frame.add_link(p.mol_index(), p.index());
    }

    // Put all particles from molecule on same frame based on COM
    // Update angle of molecule
    for (auto &m: frame.molecules){
        m.same_period();
    }
}
