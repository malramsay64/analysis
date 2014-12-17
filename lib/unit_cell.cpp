//
//  unit_cell.cpp
//  analysis
//
//  Created by Malcolm Ramsay on 18/12/2014.
//  Copyright (c) 2014 Malcolm Ramsay. All rights reserved.
//

#include "unit_cell.h"

using namespace std;

static double eps = 5*PI/180;

int unit_cell(Frame * frame) {
    vector<molecule>::iterator m;
    molecule mol = frame->molecules.front();
    vect d;
    angle_list angles;
    int count = 0;
    for (m = frame->molecules.begin(); m != frame->molecules.end(); m++){
        d = frame->direction((*m).COM(), mol.COM());
        count++;
        if (fabs(angle(&(*m),frame) - angle(&mol, frame)) < eps){
            angles.push(atan2(d),d.length());
        }
        if (fabs(fabs(atan2(d))-PI/2) < eps && fabs(angle(&(*m),frame) - angle(&mol, frame)) < eps){
            //angles.push(atan2(d),d.length());
            //cout << d.length() << endl;
        }
    }
    ostream * file = &cout;
    angles.print(file);
    //cout << count << endl;
    return 0;
}