//
//  distribution.cpp
//  analysis
//
//  Created by Malcolm Ramsay on 6/02/2015.
//  Copyright (c) 2015 Malcolm Ramsay. All rights reserved.
//

#include "distribution.h"

using namespace std;

int print_radial_distribution(distribution<int> *d, string filename, int nmol, double frame_area){
    ofstream file;
    file.open(filename.c_str());
    double area;
    double density = nmol/frame_area;
    for (int i = 0; i < d->get_size(); i++){
        area = PI*(i*d->get_delta_r())*d->get_delta_r();
        file << i*d->get_delta_r() << " " << d->at(i)/(area*nmol*density) << endl;
    }
    return 0;
}


