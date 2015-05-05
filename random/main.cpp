//
//  main.cpp
//  random
//
//  Created by Malcolm Ramsay on 27/04/2015.
//  Copyright (c) 2015 Malcolm Ramsay. All rights reserved.
//

#include "swap.h"

using namespace std;

double dtheta = 0;

int main(int argc, const char * argv[]) {
    string crys_fname = "Snowman-0.637556-1.637556-p2.lammpstrj";
    ifstream crys_file;
    crys_file.open(crys_fname.c_str());
    Frame * crystal = new Frame;
    read_data(&crys_file,crystal);
    vector<vector<int>> neigh_list = vector<vector<int>>(crystal->num_mol(), vector<int>());
    find_neighbours(crystal, &neigh_list);
    
    int steps = 10000;
    int swaps = 0;
    srand(0);
    for (int i = 0; i < steps; i++){
        int r = rand() % crystal->num_mol();
        //cout << r << endl;
        swaps += swap_neighbours(crystal, r);
    }
    print_frame(crystal);
    cout << swaps << " " << steps << endl;
    return 0;
}
