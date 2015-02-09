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
