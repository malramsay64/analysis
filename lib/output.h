//
//  output.h
//  analysis
//
//  Created by Malcolm Ramsay on 9/02/2015.
//  Copyright (c) 2015 Malcolm Ramsay. All rights reserved.
//

#ifndef __analysis__output__
#define __analysis__output__

#include <stdio.h>
#include <map>
#include <limits>
#include <iomanip>
#include "my_mean.h"
#include "neighbours.h"

int print_map(std::map<int, my_mean>, std::ofstream *);
int print_short_order(std::ofstream * file, molecule * mol,  Frame * frame);
int print_mol(std::ostream *os, molecule *mol, Frame *frame);
int print_frame(Frame * frame);
int print_radial_distribution(distribution<int> *, std::string, int, double);
int print_relax_time(std::string s, int t);
int relax_time(int t);
std::string print_relax_time(int t);
int print_moved(Frame * init, Frame * final);
std::vector<double> get_radial_distribution(distribution<int> *, int, double);
int print_rot_diff(std::vector<Frame *> key_frames, Frame * frame);

#endif /* defined(__analysis__output__) */
