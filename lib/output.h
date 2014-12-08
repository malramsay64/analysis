//
//  output.h
//  analysis
//
//  Created by Malcolm Ramsay on 7/12/2014.
//  Copyright (c) 2014 Malcolm Ramsay. All rights reserved.
//

#ifndef MY_OUTPUT_H
#define MY_OUTPUT_H

#include "frame.h"
#include "parallel.h"
#include "neighbours.h"

int print_prop(std::vector<Frame *> frames, double (*func)(Frame *, Frame *), std::string fname);
int MSD(std::vector<Frame *> frames);
int diffusion(std::vector<Frame *> frames);
int rotation(std::vector <Frame *> frames);
int order_parameter(std::string o, int reference, Frame *frame);
//void print(double (Frame::* f)(vector<int> *), Frame *frame, char *fname);
int stats(Frame *frame);
int print(Frame *frame, std::ofstream *file);
int print_angles(Frame * frame);
int print_gnuplot(Frame * frame);

#endif
