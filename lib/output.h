#include "frame.h"
#include "parallel.h"
#include "neighbours.h"

#ifndef MY_OUTPUT_H
#define MY_OUTPUT_H

int MSD(int num_frames, std::vector<Frame *> frames);
int diffusion(int num_frames, std::vector<Frame *> frames);
int rotation(int num_frames, std::vector <Frame *> frames);
int order_parameter(std::string o,int reference, Frame *frame);
//void print(double (Frame::* f)(vector<int> *), Frame *frame, char *fname);
int stats(Frame *frame);
int print(Frame *frame, std::ofstream *file);
int print_angles(Frame * frame);
int print_gnuplot(Frame * frame);

#endif 
