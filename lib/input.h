#include <iostream>
#include <string>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <stdio.h>

#include "particle.h"
#include "frame.h"

#define FALSE 1
#define TRUE 0

#ifndef PI
#define PI 3.14159265358979323846
#endif

#ifndef INPUT_H
#define INPUT_H

// Variables
int read_data(std::ifstream *myfile, Frame *frame);
int read_data(std::ifstream *myfile);

#endif 
