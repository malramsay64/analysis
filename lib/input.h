//
//  input.h
//  analysis
//
//  Created by Malcolm Ramsay on 7/12/2014.
//  Copyright (c) 2014 Malcolm Ramsay. All rights reserved.
//

#include <iostream>
#include <string>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <stdio.h>

#include "particle.h"
#include "frame.h"
#include "constants.h"



#ifndef INPUT_H
#define INPUT_H

// Variables
int read_data(std::ifstream *myfile, Frame *frame);

int skip_frame(std::ifstream *myfile);

#endif /* defined(INPUT_H) */

