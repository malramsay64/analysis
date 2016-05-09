//
//  input.h
//  analysis
//
//  Created by Malcolm Ramsay on 7/12/2014.
//  Copyright (c) 2014 Malcolm Ramsay. All rights reserved.
//

#include <iostream>
#include <string>
#include <getopt.h>

#include "Particle.h"
#include "Frame.h"



#ifndef INPUT_H
#define INPUT_H

// Variables

Frame read_frame(std::istream &is);
void useage();

static std::string options{"i:s:k:rtoqfmd"};
#endif /* defined(INPUT_H) */

