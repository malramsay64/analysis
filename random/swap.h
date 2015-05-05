//
//  swap.h
//  analysis
//
//  Created by Malcolm Ramsay on 27/04/2015.
//  Copyright (c) 2015 Malcolm Ramsay. All rights reserved.
//

#ifndef __analysis__swap__
#define __analysis__swap__

#include <stdio.h>
#include <stdlib.h>
#include "frame.h"
#include "neighbours.h"
#include "neighbours.h"
#include "frame.h"
#include "input.h"
#include "output.h"
#include <string>
#include <sstream>
#include <iostream>

#include <set>

int swap_neighbours(Frame * frame, int molid);
int update_neighbours(molecule * mol, Frame * frame);

#endif /* defined(__analysis__swap__) */
