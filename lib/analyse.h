//
//  analyse.h
//  analysis
//
//  Created by Malcolm Ramsay on 14/01/2015.
//  Copyright (c) 2015 Malcolm Ramsay. All rights reserved.
//

#ifndef MY_ANALYSE_H
#define MY_ANALYSE_H

#include <stdio.h>
#include <map>
#include "mol_functions.h"
#include "angle_list.h"
#include "neighbours.h"

int analyse(Frame *, std::vector<Frame *>, int print = 0, int movie = 0, int dist = 0);

#endif /* defined(MY_ANALYSE_H) */
