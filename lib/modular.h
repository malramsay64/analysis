//
//  modular.h
//  analysis
//
//  Created by Malcolm Ramsay on 6/02/2015.
//  Copyright (c) 2015 Malcolm Ramsay. All rights reserved.
//

#ifndef MY_MODULAR
#define MY_MODULAR

#include <map>
#include "neighbours.h"
#include "distribution.h"
#include "movie.h"
#include "output.h"


extern bool time_structure;

int mod_analyse(Frame *, std::vector<Frame *>, int regio = 0, int print = 0, int movie = 0, \
                int dist = 0);

#endif /* defined(MY_MODULAR) */