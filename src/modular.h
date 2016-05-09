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
#include <iomanip>
#include "neighbours.h"
//#include "distribution.h"
#include "movie.h"
#include "output.h"
#include "functions.h"


extern bool time_structure, regio, movie, moved, m_orient;

int mod_analyse(Frame *, std::vector<Frame *>, int print = 0, int dist = 0);

#endif /* defined(MY_MODULAR) */
