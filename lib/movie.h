//
//  movie.h
//  analysis
//
//  Created by Malcolm Ramsay on 9/02/2015.
//  Copyright (c) 2015 Malcolm Ramsay. All rights reserved.
//

#ifndef __analysis__movie__
#define __analysis__movie__

#include <stdio.h>
#include "neighbours.h"

double movie_colour(molecule * mol, Frame *);
int print_movie(std::ofstream * file, Frame * frame, molecule * mol = 0);

#endif /* defined(__analysis__movie__) */
