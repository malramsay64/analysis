//
//  mol_functions.h
//  analysis
//
//  Created by Malcolm Ramsay on 7/12/2014.
//  Copyright (c) 2014 Malcolm Ramsay. All rights reserved.
//

#ifndef MY_MOL_FUNCTIONS
#define MY_MOL_FUNCTIONS

#include "particle.h"
#include "frame.h"
#include "my_mean.h"
#include "dyn_queue.h"

vect orientation(molecule *, Frame *);
double angle(molecule *, Frame *);
std::ostream& print_mol(std::ostream &, molecule *,Frame *);


#endif /* defined(MY_MOL_FUNCTIONS) */

