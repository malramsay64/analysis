//
//  dyn_queue.h
//  analysis
//
//  Created by Malcolm Ramsay on 7/12/2014.
//  Copyright (c) 2014 Malcolm Ramsay. All rights reserved.
//

#ifndef DYN_QUEUE_H
#define DYN_QUEUE_H

#include <list>
#include <vector>
#include <set>
#include "Molecule.h"

class dyn_queue{
    std::list<molecule *> q;
    std::list<int> depth;
    std::set<molecule *> traversed;
public:
    int push(molecule *, int);
    dyn_queue(molecule *);
    dyn_queue();
    molecule * pop_same();
    molecule * pop();
    int remove(molecule *);
    int get_depth(); 
    
};

#endif /* defined(DYN_QUEUE_H) */
