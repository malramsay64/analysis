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
    std::list<Molecule *> q;
    std::list<int> depth;
    std::set<Molecule *> traversed;
public:
    int push(Molecule *, int);
    dyn_queue(Molecule *);
    dyn_queue();
    Molecule * pop_same();
    Molecule * pop();
    int remove(Molecule *);
    int get_depth(); 
    
};

#endif /* defined(DYN_QUEUE_H) */
