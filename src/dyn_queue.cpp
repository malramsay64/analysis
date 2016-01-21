//
//  dyn_queue.cpp
//  analysis
//
//  Created by Malcolm Ramsay on 7/12/2014.
//  Copyright (c) 2014 Malcolm Ramsay. All rights reserved.
//

#include "dyn_queue.h"

using namespace std;

dyn_queue::dyn_queue(){
}

dyn_queue::dyn_queue(Molecule * m){
    push(m, 0);
}

int dyn_queue::push(Molecule * t, int d){
    if (!traversed.count(t)){
        traversed.insert(t);
        q.push_back(pair<Molecule *, int>{t, d});
    }
    return 0;
}

Molecule * dyn_queue::pop(){
    if (q.size()){
        Molecule *t = q.front().first;
        int depth = q.front().second;
        for (auto &i: t->my_neighbours){
            push(i.first, depth+1);
        }
        q.pop_front();
        return t;
    }
    return nullptr;
}

int dyn_queue::get_depth(){
    return q.front().second;
}
