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

template <class type> class dyn_queue{
    std::list<type *> q;
    std::list<int> depth;
    int push(type *, int);
public:
    dyn_queue(type *);
    type * pop();
    int remove(type *);
    int get_depth();
};

template <class type>
dyn_queue<type>::dyn_queue(type * t){
    push(t,0);
}

template <class type>
int dyn_queue<type>::push(type * t, int d){
    if (!t->get_traversed()){
        t->traverse();
        q.push_back(t);
        depth.push_back(d);
    }
    return 0;
}

template <class type>
type * dyn_queue<type>::pop(){
    if (q.size()){
        type * t = q.front();
        q.pop_front();
        typename std::vector<type *>::iterator i;
        for (i = t->my_neighbours.begin(); i != t->my_neighbours.end(); i++){
            if ((*i)->type == t->type){
                push(*i, depth.front()+1);
                
            }
        }
        depth.pop_front();
        return t;
    }
    return 0;
}

template <class type>
int dyn_queue<type>::remove(type * t){
    q.remove(t);
    return 0;
}

template <class type>
int dyn_queue<type>::get_depth(){
    return depth.front();
}

#endif /* defined(DYN_QUEUE_H) */
