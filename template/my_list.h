//
//  my_list.h
//  analysis
//
//  Created by Malcolm Ramsay on 8/02/2015.
//  Copyright (c) 2015 Malcolm Ramsay. All rights reserved.
//

#ifndef __analysis__my_list__
#define __analysis__my_list__

#include <stdio.h>
#include <vector>

template <class type>
class my_list {
    std::vector<type> list;
    type sum;
    
public:
    my_list();
    
    type get_sum();
    void add(type);
};

template <class type>
my_list<type>::my_list(){
    sum = 0;
}

template <class type>
void my_list<type>::add(type d) {
    list.push_back(d);
    sum += d;
}

template <class type>
type my_list<type>::get_sum() {
    return sum;
}

#endif /* defined(__analysis__my_list__) */
