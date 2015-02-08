//
//  main.cpp
//  template
//
//  Created by Malcolm Ramsay on 8/02/2015.
//  Copyright (c) 2015 Malcolm Ramsay. All rights reserved.
//

#include <iostream>
#include "my_list.h"

using namespace std;

int main(int argc, const char * argv[]) {
    my_list<int> my_int;
    my_int.add(1);
    my_int.add(2);
    my_int.add(3);
    
    cout << my_int.get_sum() << endl;
    
    my_list<double> my_double;
    my_double.add(1.1);
    my_double.add(2.2);
    my_double.add(3.3);
                  
    cout << my_int.get_sum() << endl;
    
    return 0;
}
