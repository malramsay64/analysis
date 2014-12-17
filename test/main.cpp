//
//  main.cpp
//  test
//
//  Created by Malcolm Ramsay on 18/12/2014.
//  Copyright (c) 2014 Malcolm Ramsay. All rights reserved.
//

#include <iostream>
#include <cmath>
#include "functions.h"

using namespace std;

//static double PI = 3.141526535;

int main(int argc, const char * argv[]) {
    // insert code here...
    std::cout << "Hello, World!\n";
    cout << my_mod(-120*PI/180,PI) << endl;
    cout << atan2(sin(-2*PI/3), cos(-2*PI/3)) << endl;
    cout << acos(cos(-2*PI/3)) << endl;
    cout << 3 % 4 << endl;
    cout << -3 % 4 << endl;
    cout << 4 % 3 << endl;
    cout << (-4) % 3 << endl;
    return 0;
}
