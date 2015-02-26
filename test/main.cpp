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
#include "frame.h"
#include "angle_list.h"
#include <assert.h>
#include "output.h"

using namespace std;

//static double PI = 3.141526535;

int main(int argc, const char * argv[]) {
    // insert code here...
    //std::cout << "Hello, World!\n";
    //cout << my_mod(-120*PI/180,PI) << endl;
    //cout << atan2(sin(-2*PI/3), cos(-2*PI/3)) << endl;
    //cout << acos(cos(-2*PI/3)) << endl;
    cout << 3 % 4 << endl;
    cout << -3 % 4 << endl;
    cout << 4 % 3 << endl;
    cout << (-4) % 3 << endl;
    
    // Testing Coordinates
    /*
    Frame frame;
    frame.set_crys(2.5, 1.8, 2*PI/3);
    vect v = vect(1.1,0.1);
    vect v1 = frame.cartesian(v);
    cout << "Initial " << v << endl;
    cout << "Fractional " << v1 << endl;
    cout << "Final " << frame.fractional(v1) << endl;
     */
    // Testing Angle List
    /*
    angle_list a;
    a.push(0, 1);
    a.push(0.0005, 4);
    a.push(2*PI-0.0005, -4);
    a.push(PI, 1);
    ostream * file = &cout;
    a.print(file);
    */
    // Testing vect functions
    
    /*
    // Testing my_mean
    vector<int> sample = { 2, 4, 4, 4, 5, 5, 7, 9 };
    my_mean mean_test;
    for (auto &i: sample){
        mean_test.add(i);
    }
    cout << "Mean - Test: " <<  mean_test.get_mean() << " Real: 5" << endl;
    cout << "Stdev - Test: " <<  mean_test.get_stdev() << " Real: 2" << endl;
    */
    cout << " relax " << print_relax_time(0) << endl;
    cout << "NAN " <<  NAN << endl;
    return 0;
}
