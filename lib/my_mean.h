//
//  my_mean.h
//  analysis
//
//  Created by Malcolm Ramsay on 7/12/2014.
//  Copyright (c) 2014 Malcolm Ramsay. All rights reserved.
//

#ifndef MY_MEAN
#define MY_MEAN

class my_mean{
    int n;
    double m;
    double stdev;
public:
    my_mean();
    double add(double);
    double get_mean();
    double get_stdev();
    double combine(my_mean);
    
};

#endif /* defined(MY_MEAN) */
