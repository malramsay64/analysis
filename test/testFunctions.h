//
// Created by malcolm on 18/01/16.
//

#ifndef ANALYSIS_TESTFUNCTIONS_H
#define ANALYSIS_TESTFUNCTIONS_H

#include <gtest/gtest.h>
#include "../src/functions.h"

TEST(Functions, Defs){
    EXPECT_TRUE(XX == 0);
    EXPECT_TRUE(YY == 1);
    EXPECT_DOUBLE_EQ(M_PI, PI);
}

TEST(Functions, DotProduct){
    double x1[5]{0.1, 0.1, 0.1, 0.1, 0.1}, x2[5]{2.,3.,4.,5.,6.};
    EXPECT_DOUBLE_EQ(2, dot_product(x1,x2,5));
}

TEST(Functions, Legendre){
    EXPECT_DOUBLE_EQ(0, legendre(0, 0));
    EXPECT_DOUBLE_EQ(5, legendre(1,5));
    EXPECT_DOUBLE_EQ(0.1,legendre(1,0.1));
    EXPECT_DOUBLE_EQ(2.3, legendre(1,2.3));
    EXPECT_DOUBLE_EQ(-1, legendre(2,0));
    EXPECT_NEAR(0, legendre(2,1/sqrt(2)), 1e-15);
    EXPECT_DOUBLE_EQ(7, legendre(2,2));
}

TEST(Functions, MyMod){
    EXPECT_DOUBLE_EQ(0,my_mod(0,5));
    EXPECT_DOUBLE_EQ(1,my_mod(1,5));
    EXPECT_DOUBLE_EQ(0,my_mod(5,5));
    EXPECT_DOUBLE_EQ(3,my_mod(8,5));
    EXPECT_DOUBLE_EQ(4,my_mod(-1,5));
    EXPECT_DOUBLE_EQ(1,my_mod(-4,5));
    EXPECT_DOUBLE_EQ(2,my_mod(-8,5));
    EXPECT_DOUBLE_EQ(3,my_mod(-7,5));
}

TEST(Functions, PosDefMod){
    EXPECT_DOUBLE_EQ(0,pos_def_mod(0,5));
    EXPECT_DOUBLE_EQ(1,pos_def_mod(1,5));
    EXPECT_DOUBLE_EQ(0,pos_def_mod(5,5));
    EXPECT_DOUBLE_EQ(3,pos_def_mod(8,5));
    EXPECT_DOUBLE_EQ(4,pos_def_mod(-1,5));
    EXPECT_DOUBLE_EQ(1,pos_def_mod(-4,5));
    EXPECT_DOUBLE_EQ(2,pos_def_mod(-8,5));
    EXPECT_DOUBLE_EQ(3,pos_def_mod(-7,5));
}

TEST(Functions, Dist){
    EXPECT_DOUBLE_EQ(0, dist(0,0));
    EXPECT_DOUBLE_EQ(1, dist(0,1));
    EXPECT_DOUBLE_EQ(PI, dist(0,PI));
    EXPECT_DOUBLE_EQ(PI-0.1, dist(-0.1, PI));
    EXPECT_DOUBLE_EQ(PI-0.5, dist(0.5, PI+1));
    EXPECT_DOUBLE_EQ(1.5, dist(2*PI-1, 0.5));
}

#endif //ANALYSIS_TESTFUNCTIONS_H
