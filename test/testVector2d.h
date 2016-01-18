//
// Created by malcolm on 18/01/16.
//

#ifndef ANALYSIS_TESTVECTOR2D_H
#define ANALYSIS_TESTVECTOR2D_H

#include <gtest/gtest.h>
#include "../src/Vector2d.h"

TEST(Vector2d, Initialisation){
    Vector2d v1{};
    EXPECT_DOUBLE_EQ(0, v1.x);
    EXPECT_DOUBLE_EQ(0, v1.y);
    Vector2d v2(0.1,0.2);
    EXPECT_DOUBLE_EQ(0.1, v2.x);
    EXPECT_DOUBLE_EQ(0.2, v2.y);
    double d[2]{0.1,0.3};
    Vector2d v3{d};
    EXPECT_DOUBLE_EQ(0.1, v3.x);
    EXPECT_DOUBLE_EQ(0.3, v3.y);
}
TEST(Vector2d, Length){
    Vector2d v1{};
    EXPECT_DOUBLE_EQ(0, v1.length());
    Vector2d v2{1,0};
    EXPECT_DOUBLE_EQ(1, v2.length());
    Vector2d v3{0, 0.1};
    EXPECT_DOUBLE_EQ(0.1, v3.length());
    Vector2d v4{1,1};
    EXPECT_DOUBLE_EQ(sqrt(2), v4.length());
    Vector2d v5{0.3,0.3};
    EXPECT_DOUBLE_EQ(sqrt(2*0.3*0.3), v5.length());
}

TEST(Vector2d, Normalise){
    Vector2d v1{};
    EXPECT_DOUBLE_EQ(0,v1.x);
    EXPECT_DOUBLE_EQ(0,v1.y);
    Vector2d v2{0.1,0};
    v2.normalise();
    EXPECT_DOUBLE_EQ(1,v2.x);
    EXPECT_DOUBLE_EQ(0,v2.y);
    Vector2d v3{0, 0.1};
    v3.normalise();
    EXPECT_DOUBLE_EQ(0,v3.x);
    EXPECT_DOUBLE_EQ(1,v3.y);
    Vector2d v4{1,1};
    v4.normalise();
    EXPECT_DOUBLE_EQ(1/sqrt(2),v4.x);
    EXPECT_DOUBLE_EQ(1/sqrt(2),v4.y);
}

TEST(Vector2d, Orthogonalise){
    Vector2d v1{};
    EXPECT_DOUBLE_EQ(0,v1.x);
    EXPECT_DOUBLE_EQ(0,v1.y);
    Vector2d v2{1,0};
    v2.orthogonalise();
    EXPECT_DOUBLE_EQ(0,v2.x);
    EXPECT_DOUBLE_EQ(1,v2.y);
    Vector2d v3{0.4,0.5};
    v3.orthogonalise();
    EXPECT_DOUBLE_EQ(-0.5, v3.x);
    EXPECT_DOUBLE_EQ(0.4, v3.y);
}

TEST(Vector2d, Angle){
    Vector2d v1{};
    EXPECT_DOUBLE_EQ(0, v1.angle());
    Vector2d v2{1,0};
    EXPECT_DOUBLE_EQ(0, v2.angle());
    Vector2d v3{0,1};
    EXPECT_DOUBLE_EQ(PI/2, v3.angle());
    Vector2d v4{-1,0};
    EXPECT_DOUBLE_EQ(PI, v4.angle());
    Vector2d v5{0,-1};
    EXPECT_DOUBLE_EQ(-PI/2, v5.angle());
    Vector2d v6{1,1};
    EXPECT_DOUBLE_EQ(PI/4, v6.angle());
}

TEST(Vector2d, Summation){

}

#endif //ANALYSIS_TESTVECTOR2D_H
