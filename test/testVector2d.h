//
// Created by malcolm on 18/01/16.
//

#ifndef ANALYSIS_TESTVECTOR2D_H
#define ANALYSIS_TESTVECTOR2D_H

#include <gtest/gtest.h>
#include "../src/Vector2d.h"

using namespace LAlgebra;

TEST(Vector2d, Initialisation){
    Vector2d v1{};
    EXPECT_DOUBLE_EQ(0, v1.x);
    EXPECT_DOUBLE_EQ(0, v1.y);
    Vector2d v2(0.1, 0.2);
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
    Vector2d v2{1, 0};
    EXPECT_DOUBLE_EQ(1, v2.length());
    Vector2d v3{0, 0.1};
    EXPECT_DOUBLE_EQ(0.1, v3.length());
    Vector2d v4{1, 1};
    EXPECT_DOUBLE_EQ(sqrt(2), v4.length());
    Vector2d v5{0.3, 0.3};
    EXPECT_DOUBLE_EQ(sqrt(2*0.3*0.3), v5.length());
}

TEST(Vector2d, Normalise){
    Vector2d v1{};
    EXPECT_DOUBLE_EQ(0,v1.x);
    EXPECT_DOUBLE_EQ(0,v1.y);
    Vector2d v2{0.1, 0};
    v2.normalise();
    EXPECT_DOUBLE_EQ(1,v2.x);
    EXPECT_DOUBLE_EQ(0,v2.y);
    Vector2d v3{0, 0.1};
    v3.normalise();
    EXPECT_DOUBLE_EQ(0,v3.x);
    EXPECT_DOUBLE_EQ(1,v3.y);
    Vector2d v4{1, 1};
    v4.normalise();
    EXPECT_DOUBLE_EQ(1/sqrt(2),v4.x);
    EXPECT_DOUBLE_EQ(1/sqrt(2),v4.y);
}

TEST(Vector2d, Orthogonalise){
    Vector2d v1{};
    EXPECT_DOUBLE_EQ(0,v1.x);
    EXPECT_DOUBLE_EQ(0,v1.y);
    Vector2d v2{1, 0};
    v2.orthogonalise();
    EXPECT_DOUBLE_EQ(0,v2.x);
    EXPECT_DOUBLE_EQ(1,v2.y);
    Vector2d v3{0.4, 0.5};
    v3.orthogonalise();
    EXPECT_DOUBLE_EQ(-0.5, v3.x);
    EXPECT_DOUBLE_EQ(0.4, v3.y);
}

TEST(Vector2d, Angle){
    Vector2d v1{};
    EXPECT_DOUBLE_EQ(0, v1.angle());
    Vector2d v2{1, 0};
    EXPECT_DOUBLE_EQ(0, v2.angle());
    Vector2d v3{0, 1};
    EXPECT_DOUBLE_EQ(PI/2, v3.angle());
    Vector2d v4{-1, 0};
    EXPECT_DOUBLE_EQ(PI, v4.angle());
    Vector2d v5{0, -1};
    EXPECT_DOUBLE_EQ(-PI/2, v5.angle());
    Vector2d v6{1, 1};
    EXPECT_DOUBLE_EQ(PI/4, v6.angle());
}

TEST(Vector2d, Summation){
    Vector2d v1{}, v2{1,1}, v3{1.1,1.1}, v4{2.5,1.5}, v5{1.5,0.5};
    EXPECT_EQ(v2, v1+1);
    EXPECT_EQ(v3, v1+1.1);
    EXPECT_EQ(v3, v2+0.1);
    EXPECT_EQ(v4, v5+1);
    v5 += 1;
    v1 += 1.1;
    EXPECT_EQ(v3, v1);
    EXPECT_EQ(Vector2d(), v1-1.1);
    EXPECT_EQ(Vector2d(), v1 -= 1.1);
    EXPECT_EQ(v5-=1, v4-1);
    EXPECT_EQ(v5, v4-1);
    EXPECT_EQ(v5, v4-v2);
    EXPECT_EQ(v4, v5+v2);
    v5 += v2;
    EXPECT_EQ(v4, v5);
    v5 -= v3;
    v4 -= v3;
    EXPECT_EQ(v5,v4);
    v3 -= v3;
    EXPECT_EQ(v1, v3);
}

TEST(Vector2d, Multiplication){
    Vector2d v1{}, v2{1,1}, v3{1.1,1.1}, v4{4.5,1.5}, v5{1.5,0.5};
    EXPECT_EQ(v3, v2*1.1);
    EXPECT_EQ(v4, v5*3);
    EXPECT_EQ(v1, v2*0);
    EXPECT_EQ(v1,v3*v1);
    EXPECT_EQ(v1,v5*v1);
    EXPECT_EQ(v3, v2*v3);
    EXPECT_EQ(Vector2d(1.5*4.5,0.5*1.5), v4*v5);
    EXPECT_EQ(v2, v3/1.1);
    EXPECT_EQ(v5, v4/3);
    EXPECT_EQ(v2/1.1, v2/v3);
    v1 += 1;
    EXPECT_EQ(v2,v1);
    v1 *= 1.1;
    EXPECT_EQ(v3,v1);
    v1 /= 1.1;
    v1 += Vector2d{0.5,-0.5};
    EXPECT_EQ(v1, v5);
    v1 *= 3;
    EXPECT_EQ(v4, v1);
    v1 /= 3;
    EXPECT_EQ(v5,v1);
}

TEST(Vector2d, Equality){
    Vector2d v1{}, v2{0, 1}, v3{1,0}, v4{1,1}, v5{1,1};
    EXPECT_TRUE(v1 == v1);
    EXPECT_TRUE(v1 != v2);
    EXPECT_TRUE(v1 != v3);
    EXPECT_TRUE(v1 != v4);
    EXPECT_FALSE(v1 == v5);
    EXPECT_FALSE(v4 != v5);
}

TEST(Vector2d, Direction){
    Vector2d v1{}, v2{0,1}, v3{1,0}, v4{1,1};
    EXPECT_EQ(v1, direction(v1,v1));
    EXPECT_EQ(v1, direction(v4,v4));
    EXPECT_EQ(v2, direction(v1,v2));
    EXPECT_EQ(-v2, direction(v2,v1));
    EXPECT_EQ(Vector2d(-1,1), direction(v3,v2));
    EXPECT_EQ(v2,direction(v3,v4));
    EXPECT_EQ(1, dist(v1,v2));
    EXPECT_EQ(0, dist(v2,v2));
    EXPECT_EQ(sqrt(2), dist(v1,v4));
    EXPECT_EQ(sqrt(2), dist(v4,v1));

    // Trying to break it!!
    Vector2d v5{2*PI, 2*PI}, v6{PI, PI}, v7{PI, 0}, v8{0,PI}, v9{0.25,0.45}, v10{2*PI-0.5, 1.2};
    EXPECT_NEAR(0, direction(v5,v1).x, 1e-15);
    EXPECT_NEAR(0, direction(v5,v1).y, 1e-15);
    EXPECT_NEAR(-PI+1, direction(v6,v4).x, 1e-15);
    EXPECT_NEAR(-PI+1, direction(v6,v4).y, 1e-15);
    EXPECT_NEAR(1, direction(v5,v3).x, 1e-15);
    EXPECT_NEAR(0, direction(v5,v3).y, 1e-15);
    EXPECT_NEAR(1, direction(v5,v3).x, 1e-15);
    EXPECT_NEAR(0, direction(v5,v3).y, 1e-15);
    EXPECT_NEAR(-0.75, direction(v9,v10).x, 1e-15);
    EXPECT_NEAR(0.75, direction(v9,v10).y, 1e-15);
    EXPECT_NEAR(0.75, direction(v10,v9).x, 1e-15);
    EXPECT_NEAR(-0.75, direction(v10,v9).y, 1e-15);

}

TEST(Vector2d, Math){
    Vector2d v0{}, v1{1,0}, v2{0,1}, v3{1,1}, v4{2,1};
    EXPECT_EQ(0,dot_product(v0,v0));
    EXPECT_EQ(1, dot_product(v1,v1));
    EXPECT_EQ(0, dot_product(v1,v2));
    EXPECT_EQ(2, dot_product(v3,v3));
    EXPECT_EQ(3, dot_product(v4,v3));
    EXPECT_EQ(-1, dot_product(-v4, v2));
    EXPECT_NEAR(1, sin(v1*PI/2).x, 1e-15);
    EXPECT_NEAR(0, sin(v1*PI/2).y, 1e-15);
    EXPECT_NEAR(0, cos(v1*PI/2).x, 1e-15);
    EXPECT_NEAR(1, cos(v1*PI/2).y,1e-15);
    EXPECT_NEAR(0, sin(v4*PI/2).x,1e-15);
    EXPECT_NEAR(1, sin(v4*PI/2).y,1e-15);
    EXPECT_NEAR(-1, cos(v4*PI/2).x,1e-15);
    EXPECT_NEAR(0, cos(v4*PI/2).y,1e-15);
    EXPECT_NEAR(0, sin(v4*PI).x,1e-15);
    EXPECT_NEAR(0, sin(v4*PI).y,1e-15);
    EXPECT_NEAR(1, cos(v4*PI).x,1e-15);
    EXPECT_NEAR(-1, cos(v4*PI).y,1e-15);
    EXPECT_NEAR(0, tan(v4*v4*PI/4).x,1e-15);
    EXPECT_NEAR(1, tan(v4*v4*PI/4).y,1e-15);
    EXPECT_DOUBLE_EQ(0, atan2(v0));
    EXPECT_DOUBLE_EQ(0, atan2(v1));
    EXPECT_DOUBLE_EQ(-PI, atan2(-v1));
    EXPECT_DOUBLE_EQ(PI/2, atan2(v2));
    EXPECT_DOUBLE_EQ(-PI/2, atan2(-v2));
    EXPECT_DOUBLE_EQ(PI/4, atan2(v3));
    EXPECT_DOUBLE_EQ(PI/3, atan2(Vector2d(1,sqrt(3))));

    /* TODO
     * map_angle
     */
}

TEST(Vector2d, Output){
    /* TODO
     * operator<<
     */
}

#endif //ANALYSIS_TESTVECTOR2D_H
