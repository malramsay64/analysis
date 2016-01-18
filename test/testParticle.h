//
// Created by malcolm on 18/01/16.
//

#ifndef ANALYSIS_TESTPARTICLE_H
#define ANALYSIS_TESTPARTICLE_H

#include <gtest/gtest.h>
#include "../src/Particle.h"

TEST(Particle, DefaultInitialiser){
    Particle p{};
    EXPECT_EQ(0, p.pos.x);
    EXPECT_EQ(0, p.pos.y);
    EXPECT_DOUBLE_EQ(1, p.mass);
    EXPECT_EQ(0, p.molid);
    EXPECT_DOUBLE_EQ(0, p.radius);
    EXPECT_EQ(0, p.type);
    EXPECT_EQ(0, p.my_neighbours.size());
}

Particle default_p(){
    Particle p{};
    p.mass = 0.8;
    p.id = 3;
    p.molid = 5;
    p.type = 2;
    p.radius = 1.2;
    p.pos = Vector2d{0.1,0.2};
    p.my_neighbours = std::vector<Particle *>{&p};
    return p;
}


TEST(Particle, Initialising){
    Particle p = default_p();
    EXPECT_DOUBLE_EQ(0.8, p.mass);
    EXPECT_EQ(3, p.id);
    EXPECT_EQ(5, p.molid);
    EXPECT_EQ(2, p.type);
    EXPECT_DOUBLE_EQ(1.2, p.radius);
    EXPECT_TRUE(p.pos == Vector2d(0.1, 0.2));
}

TEST(Particle, CopyConstructor){
    Particle p1 = default_p();
    Particle p2{p1};
    p1.mass = 1;
    p1.pos.x = 3;
    p1.radius = 0.5;
    EXPECT_DOUBLE_EQ(1, p1.mass);
    EXPECT_DOUBLE_EQ(0.8, p2.mass);
    EXPECT_TRUE(p1.pos == Vector2d(3,0.2));
    EXPECT_TRUE(p2.pos == Vector2d(0.1,0.2));
    EXPECT_DOUBLE_EQ(0.5, p1.radius);
    EXPECT_DOUBLE_EQ(1.2, p2.radius);
}

TEST(Particle, StructConstructor){
    Particle p{particle_vars{}};
    EXPECT_DOUBLE_EQ(1, p.mass);
    EXPECT_EQ(0, p.id);
    EXPECT_EQ(0, p.molid);
    EXPECT_EQ(0, p.type);
    EXPECT_DOUBLE_EQ(1, p.radius);
    EXPECT_TRUE(p.pos == Vector2d(0, 0));
}

TEST(Particle, OperatorsBoolean){
    Particle p1{}, p2{}, p3{};
    p1.id = 1;
    p2.id = 2;
    p3.id = 1;
    EXPECT_TRUE(p1 == p1);
    EXPECT_FALSE(p1 == p2);
    EXPECT_TRUE(p1 == p3);
    EXPECT_TRUE(p1 != p2);
    EXPECT_FALSE(p1 != p3);
    EXPECT_TRUE(p1 <= p2);
    EXPECT_TRUE(p1 <= p3);
    EXPECT_TRUE(p1 < p2);
    EXPECT_FALSE(p2 < p1);
    EXPECT_FALSE(p1 < p3);
    EXPECT_TRUE(p2 >= p1);
    EXPECT_TRUE(p3 >= p1);
    EXPECT_TRUE(p2 > p1);
    EXPECT_FALSE(p3 > p1);
}

#endif //ANALYSIS_TESTPARTICLE_H
