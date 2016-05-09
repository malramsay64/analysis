//
// Created by malcolm on 29/01/16.
//

#ifndef ANALYSIS_TESTCOMPUTE_H
#define ANALYSIS_TESTCOMPUTE_H

#include <gtest/gtest.h>
#include "../src/Compute.h"

/* TODO
 * Create computeTest class
 *  - This class will set up a frame with a number of particles and molecules
 *  - It will be repsonsible for the construction and destruction
 *  - It will have functions for testing the resulting output of the test to a 'known' value.
 *
 *  This will possibly be a subclass of the frame testing class as many of the features are going to be the same.
 *
 *  This class will then test all the compute styles, checking they get sensible results.
 */

class ComputeTest : public testing::Test {
protected:
    Frame frame1;
    Frame frame2;
    int natoms{10};
    virtual void SetUp(){
        frame1.set_atoms(natoms);
        frame2.set_atoms(natoms);
        for (int i=0; i<natoms; i++){
            frame1.particles[i]=Particle{};
            frame2.particles[i]=Particle{};
            frame1.particles[i].set_pos(Vector{});
            frame2.particles[i].set_pos(Vector{1,1,1});
        }
    }

    virtual void TearDown() {};

    void EXPECT_VECTOR_DOUBLE_EQ(std::vector<double> &expected, std::vector<double> &actual){
        ASSERT_EQ(expected.size(), actual.size());
        for (int i=0; i<expected.size(); i++){
            EXPECT_DOUBLE_EQ(expected[i], actual[i]);
        }
    }

};


TEST(Compute, ComputeAngle){
    /* TODO
     */
}

TEST(Compute, ComputeCoord){
    /* TODO
     */
}

TEST(Compute, ComputeMolHexatic){
    /* TODO
     */
}

TEST(Compute, ComputeCircleOrder){
    /* TODO
     */
}

TEST(Compute, ComputeOrientationOrder){
    /* TODO
     */
}

TEST_F(ComputeTest, ComputeDisplacement){
    ComputeDisplacement c1(frame1, frame1);
    ComputeDisplacement c2(frame1, frame2);
    c1.compute_array();
    c2.compute_array();
    std::vector<double> c1_expect(natoms,0.);
    std::vector<double> c2_expect(natoms, std::sqrt(3));
    EXPECT_VECTOR_DOUBLE_EQ(c1_expect, c1.array);
    EXPECT_VECTOR_DOUBLE_EQ(c2_expect, c2.array);
}

#endif //ANALYSIS_TESTCOMPUTE_H
