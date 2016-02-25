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

TEST(Cmpute, ComputeOrientationOrder){
    /* TODO
     */
}

#endif //ANALYSIS_TESTCOMPUTE_H
