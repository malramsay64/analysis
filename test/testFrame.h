//
// Created by malcolm on 20/01/16.
//

#ifndef ANALYSIS_TESTFRAME_H
#define ANALYSIS_TESTFRAME_H

#include <gtest/gtest.h>
#include "../src/Frame.h"

TEST(Frame, Initialisation){
    /* TODO
     * Check all values
     * Check vectors initialised
     */
}

TEST(Frame, StructInitialiser){
    /* TODO
     * Test initialisation with struct
     */
}

TEST(Frame, Getters){
    /* TODO
     *
     * get_timestep()
     * get_stepsize()
     * get_time()
     * get_area()
     * get_density()
     * get_a()
     * get_b()
     * get_theta()
     * get_height()
     * get_tilt()
     */
}

TEST(Frame, Linking){
    /* TODO
     *
     * add_link(int, int)
     * add_link(int, int, bool)
     * update_links()
     */
}

TEST(Frame, Conversions){
    /* TODO
     *
     * cartesian
     * fractional
     * dist
     * direction
     * angle
     */
}

#endif //ANALYSIS_TESTFRAME_H
