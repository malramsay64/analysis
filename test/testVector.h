//
// Created by malcolm on 21/01/16.
//

#ifndef ANALYSIS_TESTVECTOR_H
#define ANALYSIS_TESTVECTOR_H

#include <gtest/gtest.h>
#include "../src/Vector2d.h"
#include "testVector2d.h"

using namespace LAlgebra;

//typedef Vector* CreateVect();

template <int i>
Vector<i> * CreateV(){
    return new Vector<i>{};
}

class VectorTest : public testing::TestWithParam<Vector*> {
public:
    virtual void SetUp(){
        v0 = (*GetParam())();
        v1 = (*GetParam())();
        v2 = (*GetParam())();
        v3 = (*GetParam())();
        v4 = (*GetParam())();
        v5 = (*GetParam())();

        for (auto a=0; a<v3->size(); a++){
            v3->r[a] = a;
        }
    }

    Vector *v0, *v1, *v2, *v3, *v4, *v5;

};
TEST_P(VectorTest, Initialisation){
    EXPECT_EQ(0, v0->sum());
}

INSTANTIATE_TEST_CASE_P(
        Vector,
        VectorTest,
        testing::Values(CreateV<2>, CreateV<3>));

#endif //ANALYSIS_TESTVECTOR_H
