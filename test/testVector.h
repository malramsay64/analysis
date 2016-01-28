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

template <size_t N>
class TypeValue{
public:
    static const size_t value = N;
};
template <size_t N> const size_t TypeValue<N>::value;



template <typename T>
class VectorTest : public testing::Test {
public:
    static const size_t N = T::value;
    virtual void SetUp(){
        v0 = Vector<N>{}; // Vector of zeros
        v1 = Vector<N>{}; // Vector of 1s
        v2 = Vector<N>{}; // Vector of 1.1
        v3 = Vector<N>{}; // Vector from 0 to i
        v4 = Vector<N>{}; // Vector from 0 to 2*i
        v5 = Vector<N>{}; // Vector from 0.1 to 2*i+0.1

        v1.r = std::valarray<double>(1.,N);
        v2.r = std::valarray<double>(1.1,N);
        for (int a=1; a<=N; a++){
            v3.r[a-1] = a;
            v4.r[a-1] = 2*a;
            v5.r[a-1] = 2*a+0.1;
        }
    }

    virtual size_t size(){
        return N;
    }

    virtual int cum_sum(){
        int sum{0};
        for (int i=1; i<=N;i++){
            sum += i;
        }
        return sum;
    }

    virtual int cum_sum2(){
        int sum{0};
        for (int i=1; i<=N;i++){
            sum += i*i;
        }
        return sum;
    }

    Vector<N> v0, v1, v2, v3, v4, v5;
};

using testing::Types;

typedef Types<TypeValue<1>, TypeValue<2>,TypeValue<3>, TypeValue<4>> VectorSizes;

TYPED_TEST_CASE(VectorTest, VectorSizes);

TYPED_TEST(VectorTest, Initialisation){
    EXPECT_EQ(0, this->v0.sum());
    EXPECT_EQ(this->size(), this->v1.r.sum());
    EXPECT_DOUBLE_EQ(this->size()*1.1, this->v2.sum());
    int sum{0};
    for (int i=1;i<=this->size(); i++) sum += i;
    EXPECT_DOUBLE_EQ(sum, this->v3.sum());
    EXPECT_DOUBLE_EQ(sum*2, this->v4.sum());
    EXPECT_DOUBLE_EQ(sum*2+this->size()*0.1, this->v5.sum());
}

TYPED_TEST(VectorTest, Length){
    EXPECT_DOUBLE_EQ(0, this->v0.length());
    EXPECT_DOUBLE_EQ(sqrt(this->size()), this->v1.length());
    EXPECT_DOUBLE_EQ(sqrt(this->size())*1.1, this->v2.length());
    EXPECT_DOUBLE_EQ(sqrt(this->cum_sum2()), this->v3.length());
    EXPECT_DOUBLE_EQ(2*sqrt(this->cum_sum2()), this->v4.length());
    double sum{0};
    for (int i=1;i<=this->size(); i++) sum += (2*i+0.1)*(2*i+0.1);
    EXPECT_DOUBLE_EQ(sqrt(sum), this->v5.length());
}

TYPED_TEST(VectorTest, Normalise){
    this->v0.normalise();
    this->v1.normalise();
    this->v2.normalise();
    this->v3.normalise();
    this->v4.normalise();
    this->v5.normalise();
    EXPECT_DOUBLE_EQ(0, this->v0.length());
    EXPECT_DOUBLE_EQ(1, this->v1.length());
    EXPECT_DOUBLE_EQ(1, this->v2.length());
    EXPECT_DOUBLE_EQ(1, this->v3.length());
    EXPECT_DOUBLE_EQ(1, this->v4.length());
    EXPECT_DOUBLE_EQ(1, this->v5.length());
}

TYPED_TEST(VectorTest, ){

}

#endif //ANALYSIS_TESTVECTOR_H
