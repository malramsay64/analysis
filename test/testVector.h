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

    virtual void expect_near(const Vector<N> &expected, const Vector<N> &actual, double abs_error){
        for (auto i=0; i<N; i++){
            EXPECT_NEAR(expected[i], actual[i], abs_error);
        }
    }

    virtual void expect_double_eq(const Vector<N> &expected, const Vector<N> &actual){
        for (auto i=0; i<N; i++){
            EXPECT_DOUBLE_EQ(expected[i], actual[i]);
        }
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

TYPED_TEST(VectorTest, Size){
    EXPECT_EQ(this->size(), this->v0.size());
    EXPECT_EQ(this->size(), this->v1.size());
    EXPECT_EQ(this->size(), this->v2.size());
    EXPECT_EQ(this->size(), this->v3.size());
    EXPECT_EQ(this->size(), this->v4.size());
    EXPECT_EQ(this->size(), this->v5.size());
}

TYPED_TEST(VectorTest, SelectOperator){
    EXPECT_EQ(0, this->v0[0]);
    EXPECT_NO_THROW(this->v0[this->size()+1]);
    EXPECT_EQ(0, this->v0[this->size()-1]);
    EXPECT_EQ(this->size(), this->v3[this->size()-1]);
    EXPECT_EQ(1, this->v1[0]);
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

TYPED_TEST(VectorTest, Orthogonalise){
    Vector<this->N> vo0 = this->v0.orthogonal();
    EXPECT_DOUBLE_EQ(0, (this->v0*vo0).sum());
    Vector<this->N> vo1 = this->v1.orthogonal();
    EXPECT_DOUBLE_EQ(0, (this->v1*vo1).sum());
    Vector<this->N> vo2 = this->v2.orthogonal();
    EXPECT_DOUBLE_EQ(0, (this->v2*vo2).sum());
    Vector<this->N> vo3 = this->v3.orthogonal();
    EXPECT_DOUBLE_EQ(0, (this->v3*vo3).sum());
    Vector<this->N> vo4 = this->v4.orthogonal();
    EXPECT_DOUBLE_EQ(0, (this->v4*vo4).sum());
    Vector<this->N> vo5 = this->v5.orthogonal();
    EXPECT_DOUBLE_EQ(0, (this->v5*vo5).sum());
}

TYPED_TEST(VectorTest, Wrap) {
    //TODO test wrap function
    Vector<this->N> v1 = this->v1 * 2 * PI + 0.1;
    this->expect_near(this->v1 * 0.1, wrap(v1), 1e-15);

}

TYPED_TEST(VectorTest, Sum){
    // TODO
}

TYPED_TEST(VectorTest, Angular){
    // TODO
}

TYPED_TEST(VectorTest, Distance){
    // TODO distance between two vectors
}

TYPED_TEST(VectorTest, Direction){
    // TODO direction between two vectors
}

TYPED_TEST(VectorTest, DotProduct){
    // TODO
}

TYPED_TEST(VectorTest, MathFunctions){
    /* TODO
     * sin
     * cos
     * atan2
     */
}

TYPED_TEST(VectorTest, MathOperators){
    /* TODO
     * +, +=
     * -, -=
     * *, *=
     * /, /=
     */
}

TYPED_TEST(VectorTest, BoolOperators){
    EXPECT_EQ(this->v0, this->v0);
    EXPECT_EQ(Vector<this->N>(), this->v0);
    EXPECT_EQ(this->v1, Vector<this->N>(this->v1));
    EXPECT_EQ(this->v5,this->v5);
    EXPECT_NE(this->v4, this->v5);
    EXPECT_NE(this->v1, this->v0);
    EXPECT_NE(this->v1, this->v5);
}

TYPED_TEST(VectorTest, FstreamOperators){
    /* TODO
     */
}

#endif //ANALYSIS_TESTVECTOR_H
