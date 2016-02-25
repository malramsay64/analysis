//
// Created by malcolm on 21/01/16.
//

#ifndef ANALYSIS_TESTVECTOR_H
#define ANALYSIS_TESTVECTOR_H

#include <iostream>
#include <gtest/gtest.h>
#include "../src/Vector.h"

using namespace LAlgebra;

class VectorTest : public testing::Test {
public:
    static const size_t N = 3;
    const double accuracy{1e-14};
    virtual void SetUp(){
        v0 = Vector{}; // Vector of zeros
        v1 = Vector{}; // Vector of 1s
        v2 = Vector{}; // Vector of 1.1
        v3 = Vector{}; // Vector from 0 to i
        v4 = Vector{}; // Vector from 0 to 2*i
        v5 = Vector{}; // Vector from 0.1 to 2*i+0.1

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

    virtual void expect_near(const Vector &expected, const Vector &actual, double abs_error){
        for (auto i=0; i<N; i++){
            EXPECT_NEAR(expected[i], actual[i], abs_error);
        }
    }

    virtual void expect_double_eq(const Vector &expected, const Vector &actual){
        for (auto i=0; i<N; i++){
            EXPECT_DOUBLE_EQ(expected[i], actual[i]);
        }
    }

    Vector v0, v1, v2, v3, v4, v5;
};


TEST_F(VectorTest, Initialisation){
    EXPECT_EQ(0, this->v0.sum());
    EXPECT_EQ(this->size(), this->v1.r.sum());
    EXPECT_DOUBLE_EQ(this->size()*1.1, this->v2.sum());
    int sum{0};
    for (int i=1;i<=this->size(); i++) sum += i;
    EXPECT_DOUBLE_EQ(sum, this->v3.sum());
    EXPECT_DOUBLE_EQ(sum*2, this->v4.sum());
    EXPECT_DOUBLE_EQ(sum*2+this->size()*0.1, this->v5.sum());
}

TEST_F(VectorTest, Size){
    EXPECT_EQ(this->size(), this->v0.size());
    EXPECT_EQ(this->size(), this->v1.size());
    EXPECT_EQ(this->size(), this->v2.size());
    EXPECT_EQ(this->size(), this->v3.size());
    EXPECT_EQ(this->size(), this->v4.size());
    EXPECT_EQ(this->size(), this->v5.size());
}

TEST_F(VectorTest, SelectOperator){
    EXPECT_EQ(0, this->v0[0]);
    EXPECT_NO_THROW(this->v0[this->size()+1]);
    EXPECT_EQ(0, this->v0[this->size()-1]);
    EXPECT_EQ(this->size(), this->v3[this->size()-1]);
    EXPECT_EQ(1, this->v1[0]);
}

TEST_F(VectorTest, Length){
    EXPECT_DOUBLE_EQ(0, this->v0.length());
    EXPECT_DOUBLE_EQ(sqrt(this->size()), this->v1.length());
    EXPECT_DOUBLE_EQ(sqrt(this->size())*1.1, this->v2.length());
    EXPECT_DOUBLE_EQ(sqrt(this->cum_sum2()), this->v3.length());
    EXPECT_DOUBLE_EQ(2*sqrt(this->cum_sum2()), this->v4.length());
    double sum{0};
    for (int i=1;i<=this->size(); i++) sum += (2*i+0.1)*(2*i+0.1);
    EXPECT_DOUBLE_EQ(sqrt(sum), this->v5.length());
}

TEST_F(VectorTest, Normalise){
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

TEST_F(VectorTest, Orthogonalise){
    Vector vo0 = this->v0.orthogonal();
    EXPECT_DOUBLE_EQ(0, (this->v0*vo0).sum());
    Vector vo1 = this->v1.orthogonal();
    EXPECT_DOUBLE_EQ(0, (this->v1*vo1).sum());
    Vector vo2 = this->v2.orthogonal();
    EXPECT_DOUBLE_EQ(0, (this->v2*vo2).sum());
    Vector vo3 = this->v3.orthogonal();
    EXPECT_DOUBLE_EQ(0, (this->v3*vo3).sum());
    Vector vo4 = this->v4.orthogonal();
    EXPECT_DOUBLE_EQ(0, (this->v4*vo4).sum());
    Vector vo5 = this->v5.orthogonal();
    EXPECT_DOUBLE_EQ(0, (this->v5*vo5).sum());
}

TEST_F(VectorTest, Wrap) {
    Vector vw1 = v1 * 2 * PI + 0.1;
    // Default
    expect_near(v1 * 0.1, wrap(vw1), accuracy);
    expect_near(v1, wrap(v1 + 2 * PI), accuracy);
    expect_near(v3, wrap(v3 - 2 * PI), accuracy);
    expect_near(v1, wrap(v1 + 4 * PI), accuracy);
    expect_near(v1, wrap(v1 + 16 * PI), accuracy);
    // Double limit
    expect_near(v1, wrap(v1,1), accuracy);
    expect_near(v2-1, wrap(v2,1), accuracy);
    // Vector limit
    expect_near(v2-1, wrap(v2,v1), accuracy);

}

TEST_F(VectorTest, Sum){
    // TODO
    EXPECT_EQ(size(), v1.sum());
    EXPECT_EQ(size()*1.1, v2.sum());
    EXPECT_EQ(0, v0.sum());
    EXPECT_EQ(this->cum_sum(), v3.sum());
}

TEST_F(VectorTest, Angular){
    EXPECT_EQ(v0, v0.angular());
    expect_near(Vector{sqrt(3), asin(sqrt(2./3)), PI/4}, v1.angular(), accuracy);
    expect_near(Vector{2,PI/2,0}, Vector{2,0,0}.angular(), accuracy);
    expect_near(Vector{2, PI/2, PI/2}, Vector{0,2,0}.angular(), accuracy);
    expect_near(Vector{2, 0, 0}, Vector{0,0,2}.angular(), accuracy);
}

TEST_F(VectorTest, Distance){
    // TODO distance between two vectors
}

TEST_F(VectorTest, Direction){
    // TODO direction between two vectors
    this->expect_double_eq(this->v1, direction(this->v0, this->v1));
    this->expect_double_eq(-this->v1, direction(this->v1, this->v0));
}

TEST_F(VectorTest, DotProduct){
    // TODO
}

TEST_F(VectorTest, MathFunctions){
    /* TODO
     * sin
     * cos
     * atan2
     */
}

TEST_F(VectorTest, MathOperators){
    /* TODO
     * +, +=
     * -, -=
     * *, *=
     * /, /=
     */
    EXPECT_EQ(this->v0, this->v1-this->v1);
    EXPECT_EQ(this->v0, this->v2-this->v2);
    EXPECT_EQ(this->v0, this->v3-this->v3);
    EXPECT_EQ(this->v0, this->v4-this->v4);
    this->expect_double_eq(this->v1, (this->v2-this->v1)*10);
}

TEST_F(VectorTest, BoolOperators){
    EXPECT_EQ(this->v0, this->v0);
    EXPECT_EQ(Vector(), this->v0);
    EXPECT_EQ(this->v1, Vector(this->v1));
    EXPECT_EQ(this->v5,this->v5);
    EXPECT_NE(this->v4, this->v5);
    EXPECT_NE(this->v1, this->v0);
    EXPECT_NE(this->v1, this->v5);
}

TEST_F(VectorTest, FstreamOperators){
    /* TODO
     *
     * >>
     * <<
     */
    std::stringstream ss;
    ss << this->v0;
}

#endif //ANALYSIS_TESTVECTOR_H
