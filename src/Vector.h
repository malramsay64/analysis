//
// Created by malcolm on 20/01/16.
//

#ifndef ANALYSIS_VECTOR_H
#define ANALYSIS_VECTOR_H

#include <valarray>
#include "functions.h"

namespace LAlgebra {

    class Vector {
    private:

    public:
        std::valarray<double> r;

        Vector() : r(std::valarray<double>(0.,3)) {};
        Vector(const std::valarray<double> &d) : r(d) {};
        Vector(std::initializer_list<double> d): r(std::valarray<double>(d)) {};
        Vector(const Vector &v) : r(std::valarray<double>{v.r}) {};

        double length() const;
        void normalise();
        Vector orthogonal() const;
        size_t size() const;
        double sum() const;
        Vector angular() const;

        double& operator[] (size_t s){ return r[s]; }
        double operator[] (size_t s) const{ return r[s]; }
    };


    double dist(const Vector &from, const Vector &to);
    Vector direction(const Vector &v1, const Vector &v2);
    double dot_product(const Vector &v1, const Vector &v2);
    Vector cos(const Vector &v);
    Vector sin(const Vector &v);
    Vector atan2(const Vector &vy, const Vector &vx);
    double atan2(const Vector &v);
    Vector pos_def_mod(const Vector &v, double b);
    Vector wrap(const Vector &v);
    Vector wrap(const Vector &v, double d);
    Vector wrap(const Vector &v, const Vector &bounds);

    /* Vector Operators
     *
     * Implementations for
     * +, +=
     * -, -=
     * *, *=
     * /, /=
     * With both a double and a Vector for each
     *
     * Also
     *  ==
     *  !=
     *  where this takes into account the entire vector
     *
     */
    Vector operator+(const Vector &v, double d);
    Vector operator+(double d, const Vector &v);
    Vector operator+(const Vector &lhs, const Vector &rhs);
    Vector& operator+=(Vector &lhs, double rhs);
    Vector& operator+=(double lhs, Vector &rhs);
    Vector& operator+=(Vector &lhs, const Vector &rhs);

    Vector operator-(const Vector &lhs);
    Vector operator-(const Vector &lhs, double rhs);
    Vector operator-(double lhs, const Vector &rhs);
    Vector operator-(const Vector &lhs, const Vector &rhs);
    Vector& operator-=(Vector &lhs, double rhs);
    Vector& operator-=(Vector &lhs, const Vector &rhs);

    Vector operator*(const Vector &lhs, double rhs);
    Vector operator*(double lhs, const Vector &rhs);
    Vector operator*(const Vector &lhs, const Vector &rhs);
    Vector& operator*=(Vector &lhs, double rhs);
    Vector& operator*=(Vector &lhs, const Vector &rhs);

    Vector operator/(const Vector &lhs, double rhs);
    Vector operator/(double lhs, const Vector &rhs);
    Vector operator/(const Vector &lhs, const Vector &rhs);
    Vector& operator/=(Vector &lhs, double rhs);
    Vector& operator/=(Vector &lhs, const Vector &rhs);

    bool operator==(const Vector &lhs, const Vector &rhs);
    bool operator!=(const Vector &lhs, const Vector &rhs);

    std::ostream& operator<< (std::ostream &os, const Vector& v);
    std::ostream& operator<< (std::ostream &os, Vector& v);
};

#endif //ANALYSIS_VECTOR_H
