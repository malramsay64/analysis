//
//  Vector2d.h
//  analysis
//
//  Created by Malcolm Ramsay on 7/12/2014.
//  Copyright (c) 2014 Malcolm Ramsay. All rights reserved.
//



#ifndef MY_VECT
#define MY_VECT

#include "Vector.h"

/* Defines a two dimensional vector which can be easily used for vector arithmetic
 * and linear algebra. This abstracts the vector operations present in the analysis
 * code.
 *
 * */

namespace LAlgebra {

    class Vector2d {
    public:
        double x;
        double y;

        Vector2d();

        Vector2d(double, double);

        Vector2d(double *v);

        Vector2d(const Vector2d &);

        void normalise();

        double length() const;

        void orthogonalise();

        double angle() const;

        Vector2d operator+=(double);

        Vector2d operator+=(const Vector2d &);

        Vector2d operator+=(int);

        Vector2d operator-() const;

        Vector2d operator/=(double i);

        Vector2d operator*=(double i);

        Vector2d operator/=(const Vector2d &i);

        Vector2d operator*=(const Vector2d &i);

        Vector2d operator+(double i) const;

        Vector2d operator+(const Vector2d &v) const;

        Vector2d operator-(const Vector2d &v) const;
        Vector2d operator-=(const Vector2d &v);

        Vector2d operator-(double i) const;
        Vector2d operator-=(double i);

        Vector2d operator+(int i) const;

        Vector2d operator/(double i) const;

        Vector2d operator*(double i) const;

        Vector2d operator*(const Vector2d &v) const;

        Vector2d operator/(const Vector2d &v) const;

        bool operator==(const Vector2d &v) const;

        bool operator!=(const Vector2d &v) const;

        friend std::ostream &operator<<(std::ostream &os, const Vector2d &v) {
            os << v.x << " " << v.y;
            return os;
        }

    };

    Vector2d operator*(double i, const Vector2d &v);

    double dist(const Vector2d &, const Vector2d &);

    Vector2d direction(const Vector2d &v1, const Vector2d &v2);

    double dot_product(const Vector2d &v1, const Vector2d &v2);

    Vector2d cos(const Vector2d &v);

    Vector2d sin(const Vector2d &v);

    Vector2d tan(const Vector2d &v);

    Vector2d atan2(const Vector2d &vx, const Vector2d &vy);

    Vector2d pos_def_mod(const Vector2d &v, double b);

    double atan2(const Vector2d &v);


};

#endif

