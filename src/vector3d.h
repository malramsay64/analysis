//
// Created by malcolm on 20/01/16.
//

#ifndef ANALYSIS_VECTOR3D_H
#define ANALYSIS_VECTOR3D_H

#include "Vector.h"

namespace LAlgebra {

    class Vector3d {
        /* TODO
         * create 3d vector class with the same functionality of the
         * 2d vector class. I will probably want to make the 2D vector
         * a subclass of the 3d vector, or both the subclass of a
         * virtual vector class.
         *
         * I will also want to put these classes in their own namespace
         * to prevent collisions with the std::vector class. I could
         * come up with a different name but that is hard!
         */
    public:
        double r[3];

        Vector3d();

        Vector3d(double, double, double);

        Vector3d(double *v);

        Vector3d(const Vector3d &);

        void normalise();

        double length() const;

        void orthogonalise();

        double angle() const;

        Vector3d operator+=(double);

        Vector3d operator+=(const Vector3d &);

        Vector3d operator+=(int);

        Vector3d operator-() const;

        Vector3d operator/=(double i);

        Vector3d operator*=(double i);

        Vector3d operator/=(const Vector3d &i);

        Vector3d operator*=(const Vector3d &i);

        Vector3d operator+(double i) const;

        Vector3d operator+(const Vector3d &v) const;

        Vector3d operator-(const Vector3d &v) const;

        Vector3d operator-(double i) const;

        Vector3d operator+(int i) const;

        Vector3d operator/(double i) const;

        Vector3d operator*(double i) const;

        Vector3d operator*(const Vector3d &v) const;

        Vector3d operator/(const Vector3d &v) const;

        bool operator==(const Vector3d &v) const;

        bool operator!=(const Vector3d &v) const;

        friend std::ostream &operator<<(std::ostream &os, const Vector3d &v) {
            os << v.r[0] << " " << v.r[1] << " " << v.r[2];
            return os;
        }

    };

    Vector3d operator*(double i, const Vector3d &v);

    double dist(const Vector3d &, const Vector3d &);

    Vector3d direction(const Vector3d &v1, const Vector3d &v2);

    double dot_product(const Vector3d &v1, const Vector3d &v2);

    Vector3d cos(const Vector3d &v);

    Vector3d sin(const Vector3d &v);

    Vector3d atan2(const Vector3d &vx, const Vector3d &vy);

    Vector3d pos_def_mod(const Vector3d &v, double b);


};

#endif //ANALYSIS_VECTOR3D_H
