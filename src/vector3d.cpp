//
// Created by malcolm on 20/01/16.
//

#include "vector3d.h"

namespace LAlgebra {

    Vector3d::Vector3d() {
        r[0] = r[1] = r[2] = 0;
    }

    Vector3d::Vector3d(double new_x, double new_y, double new_z) {
        r[0] = new_x;
        r[1] = new_y;
        r[2] = new_z;
    }

    Vector3d::Vector3d(double *v) {
        r[0] = v[0];
        r[1] = v[1];
        r[2] = v[2];
    }

    Vector3d::Vector3d(const Vector3d &v) {
        r[0] = v.r[0];
        r[1] = v.r[1];
        r[2] = v.r[2];
    }


/* Resizes the vector such that is has length 1 */
    void Vector3d::normalise() {
        double l = length();
        r[0] /= l;
        r[1] /= l;
        r[2] /= l;
    }

/* Reorients the vector to it's right handed normal */
    void Vector3d::orthogonalise() {
        // TODO Implement this
    }

/* Gives the angle of the vector on the range (-PI,PI] */
    double Vector3d::angle() const {
        // TODO Implement this - will need to change output type
        return 0;

    }

    double Vector3d::length() const {
        return sqrt(dot_product(*this, *this));
    }

    double dot_product(const Vector3d &v1, const Vector3d &v2) {
        return v1.r[0] * v2.r[0] + v1.r[1] * v2.r[1] + v1.r[2] * v2.r[2];
    }


/* Finds the vector from the point defined by v1 to that defined by v2
 * This assumes the vectors have been scaled such that they describe points
 * in a box from 0,0 to 2PI,2PI with periodic boundary conditions.
 * The vector given is the shortest distance between the two points,
 * i.e. vector.length() <= box diagonal = sqrt(2*(2PI*2PI))
 */
    Vector3d direction(const Vector3d &v1, const Vector3d &v2) {
        Vector3d v = (v2 - v1);
        // TODO Box wrapping for 3d vectors
        return v; //atan2(sin(v), cos(v));
    }

/* Finds the length of the vector found in the direction function
 */
    double dist(const Vector3d &v1, const Vector3d &v2) {
        return direction(v1, v2).length();
    }


    Vector3d Vector3d::operator+=(double i) {
        r[0] += i;
        r[1] += i;
        r[2] += i;
        return *this;
    }

    Vector3d Vector3d::operator+=(int i) {
        r[0] += i;
        r[1] += i;
        r[2] += i;
        return *this;
    }

    Vector3d Vector3d::operator+(const Vector3d &v) const {
        return Vector3d{r[0] + v.r[0], r[1] + v.r[1], r[2] + v.r[2]};
    }

    Vector3d Vector3d::operator+=(const Vector3d &v) {
        r[0] += v.r[0];
        r[1] += v.r[1];
        r[2] += v.r[2];
        return *this;
    }

    Vector3d Vector3d::operator-() const {
        return Vector3d{-r[0], -r[1], -r[2]};
    }

    Vector3d Vector3d::operator-(const Vector3d &v) const {
        return Vector3d{r[0] - v.r[0], r[1] - v.r[1], r[2] - v.r[2]};
    }

    Vector3d Vector3d::operator-(double i) const {
        return Vector3d{r[0] - i, r[1] - i, r[2] - i};
    }

    Vector3d Vector3d::operator/=(double i) {
        r[0] /= i;
        r[1] /= i;
        r[2] /= i;
        return *this;
    }

    Vector3d Vector3d::operator*=(double i) {
        r[0] *= i;
        r[1] *= i;
        r[2] *= i;
        return *this;
    }

    Vector3d Vector3d::operator/=(const Vector3d &i) {
        r[0] /= i.r[0];
        r[1] /= i.r[1];
        r[2] /= i.r[2];
        return *this;
    }

    Vector3d Vector3d::operator*=(const Vector3d &i) {
        r[0] *= i.r[0];
        r[1] *= i.r[1];
        r[2] *= i.r[2];
        return *this;
    }

    Vector3d Vector3d::operator*(double i) const {
        return Vector3d{r[0] * i, r[1] * i, r[2] * i};
    }

    Vector3d Vector3d::operator*(const Vector3d &v) const {
        return Vector3d{r[0] * v.r[0], r[1] * v.r[1], r[2] * v.r[2]};
    }

    Vector3d Vector3d::operator/(const Vector3d &v) const {
        return Vector3d{r[0] / v.r[0], r[1] / v.r[1], r[2] / v.r[2]};
    }

    Vector3d Vector3d::operator+(double i) const {
        return Vector3d{r[0] + i, r[1] + i, r[2] + i};
    }

    Vector3d Vector3d::operator+(int i) const {
        return Vector3d{r[0] + i, r[1] + i, r[2] + i};
    }


    Vector3d Vector3d::operator/(double i) const {
        return Vector3d{r[0] / i, r[1] / i, r[2] / i};
    }

    bool Vector3d::operator==(const Vector3d &v) const {
        return r[0] == v.r[0] && r[1] == v.r[1] && r[2] == v.r[2];
    }

    bool Vector3d::operator!=(const Vector3d &v) const {
        return r[0] != v.r[0] || r[1] != v.r[1] || r[2] != v.r[2];
    }

    Vector3d operator*(double i, const Vector3d &v) {
        return Vector3d{v.r[0] * i, v.r[1] * i, v.r[2] * i};
    }

    Vector3d cos(const Vector3d &v) {
        return Vector3d{std::cos(v.r[0]), std::cos(v.r[1]), std::cos(v.r[2])};
    }

    Vector3d sin(const Vector3d &v) {
        return Vector3d{std::sin(v.r[0]), std::sin(v.r[1]), std::sin(v.r[2])};
    }

    Vector3d atan2(const Vector3d &vx, const Vector3d &vy) {
        return Vector3d{std::atan2(vx.r[0], vy.r[0]), std::atan2(vx.r[1], vy.r[1]), std::atan2(vx.r[2], vy.r[2])};
    }

    Vector3d pos_def_mod(const Vector3d &v, double b) {
        return Vector3d{my_mod(v.r[0], b), my_mod(v.r[1], b), my_mod(v.r[2], b)};
    }
};