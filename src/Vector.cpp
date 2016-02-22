//
// Created by malcolm on 22/02/16.
//

#include "Vector.h"

namespace LAlgebra {

    double dist(const Vector &from, const Vector &to) {
        return (to - from).length();
    }

    Vector direction(const Vector &v1, const Vector &v2) {
        return v2 - v1;
    }

    double dot_product(const Vector &v1, const Vector &v2) {
        return (v1*v2).sum();
    }

    Vector cos(const Vector &v) {
        return Vector(cos(v.r));
    }

    Vector sin(const Vector &v) {
        Vector vo{};
        vo.r = sin(v.r);
        return vo;
    }

    Vector atan2(const Vector &vy, const Vector &vx) {
        return Vector(atan2(vy.r, vx.r));
    }

    double atan2(const Vector &v) {
        return std::atan2(v[1], v[0]);
    }

    Vector pos_def_mod(const Vector &v, double b) {
        Vector vo{v.r};
        for (auto &a: vo.r) {
            vo.r = my_mod(a, b);
        }
        return vo;
    }

    double Vector::length() const {
        auto sq = r*r;
        double s = sq.sum();
        double sr = sqrt(s);
        return sr;
        return sqrt((r*r).sum());
    }

    void Vector::normalise() {
        double l = length();
        if (l > 0) r /= length();
    }

    size_t Vector::size() const {
        return r.size();
    }

    Vector Vector::orthogonal() const {
        Vector vo{};
        if (size() > 1) {
            vo[0] = -r[1];
            vo[1] = r[0];
        }
        else {
            vo[0] = 0;
        }
        return vo;
    }


    double Vector::sum() const {
        return r.sum();
    }

    /* Returns a vector which represents the current vector
     * in spherical coordinates, reference
     * https://en.wikipedia.org/wiki/N-sphere#Spherical_coordinates
     */
    Vector Vector::angular() const {
        Vector vo{};
        for (auto i = 0; i < vo.r.size(); i++) {
            if (i == 0) {
                vo[i] = length();
            }
            else {
                double sum_sq{0};
                for (int j = size() - 1; j >= i; --j) {
                    sum_sq += r[j] * r[j];
                }
                vo[i] = r[i] / sqrt(sum_sq);
            }
        }
        return vo;
    }

    Vector wrap(const Vector &v) {
        return Vector(v.r.apply([](double val){return map_angle(val);}));
    };

    Vector wrap(const Vector &v, double d) {
        return Vector((v*2.*PI/d).r.apply([](double val){return map_angle(val);})*(d/2.*PI));
    };

    Vector wrap(const Vector &v, const Vector &bounds) {
        return Vector((v*2.*PI/bounds).r.apply([](double val){return map_angle(val);})*(bounds.r/2.*PI));
    };

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
    Vector operator+(const Vector &v, double d) {
        return Vector(v.r + d);
    }

    Vector operator+(double d, const Vector &v) {
        return Vector(d + v.r);
    }

    Vector operator+(const Vector &lhs, const Vector &rhs) {
        return Vector(lhs.r + rhs.r);
    }

    Vector &operator+=(Vector &lhs, double rhs) {
        lhs.r -= rhs;
        return lhs;
    }

    Vector &operator+=(double lhs, Vector &rhs) {
        rhs.r += lhs;
        return rhs;
    }

    Vector &operator+=(Vector &lhs, const Vector &rhs) {
        lhs.r += rhs.r;
        return lhs;
    }

    Vector operator-(const Vector &lhs) {
        return Vector(-lhs.r);
    }

    Vector operator-(const Vector &lhs, double rhs) {
        return Vector(lhs.r - rhs);
    }

    Vector operator-(double lhs, const Vector &rhs) {
        return Vector(lhs - rhs.r);
    }

    LAlgebra::Vector operator-(const LAlgebra::Vector &lhs, const LAlgebra::Vector &rhs) {
        return LAlgebra::Vector(lhs.r - rhs.r);
    }

    Vector &operator-=(Vector &lhs, double rhs) {
        lhs.r -= rhs;
        return lhs;
    }

    Vector &operator-=(Vector &lhs, const Vector &rhs) {
        lhs.r -= rhs.r;
        return lhs;
    }

    Vector operator*(const Vector &lhs, double rhs) {
        return Vector(lhs.r * rhs);
    }

    Vector operator*(double lhs, const Vector &rhs) {
        return Vector(lhs * rhs.r);
    }

    Vector operator*(const Vector &lhs, const Vector &rhs) {
        return Vector(lhs.r * rhs.r);
    }

    Vector &operator*=(Vector &lhs, double rhs) {
        lhs.r *= rhs;
        return lhs;
    }

    Vector &operator*=(Vector &lhs, const Vector &rhs) {
        lhs.r *= rhs.r;
        return lhs;
    }

    Vector operator/(const Vector &lhs, double rhs) {
        return Vector(lhs.r / rhs);
    }

    Vector operator/(double lhs, const Vector &rhs) {
        return Vector(lhs / rhs.r);
    }

    Vector operator/(const Vector &lhs, const Vector &rhs) {
        return Vector(lhs.r / rhs.r);
    }

    Vector &operator/=(Vector &lhs, double rhs) {
        lhs.r /= rhs;
        return lhs;
    }

    Vector &operator/=(Vector &lhs, const Vector &rhs) {
        lhs.r /= rhs.r;
        return lhs;
    }

    bool operator==(const Vector &lhs, const Vector &rhs) {
        return (lhs.r == rhs.r).min();
    };

    bool operator!=(const Vector &lhs, const Vector &rhs) {
        return (lhs.r != rhs.r).min();
    };

    std::ostream &operator<<(std::ostream &os, Vector &v) {
        for (size_t j = 0; j < v.size(); j++) {
            os << v.r[j];
            j < v.size() - 1 ? os << " " : os << "";
        }
        return os;
    }

    std::ostream &operator<<(std::ostream &os, const Vector &v) {
        for (size_t j = 0; j < v.size(); j++) {
            os << v.r[j];
            j < v.size() - 1 ? os << " " : os << "";
        }
        return os;
    }
}