//
// Created by malcolm on 20/01/16.
//

#ifndef ANALYSIS_VECTOR_H
#define ANALYSIS_VECTOR_H

#include <valarray>
#include "functions.h"

namespace LAlgebra {

    template <size_t i>
    class Vector {
    private:

    public:
        std::valarray<double> r;

        Vector();
        Vector(const std::valarray<double> &d) : r(d) {};
        Vector(std::initializer_list<double> d): r(d) {};
        Vector(const Vector<i> &v) : r(v.r) {};

        double length() const;
        void normalise();
        Vector<i> orthogonal() const;
        size_t size() const;
        double sum() const;
        Vector<i> angular() const;
        double& operator[] (size_t s){
            return r[s];
        }
        double operator[] (size_t s) const{
            return r[s];
        }
    };

    template <size_t i>
    double dist(const Vector<i> &from, const Vector<i> &to){
        return (to - from).length();
    }

    template <size_t i>
    Vector<i> direction(const Vector<i> &v1, const Vector<i> &v2){
        return v2-v1;
    }

    template <size_t i>
    double dot_product(const Vector<i> &v1, const Vector<i> &v2){
        return (v1*v2).sum();
    }

    template <size_t N>
    Vector<N> cos(const Vector<N> &v){
        return Vector<N>(cos(v.r));
    }

    template <size_t i>
    Vector<i> sin(const Vector<i> &v){
        Vector<i> vo{};
        vo.r = sin(v.r);
        return vo;
    }

    template <size_t i>
    Vector<i> atan2(const Vector<i> &vy, const Vector<i> &vx){
        return Vector<i>(atan2(vy.r, vx.r));
    }

    template <size_t N>
    double atan2(const Vector<N> &v){
        return std::atan2(v[1], v[0]);
    }

    template <size_t i>
    Vector<i> pos_def_mod(const Vector<i> &v, double b){
        Vector<i> vo{v.r};
        for (auto &a: vo.r){
            vo.r = my_mod(a,b);
        }
        return vo;
    }

    template <size_t i>
    Vector<i>::Vector(){
        r = std::valarray<double>(0.,i);
    }


    template <size_t i>
    double Vector<i>::length() const {
        return sqrt((r*r).sum());
    }


    template  <size_t i>
    void Vector<i>::normalise() {
        if (length() > 0) r /= length();
    }

    template <size_t i>
    size_t Vector<i>::size() const {
        return i;
    }

    template <size_t i>
    Vector<i> Vector<i>::orthogonal() const{
        Vector<i> vo{};
        if(size() > 1 ) {
            vo[0] = -r[1];
            vo[1] = r[0];
        }
        else {
            vo[0] = 0;
        }
        return vo;
    }

    template <size_t i>
    double Vector<i>::sum() const{
        return r.sum();
    }

    /* Returns a vector which represents the current vector
     * in spherical coordinates, reference
     * https://en.wikipedia.org/wiki/N-sphere#Spherical_coordinates
     */
    template <size_t N>
    Vector<N> Vector<N>::angular() const {
        Vector<N> vo{};
        for (auto i=0; i<vo.r.size(); i++){
            if (i == 0){
                vo[i] = length();
            }
            else {
                double sum_sq{0};
                for (int j = N-1; j >= i; --j) {
                    sum_sq += r[j]*r[j];
                }
                vo[i] = r[i]/sqrt(sum_sq);
            }
        }
        return vo;
    }

    template <size_t i>
    Vector<i> wrap(const Vector<i> &v){
        return Vector<i>(v.r.apply([](double val){return map_angle(val);}));
    };

    template <size_t i>
    Vector<i> wrap(const Vector<i> &v, double d){
        return Vector<i>((v*2*PI/d).r.apply([](double val){return map_angle(val);})*(d/2*PI));
    };

    template <size_t i>
    Vector<i> wrap(const Vector<i> &v, const Vector<i> &bounds){
        return Vector<i>((v*2*PI/bounds).r.apply([](double val){return map_angle(val);})*(bounds.r/2*PI));
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
    template <size_t N>
    Vector<N> operator+(const Vector<N> &v, double d){
        return Vector<N>(v.r+d);
    }

    template <size_t N>
    Vector<N> operator+(double d, const Vector<N> &v){
        return Vector<N>(d+v.r);
    }
    template <size_t N>
    Vector<N> operator+(const Vector<N> &lhs, const Vector<N> &rhs) {
        return Vector<N>(lhs.r+rhs.r);
    }
    template <size_t N>
    Vector<N>& operator+=(Vector<N> &lhs, double rhs){
        lhs.r -= rhs;
        return lhs;
    }

    template <size_t N>
    Vector<N>& operator+=(double lhs, Vector<N> &rhs){
        rhs.r += lhs;
        return rhs;
    }
    template <size_t N>
    Vector<N>& operator+=(Vector<N> &lhs, const Vector<N> &rhs){
        lhs.r += rhs.r;
        return lhs;
    }

    template <size_t N>
    Vector<N> operator-(const Vector<N> &lhs){
        return Vector<N>(-lhs.r);
    }

    template <size_t N>
    Vector<N> operator-(const Vector<N> &lhs, double rhs){
        return Vector<N>(lhs.r-rhs);
    }

    template <size_t N>
    Vector<N> operator-(double lhs, const Vector<N> &rhs){
        return Vector<N>(lhs-rhs.r);
    }

    template <size_t N>
    Vector<N> operator-(const Vector<N> &lhs, const Vector<N> &rhs){
        return Vector<N>(lhs.r-rhs.r);
    }

    template <size_t N>
    Vector<N>& operator-=(Vector<N> &lhs, double rhs){
        lhs.r -= rhs;
        return lhs;
    }

    template <size_t N>
    Vector<N>& operator-=(Vector<N> &lhs, const Vector<N> &rhs){
        lhs.r -= rhs.r;
        return lhs;
    }

    template <size_t N>
    Vector<N> operator*(const Vector<N> &lhs, double rhs){
        return Vector<N>(lhs.r*rhs);
    }

    template <size_t N>
    Vector<N> operator*(double lhs, const Vector<N> &rhs){
        return Vector<N>(lhs*rhs.r);
    }

    template <size_t N>
    Vector<N> operator*(const Vector<N> &lhs, const Vector<N> &rhs){
        return Vector<N>(lhs.r*rhs.r);
    }

    template <size_t N>
    Vector<N>& operator*=(Vector<N> &lhs, double rhs){
        lhs.r *= rhs;
        return lhs;
    }

    template <size_t N>
    Vector<N>& operator*=(Vector<N> &lhs, const Vector<N> &rhs){
        lhs.r *= rhs.r;
        return lhs;
    }

    template <size_t N>
    Vector<N> operator/(const Vector<N> &lhs, double rhs){
        return Vector<N>(lhs.r/rhs);
    }

    template <size_t N>
    Vector<N> operator/(double lhs, const Vector<N> &rhs){
        return Vector<N>(lhs/rhs.r);
    }

    template <size_t N>
    Vector<N> operator/(const Vector<N> &lhs, const Vector<N> &rhs){
        return Vector<N>(lhs.r/rhs.r);
    }

    template <size_t N>
    Vector<N>& operator/=(Vector<N> &lhs, double rhs){
        lhs.r /= rhs;
        return lhs;
    }

    template <size_t N>
    Vector<N>& operator/=(Vector<N> &lhs, const Vector<N> &rhs){
        lhs.r /= rhs.r;
        return lhs;
    }

    template <size_t N>
    bool operator==(const Vector<N> &lhs, const Vector<N> &rhs){
        return (lhs.r == rhs.r).min();
    };

    template <size_t N>
    bool operator!=(const Vector<N> &lhs, const Vector<N> &rhs){
        return (lhs.r != rhs.r).min();
    };

    template <size_t i>
    std::ostream& operator<< (std::ostream &os, Vector<i>& v){
        for (size_t j=0; j < i; j++){
            os << v.r[j];
            j < i-1 ? os << " ": os << "";
        }
        return os;
    }
};

#endif //ANALYSIS_VECTOR_H
