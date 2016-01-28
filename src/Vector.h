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
    public:
        std::valarray<double> r;

        Vector();

        double length() const;

        void normalise();

        void orthogonalise();

        size_t size() const;

        double sum() const;

        Vector<i> apply (double func(double)) const;

        Vector<i> operator+(double) const;

        Vector<i> operator+(const Vector<i> &) const;

        Vector<i>& operator+=(double);

        Vector<i>& operator+=(const Vector<i> &);

        Vector<i> operator-() const;

        Vector<i> operator-(double) const;

        Vector<i> operator-(const Vector<i> &) const;

        Vector<i>& operator-=(double);

        Vector<i>& operator-=(const Vector &v);

        Vector<i> operator*(double) const;

        Vector<i> operator*(const Vector<i> &) const;

        Vector<i>& operator*=(double);

        Vector<i>& operator*=(const Vector<i> &);

        Vector<i> operator/(double) const ;

        Vector<i> operator/(const Vector<i> &) const ;

        Vector<i>& operator/=(double);

        Vector<i> operator/=(const Vector<i> &);

        bool operator==(const Vector &) const ;

        bool operator!=(const Vector &) const;

        Vector<i> cos() const;

        Vector<i> sin() const;

    };


    // Vector operator*(double i, const Vector &v);

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

    template <size_t i>
    Vector<i> Vector<i>::cos()const {
        Vector<i> vo{};
        vo.r = cos(r);
        return vo;
    }

    template <size_t i>
    Vector<i> atan2(const Vector<i> &vx, const Vector<i> &vy);

    template <size_t i>
    Vector<i> operator%(const Vector<i> &v, double b){
        return v.apply([b](double x){return my_mod(x, b);});
    }

    template <size_t i>
    Vector<i> pos_def_mod(const Vector<i> &v, double b){
        return v.apply([b](double x){return my_mod(x, b);});
    }

    template <size_t i>
    Vector<i>::Vector(){
        r = std::valarray<double>(0.,i);
    }

    template <size_t i>
    double Vector<i>::length() const {
        return sqrt(r.apply([](double val) {return val*val;}).sum());
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
    void Vector<i>::orthogonalise() {
        if(size() > 1 ) {
            double tmp = r[0];
            r[0] = -r[1];
            r[1] = tmp;
        }
    }

    template <size_t i>
    Vector<i> Vector<i>::apply(double func(double)) const {
        Vector<i> v{*this};
        v.r.apply(func);
        return v;
    }

    template <size_t i>
    double Vector<i>::sum() const{
        return r.sum();
    }

    template <size_t i>
    Vector<i> Vector<i>::operator+(double d) const {
        return apply([d](double val) {return val+d;});
    }

    template <size_t i>
    Vector<i> Vector<i>::operator+(const Vector<i> &v) const {
        Vector<i> vo{*this};
        vo.r += v.r;
        return vo;
    }

    template <size_t i>
    Vector<i>& Vector<i>::operator+=(double d) {
        r.apply([d](double val) {return val-d;});
        return *this;
    }

    template <size_t i>
    Vector<i>& Vector<i>::operator+=(const Vector<i> &v) {
        r += v.r;
        return *this;
    }

    template <size_t i>
    Vector<i> Vector<i>::operator-() const {
        return apply([](double val) {return -val;});
    }

    template <size_t i>
    Vector<i> Vector<i>::operator-(double d) const {
        return apply([d](double val) {return val-d;});
    }

    template <size_t i>
    Vector<i> Vector<i>::operator-(const Vector<i> &v) const {
        Vector<i> vo{*this};
        vo.r -= v.r;
        return vo;
    }

    template <size_t i>
    Vector<i>& Vector<i>::operator-=(double d) {
        r.apply([d](double val) {return val-d;});
        return *this;
    }

    template <size_t i>
    Vector<i>& Vector<i>::operator-=(const Vector<i> &v) {
        r -= v.r;
        return *this;
    }

    template <size_t i>
    Vector<i> Vector<i>::operator*(double d) const {
        return apply([d](double val) {return val*d;});
    }

    template <size_t i>
    Vector<i> Vector<i>::operator*(const Vector<i> &v) const {
        Vector<i> vo{*this};
        vo.r *= v.r;
        return vo;
    }

    template <size_t i>
    Vector<i>& Vector<i>::operator*=(double d) {
        r.apply([d](double val) {return val*d;});
        return *this;
    }

    template <size_t i>
    Vector<i>& Vector<i>::operator*=(const Vector<i> &v) {
        r *= v.r;
        return *this;
    }

    template <size_t i>
    Vector<i> Vector<i>::operator/(double d) const {
        return apply([d](double val) {return val/d;});
    }

    template <size_t i>
    Vector<i> Vector<i>::operator/(const Vector<i> &v) const {
        Vector<i> vo{*this};
        vo.r /= v.r;
        return vo;
    }

    template <size_t i>
    Vector<i>& Vector<i>::operator/=(double d) {
        r /= d;
        return *this;
    }

    template <size_t i>
    Vector<i> Vector<i>::operator/=(const Vector<i> &v) {
        r /= v.r;
        return *this;
    }

    template <size_t i>
    bool Vector<i>::operator==(const Vector<i> &v) const{
        return (r == v.r).sum() == v.size();
    }

    template <size_t i>
    bool Vector<i>::operator!=(const Vector<i> &v) const{
        return (r == v.r).sum() != v.size();
    }

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