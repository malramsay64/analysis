//
//  Vector3D.cpp
//  analysis
//
//  Created by Malcolm Ramsay on 7/12/2014.
//  Copyright (c) 2014 Malcolm Ramsay. All rights reserved.
//

#include "Vector2d.h"

using namespace std;

Vector2d::Vector2d(){
    x = 0;
    y = 0;
}

Vector2d::Vector2d(double new_x, double new_y){
    x = new_x;
    y = new_y;
}

Vector2d::Vector2d(double *v){
    x = v[XX];
    y = v[YY];
}

Vector2d::Vector2d(const Vector2d &v) {
    x = v.x;
    y = v.y;
}


/* Resizes the vector such that is has length 1 */
void Vector2d::normalise(){
    double l = length();
    x /= l;
    y /= l;
}

/* Reorients the vector to it's right handed normal */
void Vector2d::orthogonalise(){
    double temp = x;
    x = -y;
    y = temp;
}

/* Gives the angle of the vector on the range (-PI,PI] */
double Vector2d::angle() const{
    return atan2(*this);
}

double Vector2d::length() const{
    return sqrt(dot_product(*this, *this));
}

double dot_product(const Vector2d &v1, const Vector2d &v2) {
    return v1.x*v2.x + v1.y*v2.y;
}


/* Finds the vector from the point defined by v1 to that defined by v2
 * This assumes the vectors have been scaled such that they describe points
 * in a box from 0,0 to 2PI,2PI with periodic boundary conditions.
 * The vector given is the shortest distance between the two points,
 * i.e. vector.length() <= box diagonal = sqrt(2*(2PI*2PI))
 */
Vector2d direction(const Vector2d &v1, const Vector2d &v2) {
    Vector2d v = (v2 - v1);
    return atan2(sin(v), cos(v));
}

/* Finds the length of the vector found in the direction function
 */
double dist(const Vector2d &v1, const Vector2d &v2) {
    return direction(v1,v2).length();
}


Vector2d Vector2d::operator+= (double i) {
    x += i;
    y += i;
    return *this;
}

Vector2d Vector2d::operator+= (int i){
    x += i;
    y += i;
    return *this;
}

Vector2d Vector2d::operator+ (const Vector2d &v) const{
    return Vector2d{x + v.x, y + v.y};
}

Vector2d Vector2d::operator+= (const Vector2d &v){
    x += v.x;
    y += v.y;
    return *this;
}

Vector2d Vector2d::operator- () const{
    return Vector2d{-x, -y};
}

Vector2d Vector2d::operator- (const Vector2d &v) const{
    return Vector2d{x - v.x, y - v.y};
}

Vector2d Vector2d::operator- (double i) const{
    return Vector2d{x - i, y - i};
}

Vector2d Vector2d::operator/= (double i){
    x /= i;
    y /= i;
    return *this;
}
Vector2d Vector2d::operator*= (double i){
    x *= i;
    y *= i;
    return *this;
}

Vector2d Vector2d::operator/= (const Vector2d &i){
    x /= i.x;
    y /= i.y;
    return *this;
}

Vector2d Vector2d::operator*= (const Vector2d &i){
    x *= i.x;
    y *= i.y;
    return *this;
}

Vector2d Vector2d::operator* (double i) const{
    return Vector2d{x * i, y * i};
}

Vector2d Vector2d::operator* (const Vector2d &v) const{
    return Vector2d{x * v.x, y * v.y};
}

Vector2d Vector2d::operator/ (const Vector2d &v) const{
    return Vector2d{x / v.x, y / v.y};
}

Vector2d Vector2d::operator+ (double i) const{
    return Vector2d{x + i, y + i};
}

Vector2d Vector2d::operator+ (int i) const{
    return Vector2d{x + i, y + i};
}


Vector2d Vector2d::operator/ (double i) const{
    return Vector2d{x / i, y / i};
}

bool Vector2d::operator== (const Vector2d &v) const{
    return x==v.x && y==v.y;
}
bool Vector2d::operator!= (const Vector2d &v) const{
    return x!=v.x || y!=v.y;
}

Vector2d operator* (double i, const Vector2d &v){
    return Vector2d{v.x * i, v.y * i};
}

Vector2d cos(const Vector2d &v){
    return Vector2d{cos(v.x), cos(v.y)};
}
Vector2d sin(const Vector2d &v){
    return Vector2d{sin(v.x), sin(v.y)};
}

Vector2d atan2(const Vector2d &vx, const Vector2d &vy){
    return Vector2d{atan2(vx.x, vy.x), atan2(vx.y, vy.y)};
}

Vector2d pos_def_mod(const Vector2d &v, double b){
    return Vector2d{my_mod(v.x, b), my_mod(v.y, b)};
}

double atan2(const Vector2d &v){
    return atan2(v.y,v.x);
}

double map_angle(double x){
    return atan2(sin(x), cos(x))+PI;
}



