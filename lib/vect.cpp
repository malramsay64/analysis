//
//  vect.cpp
//  analysis
//
//  Created by Malcolm Ramsay on 7/12/2014.
//  Copyright (c) 2014 Malcolm Ramsay. All rights reserved.
//

#include "vect.h"

using namespace std;

vect::vect(){
    x = 0;
    y = 0;
}

vect::vect(double new_x, double new_y){
    x = new_x;
    y = new_y;
}

vect::vect(double *v){
    x = v[XX];
    y = v[YY];
}

void vect::normalise(){
    double l = length();
    x /= l;
    y /= l;
}

void vect::orthogonalise(){
    double temp = x;
    x = -y;
    y = temp;
}

double vect::angle(){
    return atan2(*this) + PI;
}

double vect::length(){
    return sqrt(dot_product(*this, *this));
}

double dot_product(vect v1, vect v2){
    return v1.x*v2.x + v1.y*v2.y;
}

double dist(vect v1, vect v2){
    return direction(v1,v2).length();
}

vect direction(vect v1, vect v2){
    vect v = (v2-v1);
    return atan2(sin(v), cos(v));
}

vect vect::operator+= (double i){
    this->x += i;
    this->y += i;
    return *this;
}

vect vect::operator+= (int i){
    this->x += i;
    this->y += i;
    return *this;
}

vect vect::operator+ ( const vect &v){
    return vect(x+v.x, y+v.y);
}

vect vect::operator+= (const vect v){
    x += v.x;
    y += v.y;
    return *this;
}

vect vect::operator- (){
    return vect(-x,-y);
}

vect vect::operator- (const vect &v){
    return vect(x-v.x, y-v.y);
}

vect vect::operator- (double i){
    return vect(x-i, y-i);
}

vect vect::operator/= (double i){
    return vect(x/i,y/i);
}
vect vect::operator*= (double i){
    return vect(x/i,y/i);
}

vect vect::operator/= (vect &i){
    return vect(x/i.x, y/i.y);
}

vect vect::operator*= (vect &i){
    return vect(x*i.x, y*i.y);
}

vect vect::operator* (double i){
    return vect(x*i, y*i);
}

vect vect::operator* (const vect v){
    return vect(x*v.x, y*v.y);
}

vect vect::operator/ (const vect v){
    return vect(x/v.x, y/v.y);
}

vect vect::operator+ (double i){
    return vect(x+i, y+i);
}

vect vect::operator+ (int i){
    return vect(x+i, y+i);
}


vect vect::operator/ ( double i){
    return vect(x/i, y/i);
}

bool vect::operator== (const vect v){
    return x==v.x && y==v.y;
}
bool vect::operator!= (const vect v){
    return !(x==v.x && y==v.y);
}

vect operator* (double i, vect v){
    return vect(v.x*i, v.y*i);
}

vect cos(vect v){
    return vect(cos(v.x), cos(v.y));
}
vect sin(vect v){
    return vect(sin(v.x), sin(v.y));
}

vect atan2(vect vx, vect vy){
    return vect(atan2(vx.x, vy.x), atan2(vx.y,vy.y));
}

double atan2(vect v){
    return atan2(v.x,v.y);
}

double map_angle(double x){
    return atan2(sin(x), cos(x))+PI;
}




