#include "Vec.h"
#include "globParams.h"

Vec& Vec::operator+=(Vec v){
    this->x += v.getX();
    this->y += v.getY();
    return *this;
}

Vec& Vec::operator-=(Vec v){
    this->x -= v.getX();
    this->y -= v.getY();
    return *this;
}


Vec& Vec::operator*=(float s){
    this->x *= s;
    this->y *= s;
    return *this;
}

Vec& Vec::operator/=(float s){
    this->x /= s;
    this->y /= s;
    return *this;
}