#ifndef VEC_H
#define VEC_H

#include "globParams.h"

class Vec{
public:
    Vec(){x = 0.0;y = 0.0;};
    Vec(float x0, float y0){x = x0; y = y0;};

    float getX(){return x;};
    float getY(){return y;};
    void setX(float x0){x = x0;};
    void setY(float y0){y = y0;};

    Vec& operator+=(Vec v);
    Vec& operator-=(Vec v);
    Vec& operator*=(float s);
    Vec& operator/=(float s);


private:
    float x;
    float y;

};

inline Vec operator+(Vec a, Vec b){ return a+= b;};
inline Vec operator-(Vec a, Vec b){ return a -= b;};
inline Vec operator*(Vec a, double s){ return a*=s;};
inline Vec operator*(double s, Vec a){ return a*=s;};
inline Vec operator/(Vec a, double s){ return a/=s;};
//inline void print(Vec a){std::cout << a.getX() << "\t" << a.getY() << "\t" << a.getZ() << std::endl;};

#endif //VEC_H