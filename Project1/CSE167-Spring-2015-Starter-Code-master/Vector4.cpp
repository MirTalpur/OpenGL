/***************************************************
 * Mir Ali Talpur
 * Vector4.cpp
 * 4/1/2015
***************************************************/
#include "Vector4.h"
#include "Vector3.h"
#include <math.h>
#include <iostream>
#include <cstring>

/***************************************************
 * Vector4()
 * Constructor for the Vector4 class
***************************************************/
Vector4::Vector4()
{
    std::memset(m, 0, sizeof(m));
}
/***************************************************
 * Vector4(float, float, float)
 * Constructor for the Vector4 class with 3 values
***************************************************/
Vector4::Vector4(float m0, float m1, float m2)
{
    m[0] = m0;
    m[1] = m1;
    m[2] = m2;
    m[3] = 1;
}
/***************************************************
 * Vector4(float, float, float,float)
 * Constructor for the Vector4 class with 4 values
***************************************************/
Vector4::Vector4(float m0, float m1, float m2, float m3)
{
    m[0] = m0;
    m[1] = m1;
    m[2] = m2;
    m[3] = m3;
}
/***************************************************
 * ptr()
 * returns ptr of the vector4 values
***************************************************/
float* Vector4::ptr()
{
    return &m[0];
}
/***************************************************
 * set(float,float,float,float)
 * sets the vector4 array to x y and z w value 
 * passed in
***************************************************/
void Vector4::set(float x, float y, float z, float w)
{
    m[0] = x;
    m[1] = y;
    m[2] = z;
    m[3] = w;
}
/***************************************************
 * overload [] operator 
 * returns the location of the particular 
 * Vector4 component
***************************************************/
float& Vector4::operator [] (int loc)
{
    return m[loc];
}
/***************************************************
 * add(Vector4)
 * add all the vector4 components with the given
 * vector4
 * Basically add two Vector4's together 
***************************************************/
Vector4 Vector4::add(Vector4& a)
{
    Vector4 b;
    b.m[0] = a.m[0] + m[0];
    b.m[1] = a.m[1] + m[1];
    b.m[2] = a.m[2] + m[2];
    b.m[3] = a.m[3] + m[3];
    return b;
}
/***************************************************
 * Operator +
 * Adds two Vector4s together using overloading
 * function
***************************************************/
Vector4 Vector4::operator + (Vector4 a)
{
    return add(a);
}
/***************************************************
 * Operator -
 * Subtract two Vector4s together using overloading
 * function
***************************************************/
Vector4 Vector4::subtract(Vector4& a)
{
    Vector4 b;
    b.m[0] =  m[0] - a.m[0];
    b.m[1] =  m[1] - a.m[1];
    b.m[2] =  m[2] - a.m[2];  
    b.m[3] =  m[3] - a.m[3]; 
    return b;
}
/***************************************************
 * Operator -
 * Subtract two Vector4s together using overloading
 * function
***************************************************/
Vector4 Vector4::operator - (Vector4 a)
{
    return subtract(a);
}
/***************************************************
 * Dehomogenize()
 * Scale the vector so that the fourth component is
 * zero
***************************************************/
Vector4 Vector4::dehomogenize()
{
    Vector4 b;
    if(m[3] == 0){
        m[3] = 1;
    }
    b.m[0] =  m[0] / m[3];
    b.m[1] =  m[1] / m[3];
    b.m[2] =  m[2] / m[3];
    b.m[3] =  m[3] / m[3];
    return b;
}
/***************************************************
 * toVector3(float)
 * Convert the vector4 to Vector3
***************************************************/
Vector3 Vector4::toVector3()
{
    Vector3 b(m[0], m[1], m[2]);
    return b;
}
/***************************************************
 * dot(Vector4)
 * dot product of two vector4s
***************************************************/
float Vector4::dot(Vector4 a)
{
    return (m[0] * a.m[0]) + (m[1] * a.m[1]) + (m[2] * a.m[2]) + (m[3] * a.m[3]);
}

void Vector4::print(std::string comment)
{
    std::cout << comment << std::endl;
    std::cout << "<x:" << m[0] <<  ", y:" << m[1] << ", z:" << m[2] << ", w:" << m[3] << ">" << std::endl;
}