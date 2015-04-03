/***************************************************
 * Mir Ali Talpur
 * Vector3.cpp
 * 4/1/2015
***************************************************/
#include "Vector3.h"
#include <math.h>
#include <iostream>
#include <cstring>
#include "Vector4.h"

/***************************************************
 * Vector3()
 * Constructor for the Vector3 class
***************************************************/
Vector3::Vector3()
{
    std::memset(m, 0, sizeof(m));
}
/***************************************************
 * Vector3(float, float, float)
 * Constructor for the Vector class with 3 values
***************************************************/
Vector3::Vector3(float m0, float m1, float m2)
{
    m[0] = m0;
    m[1] = m1;
    m[2] = m2;
}
/***************************************************
 * Virtual ~Vector3()
 * Negation for Vector3 with no values
***************************************************/
Vector3::~Vector3()
{
	m[0] *= -1;
	m[1] *= -1;
	m[2] *= -1;
}
/***************************************************
 * ptr()
 * returns ptr of the vector values
***************************************************/
float* Vector3::ptr()
{
    return &m[0];
}
/***************************************************
 * set(float,float,float)
 * sets the vector3 array to x y and z value passed
 * in 
***************************************************/
void Vector3::set(float x, float y, float z)
{
    m[0] = x;
    m[1] = y;
    m[2] = z;
}
/***************************************************
 * set(int,float)
 * sets the value of the particular vector3 
 * component with the provided index and value
***************************************************/
void Vector3::set(int index, float value)
{
    m[index] = value;
}
/***************************************************
 * overload [] operator 
 * returns the location of the particular 
 * Vector3 component
***************************************************/
float& Vector3::operator [] (int loc)
{
    return m[loc];
}
/***************************************************
 * add(Vector3)
 * add all the vector3 components with the given
 * vector3 
 * Basically add two Vector3's together 
***************************************************/
Vector3 Vector3::add(Vector3 a)
{
    Vector3 b;
    b.m[0] = a.m[0] + m[0];
    b.m[1] = a.m[1] + m[1];
    b.m[2] = a.m[2] + m[2];
    return b;
}
/***************************************************
 * Operator +
 * Adds two Vectors together using overloading
 * function
***************************************************/
Vector3 Vector3::operator + (Vector3 a)
{
    return add(a);
}
/***************************************************
 * subtract(Vector3)
 * subtract all the vector3 components with the given
 * vector3 
 * Basically subtract two Vector3's together 
***************************************************/
Vector3 Vector3::subtract(Vector3 a)
{
    Vector3 b;
    b.m[0] =  m[0] - a.m[0];
    b.m[1] =  m[1] - a.m[1];
    b.m[2] =  m[2] - a.m[2];
    return b;
}
/***************************************************
 * Operator -
 * Subtract two Vectors together using overloading
 * function
***************************************************/
Vector3 Vector3::operator - (Vector3 a)
{
    return subtract(a);
}
/***************************************************
 * Negate(void)
 * Negation for Vector3 with no values
***************************************************/
Vector3 Vector3::negate(void)
{   
    Vector3 b;
    b.m[0] = m[0] * -1;
    b.m[1] = m[1] * -1;
    b.m[2] = m[2] * -1;
	return b;
}
/***************************************************
 * scale(float)
 * Negation for Vector3 with no values
***************************************************/
Vector3 Vector3::scale(float s)
{
    Vector3 b;
    b.m[0] = m[0] * s;
    b.m[1] = m[1] * s;
    b.m[2] = m[2] * s;
	return b;
}
/***************************************************
 * multiply(float)
 * Multipies vector with float
***************************************************/
Vector3 Vector3::multiply(float a)
{
    Vector3 b;
    b.m[0] =  m[0] * a;
    b.m[1] =  m[1] * a;
    b.m[2] =  m[2] * a;
    return b;
}
/***************************************************
 * operator * (float a)
 * Multilples vector3 with a float value
***************************************************/
Vector3 Vector3::operator * (float a)
{
    return multiply(a);
}
/***************************************************
 * multiply(Vector3)
 * Multilples vector3 with another Vector3
***************************************************/
Vector3 Vector3::multiply(Vector3 a)
{
    Vector3 b;
    b.m[0] =  m[0] * a.m[0];
    b.m[1] =  m[1] * a.m[1];
    b.m[2] =  m[2] * a.m[2];
    return b;
}
/***************************************************
 * operator * overload
 * Multilples vector3 with another Vector3
***************************************************/
Vector3 Vector3::operator * (Vector3 a)
{
    return multiply(a);
}

/***************************************************
 * dot(Vector3)
 * dot product of two vectors
***************************************************/
float Vector3::dot(Vector3 a)
{   
    return ((a.m[0] * m[0]) + (a.m[1] * m[1]) +
            (a.m[2] * m[2]));
}
/***************************************************
 * cross(Vector3)
 * cross product of two vectors3
***************************************************/
Vector3 Vector3::cross(Vector3 a)
{
    Vector3 b;
    b.m[0] = ((m[1] * a.m[2]) - (m[2] * a.m[1]));
    b.m[1] = ((m[2] * a.m[0]) - (m[0] * a.m[2])); 
    b.m[2] = ((m[0] * a.m[1]) - (m[1] * a.m[0]));
    return b;
}
/***************************************************
 * angle(Vector3)
 * angle between two vectors3
***************************************************/
float Vector3::angle(Vector3 a)
{
    return (acos((this->dot(a)) / (this->magnitude() * 
                  a.magnitude()))); 
}
/***************************************************
 * magnitude(void)
 * Magnitude of the Vector3 being used
***************************************************/
float Vector3::magnitude(void)
{
    return sqrt(m[0] * m[0] + m[1] * m[1] + 
                m[2] * m[2]);
}
/***************************************************
 * normailze(void)
 * Normalize the Vector3 being used
***************************************************/
Vector3 Vector3::normalize(void)
{   
    Vector3 b;
    if(this->magnitude() == 0){
        return b;
    }else
        b.m[0] = m[0] / this -> magnitude();
        b.m[1] = m[1] / this -> magnitude();
        b.m[2] = m[2] / this -> magnitude();
    return b;
}
/***************************************************
 * toVector4(float)
 * Convert the vector3 to Vector4
***************************************************/
Vector4 Vector3::toVector4(float w)
{
    Vector4 b(m[0], m[1], m[2], w);
    return b;
}
/***************************************************
 * print(string)
 * print everything and the comment
***************************************************/
void Vector3::print(std::string comment)
{
    std::cout << comment << std::endl;
    std::cout << "<x:" << m[0] <<  ", y:" << m[1] << ", z:" << m[2] << ">" << std::endl;
}
