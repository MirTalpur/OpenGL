/***************************************************
 * Mir Ali Talpur
 * Matrix4.cpp
 * 4/1/2015
***************************************************/
#include <math.h>
#include <iostream>
#include <iomanip>
#include <cstring>
#include "Matrix4.h"
#include "Vector4.h"
#include "Vector3.h"
/***************************************************
 * Matrix4()
 * Constructor for the Matrix4 class
***************************************************/
Matrix4::Matrix4()
{
    std::memset(m, 0, sizeof(m));
}
/***************************************************
 * Matrix4(34 floats)
 * Constructor for the Matrix class with 34 values
***************************************************/
Matrix4::Matrix4(
                 float m00, float m01, float m02, float m03,
                 float m10, float m11, float m12, float m13,
                 float m20, float m21, float m22, float m23,
                 float m30, float m31, float m32, float m33 )
{ 
    this->set(m00, m01, m02, m03,
              m10, m11, m12, m13,
              m20, m21, m22, m23,
              m30, m31, m32, m33);
}
/***************************************************
 * Set(34 floats)
 * Sets all the matrix4 values correspondingly 
***************************************************/
void Matrix4::set(float m00, float m01, float m02, float m03,
                  float m10, float m11, float m12, float m13,
                  float m20, float m21, float m22, float m23,
                  float m30, float m31, float m32, float m33)
{
    m[0][0] = m00;
    m[0][1] = m01;
    m[0][2] = m02;
    m[0][3] = m03;
    m[1][0] = m10;
    m[1][1] = m11;
    m[1][2] = m12;
    m[1][3] = m13;
    m[2][0] = m20;
    m[2][1] = m21;
    m[2][2] = m22;
    m[2][3] = m23;
    m[3][0] = m30;
    m[3][1] = m31;
    m[3][2] = m32;
    m[3][3] = m33;
}
/***************************************************
 * get(int , element)
 * Receive the value asked for in the matrix4
***************************************************/
float Matrix4::get(int vector,int element)
{
    return m[vector][element];
}
/***************************************************
 * Operator = 
 * Copis all the values of MAtrix4 to other_matrix4
***************************************************/
Matrix4& Matrix4::operator=(Matrix4 a)
{
    std::memcpy(m, a.m, sizeof(m));
    return *this;
}
/***************************************************
 * prt()
 * Returns the addres location of the Matrix4
***************************************************/
float* Matrix4::ptr()
{
    return &m[0][0];
}
/***************************************************
 * identity()
 * Sets the Matrix4 as the Identity Matrix
***************************************************/
void Matrix4::identity()
{
    static const float ident[4][4]={{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}};
    std::memcpy(m, ident, sizeof(m));
}
/***************************************************
 * multiply(Matrix4)
 * multiply the matrix4 with another matrix4 
***************************************************/
Matrix4 Matrix4::multiply(Matrix4 a)
{
    Matrix4 b;
    /********************************************************************
     *essentially go thru the matrix size and keep a count of the values
     *being added together. And in the inner most loop take the element
     *in the first row multiple it by the first element in the column.
     *Keep a count of those values and do it with the second element 
     *and add both the counts together and just keep doing that for the
     *size of the matrix.
    ********************************************************************/
    for (int row = 0; row < 4; row++) {
        for (int col = 0; col < 4; col++) {
            double count = 0;
            for (int inner = 0; inner < 4; inner++){ 
                count += m[row][inner] * a.m[inner][col]; 
            }
            b.m[row][col] = count;
        }
    }
    return b;
}
/***************************************************
 * Operator * Overload for matrix multipication
 * multiply two matrices together
***************************************************/
Matrix4 Matrix4::operator * (Matrix4 a)
{
    return multiply(a);
}
/***************************************************
 * multiply(Vectror4)
 * Multiply matrix4 by a Vector4
***************************************************/
Vector4 Matrix4::multiply(Vector4 a)
{
    Vector4 b;
    for(int i = 0; i < 4; i++){
        for(int j = 0; j < 4; j++){
            b.m[i] += m[i][j] * a.m[j];
        }
    }
    return b;
}
/***************************************************
 * Operator * Overload for matrix multipication
 * multiply matrices by a Vector4
***************************************************/
Vector4 Matrix4::operator * (Vector4 a)
{
    return multiply(a);
}
/***************************************************
 * Multiply(Vector3)
 * multiply matrices with a Vectro3 
***************************************************/
Vector3 Matrix4::multiply(Vector3 a)
{
    Vector3 b;
    for(int i = 0; i < 3; i++){
        for(int j = 0; j < 3; j++){
            b.m[i] += m[i][j] * a.m[j];
        }
    }
    return b;
}
/***************************************************
 * Operator Overload *
 * multiply matrices with a Vectro3 
***************************************************/
Vector3 Matrix4::operator * (Vector3 a)
{
    return multiply(a);
}
/***************************************************
 * makeRotateX(float)
 * rotate the cube about the X axis
***************************************************/
Matrix4 Matrix4::makeRotateX(float angle)
{   
    identity();
    m[1][1] = cos(angle);
    m[1][2] = -sin(angle);
    m[2][1] = sin(angle);
    m[2][2] = cos(angle); 
    return *this;
}
/***************************************************
 * makeRotateY(float)
 * rotate the cube about the Y axis
***************************************************/
Matrix4 Matrix4::makeRotateY(float angle)
{
    identity();
    m[0][0] = cos(angle);
    m[0][2] = sin(angle);
    m[2][0] = -sin(angle);
	m[2][2] = cos(angle);
    return *this;
}
/***************************************************
 * makeRotateZ(float)
 * rotate the cube about the Z axis
***************************************************/
Matrix4 Matrix4::makeRotateZ(float angle)
{
    identity();
    m[0][0] = cos(angle);
    m[0][1] = -sin(angle);
    m[1][0] = sin(angle);
	m[1][1] = cos(angle);
    return *this;
}
/***************************************************
 * makeRotateArbitary(float)
 * rotate the cube an arbitrary axis
***************************************************/
Matrix4 Matrix4::makeRotateArbitrary(Vector3 a, float angle)
{
    identity();
    m[0][0] = 1 + (1 - cos(angle)) * (pow(a.m[0], 2) - 1);
    m[0][1] = -a.m[2] * sin(angle) + (1 - cos(angle)) * a.m[0] * a.m[1];
    m[0][2] = a.m[1] * sin(angle)  + (1 - cos(angle)) * a.m[0] * a.m[2];
    m[1][0] = a.m[2] * sin(angle) + (1 - cos(angle)) * a.m[1] * a.m[0];
    m[1][1] = 1 + (1 - cos(angle)) * (pow(a.m[1], 2) - 1);
    m[1][2] = -a.m[0] * sin(angle) + (1 - cos(angle)) * a.m[1] * a.m[2];
    m[2][0] = -a.m[1] * sin(angle) + (1 - cos(angle)) * a.m[2] * a.m[0];
    m[2][1] = a.m[0] * sin(angle) + (1 - cos(angle)) * a.m[2] * a.m[1];
    m[2][2] = 1 + (1 - cos(angle)) * (pow(a.m[2], 2) - 1);
    return *this;
}
/***************************************************
 * makeScale(float)
 * scale the cube with the given values
***************************************************/
Matrix4 Matrix4::makeScale(float s)
{
    return makeScale(s, s, s);
}
/***************************************************
 * makeScale(float,float,float)
 * scale the cube with the given values
***************************************************/
Matrix4 Matrix4::makeScale(float sx, float sy, float sz)
{
    identity();
    m[0][0] = sx;
    m[1][1] = sy;
    m[2][2] = sz;
    return *this;
}
/***************************************************
 * makeTranslate(float,float,float)
 * translate the cube with the given values
***************************************************/
Matrix4 Matrix4::makeTranslate(float x, float y, float z)
{
    identity();
    m[0][3] = x;
    m[1][3] = y;
    m[2][3] = z;
    return *this;
}
/***************************************************
 * makeTranslate(Vector3)
 * translate the cube with the given values
***************************************************/
Matrix4 Matrix4::makeTranslate(Vector3 a)
{
    return makeTranslate(a[0], a[1], a[2]);
}
/***************************************************
 * transpose(void)
 * transpose the matrix basically the row and colmn
 * interchange
***************************************************/
Matrix4 Matrix4::transpose(void)
{
    Matrix4 b;
    for(int x = 0; x < 4; ++x)
    {
        for(int y = 0; y < 4; ++y)
        {
            b.m[y][x] = m[x][y];
        }
    }
    return b;
}

//Hint: Try basing this on code by cool people on the internet
//In this class it is okay to use code from the internet
//So long as you fully understand the code and can clearly explain it if asked to
//http://stackoverflow.com/questions/2624422/efficient-4x4-matrix-inverse-affine-transform
Matrix4 Matrix4::inverse(void)
{
    Matrix4 b;
   //TODO  
    return b;
}

Matrix4 Matrix4::orthoNormalInverse(void)
{
    Matrix4 b;
    
    //Calculate the inverse of this matrix with the assumption that it is ortho-normal
    //This will be useful when implementing cameras!
    
    return b;
}

void Matrix4::print(std::string comment)
{
    //Width constants and variables
    static const int pointWidth = 1;
    static const int precisionWidth = 4;
    int integerWidth = 1;
    
    //Determine the necessary width to the left of the decimal point
    float* elementPtr = (float*)m;
    float maxValue = fabsf(*(elementPtr++));
    while(elementPtr++ < ((float*)m+16)) if(fabsf(*elementPtr) > maxValue) maxValue = fabsf(*elementPtr);
    while(maxValue >= 10.0) { ++integerWidth; maxValue /= 10.0; }
    
    //Sum up the widths to determine the cell width needed
    int cellWidth = integerWidth + pointWidth + precisionWidth;
    
    //Set the stream parameters for fixed number of digits after the decimal point
    //and a set number of precision digits
    std::cout << std::fixed;
    std::cout << std::setprecision(precisionWidth);
    
    //Print the comment
    std::cout << comment << std::endl;
    
    //Loop through the matrix elements, format each, and print them to screen
    float cellValue;
    for(int element = 0; element < 4; element++)
    {
        std::cout << std::setw(1) << (element == 0 ? "[" : " ");
        for(int vector = 0; vector < 4; vector++)
        {
            cellValue =  m[vector][element];
            std::cout << std::setw(cellWidth + (cellValue >= 0.0 ? 1 : 0)) << cellValue;
            std::cout << std::setw(0) << (vector < 3 ? " " : "");
        }
        std::cout << std::setw(1) << (element == 3 ? "]" : " ") << std::endl;
    }
}
