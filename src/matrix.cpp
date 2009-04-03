// matrix.cc - Copyright (c) 2004-2009 Reed A. Cartwright (all rights reserved)

#include "dawg.h"

#include "matrix.h"

#include <algorithm>

#ifdef HAVE_MEMORY_H
#	include <memory.h>
#endif

using namespace std;

/**********************************************************
	class Vector4
**********************************************************/

const Vector4::Val g_ZeroBaseV[4] = {0.0, 0.0, 0.0, 0.0};
const Vector4::Val g_OneBaseV[4] = {1.0, 1.0, 1.0, 1.0};
const Vector4 Vector4::s_One(&g_OneBaseV[0]);
const Vector4 Vector4::s_Zero(&g_ZeroBaseV[0]);
void  Vector4::Zero() { memcpy(&m_valArray[0], &g_ZeroBaseV[0], 4*sizeof(Val)); }
void  Vector4::One()  { memcpy(&m_valArray[0], &g_OneBaseV[0],  4*sizeof(Val)); }

Vector4::Vector4(Val a, Val b, Val c, Val d)
{
	m_valArray[0] = a;
	m_valArray[1] = b;
	m_valArray[2] = c;
	m_valArray[3] = d;
}

Vector4::Vector4() { One(); }
Vector4::Vector4(const Val *valArray)
{
	memcpy(&m_valArray[0], valArray, 4*sizeof(Val)); 
}
Vector4::Vector4(const Vector4& vec)
{
	memcpy(&m_valArray[0], &vec.m_valArray[0], 4*sizeof(Val));
}

Vector4& Vector4::operator=(const Vector4& vec)
{
	if(&vec != this)
		memcpy(&m_valArray[0], &vec.m_valArray[0], 4*sizeof(Val));
	return *this;
}

/**********************************************************
	class Matrix44
**********************************************************/

const Matrix44::Val g_OneBase[4][4] =  {
	{1.0, 0.0, 0.0, 0.0},
	{0.0, 1.0, 0.0, 0.0},
	{0.0, 0.0, 1.0, 0.0},
	{0.0, 0.0, 0.0, 1.0}};
const Matrix44::Val g_ZeroBase[4][4] =  {
	{0.0, 0.0, 0.0, 0.0},
	{0.0, 0.0, 0.0, 0.0},
	{0.0, 0.0, 0.0, 0.0},
	{0.0, 0.0, 0.0, 0.0}};
const Matrix44 Matrix44::s_One(&g_OneBase[0][0]);
const Matrix44 Matrix44::s_Zero(&g_ZeroBase[0][0]);
void  Matrix44::Zero() { memcpy(&m_valMatrix[0][0], &g_ZeroBase[0][0], 16*sizeof(Val)); }
void  Matrix44::One()  { memcpy(&m_valMatrix[0][0], &g_OneBase[0][0],  16*sizeof(Val)); }

Matrix44::Matrix44(const Val *valArray)
{
	memcpy(&m_valMatrix[0][0], valArray, 16*sizeof(Val));
}

Matrix44::Matrix44()
{
	One();
}

Matrix44::Matrix44(const Vector4 &mat)
{
	Zero();
	m_valMatrix[0][0] = mat[0];
	m_valMatrix[1][1] = mat[1];
	m_valMatrix[2][2] = mat[2];
	m_valMatrix[3][3] = mat[3];
}

Matrix44::Matrix44(const Matrix44& mat)
{
	memcpy(&m_valMatrix[0][0], &mat.m_valMatrix[0][0], 16*sizeof(Val));
}

Matrix44& Matrix44::operator=(const Matrix44& mat)
{
	if(&mat != this)
		memcpy(&m_valMatrix[0][0], &mat.m_valMatrix[0][0], 16*sizeof(Val));
	return *this;
}

Matrix44& Matrix44::Transpose()
{
	for(Pos i=0;i<4;++i)
		for(Pos j=i+1;j<4;++j)
			swap((*this)(i,j), (*this)(j,i));
	return *this;
}

/**********************************************************
	Utility Functions
**********************************************************/
void MatrixAdd(Matrix44& matA, const Matrix44& matL, const Matrix44& matR)
{
	for(Matrix44::Pos i=0;i<4;++i)
		for(Matrix44::Pos j=0;j<4;++j)
			matA(i,j) = matL(i,j)+matR(i,j);
}
void MatrixSub(Matrix44& matA, const Matrix44& matL, const Matrix44& matR)
{
	for(Matrix44::Pos i=0;i<4;++i)
		for(Matrix44::Pos j=0;j<4;++j)
			matA(i,j) = matL(i,j)-matR(i,j);
}

void MatrixScale(Matrix44& matA, const Matrix44& matL, const Matrix44& matR)
{
	for(Matrix44::Pos i=0;i<4;++i)
		for(Matrix44::Pos j=0;j<4;++j)
			matA(i,j) = matL(i,j)*matR(i,j);
}

void MatrixScale(Matrix44& matA, const Matrix44& matL, const Vector4& matR)
{
	for(Matrix44::Pos i=0;i<4;++i)
		for(Matrix44::Pos j=0;j<4;++j)
			matA(i,j) = matL(i,j)*matR(j);
}

void MatrixScale(Matrix44& matA, const Vector4& matL, const Matrix44& matR)
{
	for(Matrix44::Pos i=0;i<4;++i)
		for(Matrix44::Pos j=0;j<4;++j)
			matA(i,j) = matL(i)*matR(i,j);
}

void MatrixScale(Matrix44& matA, const Matrix44& matL, Matrix44::Val scale)
{
	for(Matrix44::Pos i=0;i<4;++i)
		for(Matrix44::Pos j=0;j<4;++j)
			matA(i,j) = matL(i,j)*scale;
}

void MatrixMult(Matrix44& matA, const Matrix44& matL, const Matrix44& matR)
{
	for(Matrix44::Pos i=0;i<4;++i)
		for(Matrix44::Pos j=0;j<4;++j)
			matA(i,j) = matL(i,0)*matR(0,j) + matL(i,1)*matR(1,j)
					  + matL(i,2)*matR(2,j) + matL(i,3)*matR(3,j);
}
