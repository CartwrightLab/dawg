// matrix.h - Copyright (c) 2004-2009 Reed A. Cartwright (all rights reserved)

#ifndef DAWG_MATRIX_H
#define DAWG_MATRIX_H

// a 4x4 matrix
class Matrix44;

// a 4 vector
class Vector4;

// Find eigensystem of a 4x4 matrix
int EigenSystem(const Matrix44& rOrig, Vector4& rValues, Matrix44& rVectors);

// A(x,y) = L(x,y)+R(x,y)
void MatrixAdd(Matrix44& matA, const Matrix44& matL, const Matrix44& matR);
// A(x,y) = L(x,y)-R(x,y)
void MatrixSub(Matrix44& matA, const Matrix44& matL, const Matrix44& matR);
// A(x,y) = L(x,y)*R(x,y)
void MatrixScale(Matrix44& matA, const Matrix44& matL, const Matrix44& matR);
// A(x,y) = L(x,y)*R(x)
void MatrixScale(Matrix44& matA, const Matrix44& matL, const Vector4& matR);
// A(x,y) = L(y)*R(x,y)
void MatrixScale(Matrix44& matA, const Vector4& matL, const Matrix44& matR);
// A(x,y) = L(x,y)*s
void MatrixScale(Matrix44& matA, const Matrix44& matL, double scale);
// A = L*R
void MatrixMult(Matrix44& matA, const Matrix44& matL, const Matrix44& matR);

// a 4 vector
class Vector4
{
public:
	typedef double Val;
	typedef unsigned int Pos;

private:
	union
	{
		Val m_valArray[4];
	};

public:
	// Vector 4 Construtors
	Vector4();
	Vector4(const Val *valArray);
	Vector4(const Vector4& vec);
	Vector4(Val a, Val b, Val c, Val d);

	// Assignment Operator
	Vector4& operator=(const Vector4& vec);

	// Special Vectors
	static const Vector4 s_Zero;
	static const Vector4 s_One;
	
	// set each element of vector to 0 
	void Zero();
	// set each element of vector to 1
	void One();
	
	// Get Value at position 
	inline Val  Value(Pos c) const { return m_valArray[c]; }
	inline Val& Value(Pos c) { return m_valArray[c]; }
	inline Val  operator()(Pos c) const { return m_valArray[c]; }
	inline Val& operator()(Pos c) { return m_valArray[c]; }
	inline Val  operator[](Pos c) const { return m_valArray[c]; }
	inline Val& operator[](Pos c) { return m_valArray[c]; }
};

// a 4x4 matrix
class Matrix44
{
public:
	typedef double Val;
	typedef unsigned int Pos;

private:
	union
	{
		Val m_valMatrix[4][4];
	};

public:
	// construtors
	Matrix44();
	Matrix44(const Val *valArray);
	Matrix44(const Matrix44& mat);
	Matrix44(const Vector4& mat);
	
	// assignment operator
	Matrix44& operator=(const Matrix44& mat);
	
	// special matricies
	static const Matrix44 s_Zero;
	static const Matrix44 s_One;
	// set each element to zero
	void Zero();
	// set each element to 1
	void One();

	// Get element at position (r,c)
	inline Val  Value(Pos r, Pos c) const { return m_valMatrix[r][c]; }
	inline Val& Value(Pos r, Pos c) { return m_valMatrix[r][c]; }
	inline Val  operator()(Pos r, Pos c) const { return m_valMatrix[r][c]; }
	inline Val& operator()(Pos r, Pos c) { return m_valMatrix[r][c]; }
	
	// Get row r
	inline Vector4& operator[](Pos r) { return *(Vector4*)&m_valMatrix[r]; }
	inline const Vector4&  operator[](Pos r) const { return *(Vector4*)&m_valMatrix[r]; }

	// Member Template Functions

	template<class L, class R>
	inline Matrix44& Multiply(L matL, R matR)
	{
		MatrixMult(*this, matL, matR);
		return *this;
	}
	
	template<class L, class R>
	inline Matrix44& Scale(L matL, R matR)
	{
		MatrixScale(*this, matL, matR);
		return *this;
	}

	template<class L, class R>
	inline Matrix44& Add(L matL, R matR)
	{
		MatrixAdd(*this, matL, matR);
		return *this;
	}

	template<class L, class R>
	inline Matrix44& Sub(L matL, R matR)
	{
		MatrixSub(*this, matL, matR);
		return *this;
	}
	
	// Transpose Matrix
	Matrix44& Transpose();

	// Construct a vector from the diagonal of the matrix
	Vector4 Diagonal() const
	{	
		return Vector4( m_valMatrix[0][0], m_valMatrix[1][1],
						m_valMatrix[2][2], m_valMatrix[3][3]);
	}
};
#endif
