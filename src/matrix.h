// matrix.h - Copyright (C) 2004 Reed A. Cartwright (all rights reserved)

#ifndef DAWG_MATRIX_H
#define DAWG_MATRIX_H

class Matrix44;
class Vector4;

int EigenSystem(const Matrix44& rOrig, Vector4& rValues, Matrix44& rVectors);

void MatrixAdd(Matrix44& matA, const Matrix44& matL, const Matrix44& matR);
void MatrixSub(Matrix44& matA, const Matrix44& matL, const Matrix44& matR);
void MatrixScale(Matrix44& matA, const Matrix44& matL, const Matrix44& matR);
void MatrixScale(Matrix44& matA, const Matrix44& matL, const Vector4& matR);
void MatrixScale(Matrix44& matA, const Vector4& matL, const Matrix44& matR);
void MatrixScale(Matrix44& matA, const Matrix44& matL, double scale);
void MatrixMult(Matrix44& matA, const Matrix44& matL, const Matrix44& matR);

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
	Vector4();
	Vector4(const Val *valArray);
	Vector4(const Vector4& vec);
	Vector4(Val a, Val b, Val c, Val d);

	Vector4& operator=(const Vector4& vec);

	static const Vector4 s_Zero;
	static const Vector4 s_One;
	void Zero();
	void One();

	inline Val  Value(Pos c) const { return m_valArray[c]; }
	inline Val& Value(Pos c) { return m_valArray[c]; }
	inline Val  operator()(Pos c) const { return m_valArray[c]; }
	inline Val& operator()(Pos c) { return m_valArray[c]; }
	inline Val  operator[](Pos c) const { return m_valArray[c]; }
	inline Val& operator[](Pos c) { return m_valArray[c]; }
};

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
	Matrix44();
	Matrix44(const Val *valArray);
	Matrix44(const Matrix44& mat);
	Matrix44(const Vector4& mat);
	
	Matrix44& operator=(const Matrix44& mat);

	static const Matrix44 s_Zero;
	static const Matrix44 s_One;
	void Zero();
	void One();

	inline Val  Value(Pos r, Pos c) const { return m_valMatrix[r][c]; }
	inline Val& Value(Pos r, Pos c) { return m_valMatrix[r][c]; }
	inline Val  operator()(Pos r, Pos c) const { return m_valMatrix[r][c]; }
	inline Val& operator()(Pos r, Pos c) { return m_valMatrix[r][c]; }
	inline Vector4& operator[](Pos r) { return *(Vector4*)&m_valMatrix[r]; }
	inline const Vector4&  operator[](Pos r) const { return *(Vector4*)&m_valMatrix[r]; }

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

	Matrix44& Transpose();
	Vector4 Diagonal() const
	{	
		return Vector4( m_valMatrix[0][0], m_valMatrix[1][1],
						m_valMatrix[2][2], m_valMatrix[3][3]);
	}
};
#endif
