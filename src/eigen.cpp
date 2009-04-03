// eigen.cc - Copyright (c) 2004-2009 Reed A. Cartwright (all rights reserved)
#ifdef HAVE_CONFIG_H
#	include "config.h"
#endif

#include "dawg.h"
#include "matrix.h"

#include <iostream>

#ifdef HAVE_FLOAT_H
#	include <float.h>
#endif
#ifdef HAVE_MATH_H
#	include <math.h>
#endif

#define LAMBDA_THREASHOLD	FLT_EPSILON
#define SUM_THREASHOLD		DBL_EPSILON
#define THETA_TH1 1.3407807929942596e+154
#define THETA_TH2 1.4916681462400413e-154

using namespace std;

// Various Tests
inline double nearZero(double d)
{
	return (-LAMBDA_THREASHOLD < d && d < LAMBDA_THREASHOLD) ? 0.0f : d;
}
inline bool epsTest(double d) { return (-DBL_EPSILON <= d && d <= DBL_EPSILON); }
inline bool thetaTest(double d) { return ( d < -THETA_TH1 || THETA_TH1 < d); }


// Jacobi Rotation for a 4x4 matrix with double precision
// Adapted from Numerical Recipies for C
void JacobiRot44(int p, int q, Matrix44& a, Matrix44& v) // q > p
{
	if(a(p,q) == 0.0)
		return;
	double dTan;
	double dTheta = 0.5*(a(q,q)-a(p,p))/a(p,q);

	if(thetaTest(dTheta))
		dTan = 0.5/dTheta;
	else if(epsTest(dTheta*dTheta))
		dTan = 1.0;
	else
		dTan = 1.0/(fabs(dTheta)+sqrt(dTheta*dTheta+1.0));
	dTan = copysign(dTan, dTheta);
	double dCos = 1.0/sqrt(dTan*dTan+1.0);
	double dSin = dCos*dTan;
	double dTau = dSin/(1.0+dCos);
	a(p,p) -= dTan*a(p,q);
	a(q,q) += dTan*a(p,q);
	a(p,q) = 0.0;
	double g,h;
	for(int r=0;r<p;++r)
	{
		g=a(r,p); h=a(r,q);
		a(r,p) = g - dSin*(h+dTau*g);
		a(r,q) = h + dSin*(g-dTau*h);
	}
	for(int r=p+1;r<q;++r)
	{
		g=a(p,r); h=a(r,q);
		a(p,r) = g - dSin*(h+dTau*g);
		a(r,q) = h + dSin*(g-dTau*h);
	}
	for(int r=q+1;r<4;++r)
	{
		g=a(p,r); h=a(q,r);
		a(p,r) = g - dSin*(h+dTau*g);
		a(q,r) = h + dSin*(g-dTau*h);
	}
	for(int r=0;r<4;++r)
	{
		g=v(r,p); h=v(r,q);
		v(r,p) = g - dSin*(h+dTau*g);
		v(r,q) = h + dSin*(g-dTau*h);
	}
}

// Calculate the Eigensystem of a 4x4 matrix
int EigenSystem(const Matrix44& rOrig, Vector4& rValues, Matrix44& rVectors)
{
	rValues.Zero();
	rVectors.One();
	int nRot = 0;

	Matrix44 a(rOrig); // Copy Input Matrix

	double sm, tresh, d;
	for(int g=0;g<50;++g)
	{
		// Test for convergence
		sm = (fabs(a(0,1)) + fabs(a(0,2)) + fabs(a(0,3))
			+ fabs(a(1,2)) + fabs(a(1,3)) + fabs(a(2,3)));
		if(sm == 0.0)
		{
			rValues[0] = nearZero(a(0,0));
			rValues[1] = nearZero(a(1,1));
			rValues[2] = nearZero(a(2,2));
			rValues[3] = nearZero(a(3,3));
			for(Matrix44::Pos i=0;i<3;++i)
				for(Matrix44::Pos j=0;j<3;++j)
					rVectors(i,j) = nearZero(rVectors(i,j));
			return nRot;
		}
		tresh = ((g < 4) ? 0.0125*sm : 0.0);

		// Rotate each off diagonal if greater than threshold
		for(int i=0;i<4;++i)
		{
			for(int j=i+1;j<4;++j)
			{
				d = 100.0*fabs(a(i,j));
				if(g > 4 && fabs(a(i,i))+d == fabs(a(i,i)) && fabs(a(j,j))+d == fabs(a(j,j)))
					a(i,j) = 0.0;
				else if( fabs(a(i,j)) > tresh)
					JacobiRot44(i, j, a, rVectors);
			}
		}
		nRot+=6;
	}
	return -1; // too many iterations
}
