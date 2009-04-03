// rand.h - Copyright (c) 2004-2009 Reed A. Cartwright (all rights reserved)

#ifndef DAWG_RAND_H
#define DAWG_RAND_H

#ifdef HAVE_CONFIG_H
#	include "config.h"
#endif

#include <cmath>
#include <cfloat>

#ifdef HAVE_STDINT_H
#	include <stdint.h>
#elif defined(_MSC_VER)
#	include "compat/stdint.h"
#endif

#ifdef HAVE_SYS_TYPES_H
#	include <sys/types.h>
#endif

#define RDBL32_MIN  2.3283064365386963e-010
#define RDBL_MIN	RDBL32_MIN

#ifndef M_E
#	define M_E        2.71828182845904523536
#endif

/************************************************
	Mersenne Twister Generator
************************************************/

// Seed random number generator from an array
void mt_srand(uint32_t uKeys[], size_t uLen);

// Draw integer from [0,0xFFFFFFFF]
uint32_t mt_rand();

/************************************************
	Generic Random number generators
************************************************/

// Draw integer from [0,0xFFFFFFFF]
inline uint32_t rand_uint() { return mt_rand(); }

// Draw integer from [0,uMax]
inline uint32_t rand_uint(uint32_t uMax)
{
		uint32_t uMask = uMax;
		uMask |= uMask >> 1;
		uMask |= uMask >> 2;
		uMask |= uMask >> 4;
		uMask |= uMask >> 8;
		uMask |= uMask >> 16;
		uint32_t u;
		do { u = rand_uint() & uMask; } while ( u > uMax);
		return u;
}

// Draw floating point from [0.0,1.0) with 32-bit precision
inline double rand_real()
{
	uint32_t u = rand_uint();
	return (double) u*(1.0/4294967296.0);
}

// Draw from Bernolli
inline bool rand_bool(double p)
{
	return (rand_real() < p);
}

// Draw Exponential w/ mean 1.0
inline double rand_exp() { return -log(1.0-rand_real()); }

// Draw Exponential w/ mean L
inline double rand_exp(double L) { return L*rand_exp(); }

inline double rand_gamma_small(double a)
{
	// Ahrens & Dieter (1974) Computer methods for sampling from gamma, 
	// beta, Poisson and binomial distributions. Computing 12: 223-246.
    double b = 1.0 + a*0.36788794412;
	double p,g,r;
	do
	{
		p = b*rand_real();
		if(p < 1.0)
		{
			g = exp(log(p)/a);
			r = -g;
		}
		else
		{
			g = -log((b-p)/a);
			r = ((a-1.0)*log(g));
		}
	}	while(log(1.0-rand_real()) > r);

	return g;
}

inline double rand_gamma_big(double a)
{
	// Cheng (1977) The generation of gamma variables with non-integral
	// shape parameter. Appl. Stat. 26(1): 71-75.     
	double aa = sqrt(a+a+1.0);
	double b = a - 1.38629436111989061883;
	double c = a + aa;
	double u1,u2,v,x,z,r;
	do
	{
		do {u1 = rand_real(); } while(u1 == 0.0);
		u2 = rand_real();
		v = log(u1/(1.0-u1))/aa;
		x = a*exp(v);
		z = u1*u1*u2;
		r = b+c*v-x;
	} while( r+2.50407739677627407337 < 4.5*z && r < log(z));
	return x;
}

// Draw from Gamma with mean 'a' and var 'a'
inline double rand_gamma(double a)
{
	if(a<=0.0)
		return 0.0;
	else if(a<1.0)
		return rand_gamma_small(a);
	else if(a>1.0)
		return rand_gamma_big(a);
	else
		return rand_exp();
}

// Draw from Gamma with mean 'ab' and var 'abb'
inline double rand_gamma(double a, double b)
{
	return b*rand_gamma(a);
}

// Draw from Gamma with mean '1' and var 'b'
inline double rand_gamma1(double b)
{
	return rand_gamma(1.0/b, b);
}

// Draw from Geometric(1-q):
//   P(X=x)=(1-q)q^x; x <=0
inline uint32_t rand_geometric(double q)
{
	return (uint32_t)(log(1.0-rand_real())/log(q));
}

// Draw from Negative Binomial(r, 1-q):
//   P(X=x) = (r+x-1 nch x)q^x(1-q)^r; x>=0
inline uint32_t rand_negbinomial(uint32_t r, double q)
{
	uint32_t u=0;
	while(r--)
		u += rand_geometric(q);
	return u;
}

// Draw Poisson(l)
inline uint32_t rand_poisson(double lambda)
{
	uint32_t u = 0;
	double d, e = exp(-lambda);
	// effecient for lambda < 12
	d = rand_real();
	while(d > e) {++u; d*=rand_real();}
	return u;
}

// Draw from Zipf distribution, with parameter a > 1.0
// Devroye Luc (1986) Non-uniform random variate generation.
//     Springer-Verlag: Berlin. p551
inline uint32_t rand_zipf(double a)
{
	double b = pow(2.0, a-1.0);
	double x,t;
	do {
	 x = floor(pow(rand_real(), -1.0/(a-1.0)));
	 t = pow(1.0+1.0/x, a-1.0);
	} while( rand_real()*x*(t-1.0)*b > t*(b-1.0));
	return (uint32_t)x;
}

#endif
