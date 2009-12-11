#pragma once
#ifndef DAWG_MUTT_H
#define DAWG_MUTT_H
/****************************************************************************
 *  Copyright (C) 2009 Reed A. Cartwright, PhD <reed@scit.us>               *
 ****************************************************************************/

#include <cmath>
#include <cfloat>
#include <limits>

#include <dawg/details/mutt.h>

namespace dawg {

class mutt {
public:
	typedef details::dsfmt_mutt_gen generator;
	
	inline void seed(uint32_t s) { gen.seed(s); }
	
	// returns a random double between [0,1)
	inline double operator()() { return rand_01(); }
	inline double rand_01() { return gen.rand_01(); }
	// between (0,1] and (0,1)
	inline double rand_01oc() { return gen.rand_01oc(); }
	inline double rand_01oo() { return gen.rand_01oo(); }
	// returns random 32-bit number
	uint32_t rand_uint32() { return gen.rand_uint32(); }
	// random boolean p prob of success
	bool rand_bool(double p = 0.5) { return (rand_01() < p); }
	// random exponential with rate r
	double rand_exp(double r = 1.0) { return -log(rand_01oc())/r; }
	
	// random geometric with succcess rate p in (0,1)
	// p(1 - p)^(k-1) for k >= 1
	inline uint32_t rand_geometric(double p) {
		return 1+static_cast<uint32_t>(rand_exp(-log(1.0-p)));
	}
	
	// random zeta distribution with slope z
	// rejection-inversion method of H\"ormann and Derflinger (1996)
	// idea borrowed from Indelible
	// optimizations for usage in Dawg
	inline uint32_t rand_zeta(double z) {
		double z1 = 1.0-z;
		double z2 = 1.0/z1;
		double s = 2.0-zHi(zH(1.5,z1,z2)-pow(2.0,-z),z1,z2);
		double Him = zH(4294967296.5,z1,z2);
		double Hx0 = zH(0.5,z1,z2)-1.0-Him;
		
		double U, X,K;
		for(;;) {
			U = rand_01(); // [0,1)
			U = Him+U*Hx0;
			X = zHi(U,z1,z2);
			K = floor(X+1.5);
			if(K-X < s)
				break;
			if(U > zH(K-0.5,z1,z2)-pow(K,-z))
				break;
		}
		return static_cast<uint32_t>(K);
	}
	
private:
	generator gen;
	
	inline double zH(double x, double z1, double z2) {
		return pow(x+1.0, z1)*z2;
	}
	inline double zHi(double x, double z1, double z2) {
		return pow(z1*x,z2)-1.0;
	}
};

} // namespace dawg
#endif

