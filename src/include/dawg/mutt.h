#pragma once
#ifndef DAWG_MUTT_H
#define DAWG_MUTT_H
/****************************************************************************
 *  Copyright (C) 2009 Reed A. Cartwright, PhD <reed@scit.us>               *
 ****************************************************************************/

#include <cmath>
#include <cfloat>
#include <limits>
#include <ctime>

#include <boost/cstdint.hpp>
#include <boost/functional/hash.hpp>

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
	boost::uint32_t rand_uint32() { return gen.rand_uint32(); }
	// random boolean p prob of success
	inline bool rand_bool(double p = 0.5) { return (rand_01() < p); }
	// random exponential with rate r
	inline double rand_exp(double r = 1.0) { return -log(rand_01oc())/r; }
	
	// random geometric with succcess rate p in (0,1)
	// p(1 - p)^(k-1) for k >= 1; q = 1-p
	inline boost::uint32_t rand_geometric(double q) {
		return 1+static_cast<uint32_t>(rand_exp(-log(q)));
	}

	// random zeta with slope z
	boost::uint32_t rand_zeta(double z);
	
	// random gamma with mean 1 and var b
	inline double rand_gamma(double b) {
		rand_gamma(1.0/b, b);
	}
	// random gamma with mean ab and var abb
	double rand_gamma(double a, double b);

	// random normal with mean and sigma
	inline double rand_normal(double m, double s) {
		return rand_normal(s)+m;
	}
	// random normal with sigma
	double rand_normal(double sigma);
	
private:
	generator gen;
	
	inline double zH(double x, double z1, double z2) {
		return pow(x+1.0, z1)*z2;
	}
	inline double zHi(double x, double z1, double z2) {
		return pow(z1*x,z2)-1.0;
	}
	
};

// George Marsaglia's quick but good generator
// Short Period: 2^32-1
class shr3 {
public:
	shr3() : jsr(1776) { }
	shr2(boost::uint32_t s) : jsr(s) { }
	
	inline boost::uint32_t rand_uint32() {
		return (jsr^=(jsr<<17), jsr^=(jsr>>13), jsr^=(jsr<<5));
	}
	inline void seed(uint32_t s) {
		jsr = s;
	}
private:
	boost::uinit_32 jsr;
};

inline boost::uint32_t create_random_seed() {
	boost::uint32_t v = static_cast<boost::uint32_t>(getpid());
	v += (v << 15) + (v >> 3); // Spread 5-decimal PID over 32-bit number
	boost::hash_combine(v, time(NULL));
	// finish with one round of shr3
	return (v^=(v<<17), v^=(v>>13), v^=(v<<5));	
}

} // namespace dawg
#endif

