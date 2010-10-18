#pragma once
#ifndef DAWG_MUTT_H
#define DAWG_MUTT_H
/****************************************************************************
 *  Copyright (C) 2009 Reed A. Cartwright, PhD <reed@scit.us>               *
 ****************************************************************************/

#include <dawg/details/mutt.h>

#include <cmath>
#include <cfloat>
#include <limits>
#include <ctime>

#include <vector>

#include <boost/config.hpp>
#include <boost/functional/hash.hpp>

#ifdef BOOST_WINDOWS
#	include <process.h>
#endif

namespace dawg {

class mutt {
public:
	typedef RANDOM_GEN_CLASS generator;
	typedef generator::native_t uint_t;

	inline void seed(boost::uint32_t s) { 
		_seed.assign(1, s);
		gen.seed(_seed.begin(), _seed.end());
	}
	
	template<typename _It>
	inline void seed(_It first, _It last) {
		_seed.assign(first, last);
		if(_seed.size() == 1)
			gen.seed(_seed[0]);
		else if(_seed.size() > 1)
			gen.seed(first, last);
	}
	
	// returns a random double between (0,1)
	//inline double operator()() { return rand_real(); }
	inline double rand_real() { return gen.rand_real(); }
	// between (0,1] and (0,1)
	// inline double rand_real_oc() { return gen.rand_real_oc(); }
	// inline double rand_real_oo() { return gen.rand_real_oo(); }
	// returns random 32 or 64-bit number
	uint_t rand_uint() { return gen.rand_native(); }
	// returns random 64-bit number
	boost::uint64_t rand_uint64() { return gen.rand_uint64(); }
	// returns random 32-bit number
	boost::uint32_t rand_uint32() { return gen.rand_uint32(); }
	// returns random 32-bit number with max n-1
	boost::uint32_t rand_uint32(boost::uint32_t n) {
		return static_cast<boost::uint32_t>(rand_real()*n);
	}
	// random boolean p prob of success
	inline bool rand_bool(double p = 0.5) { return (rand_real() < p); }
	// random exponential with rate r
	inline double rand_exp(double r = 1.0) { return rand_exp_zig()/r; }
	inline double rand_exp_zig();
	inline double rand_exp_inv();
	
	// random geometric with succcess rate p in (0,1)
	// p(1 - p)^(k-1) for k >= 1; q = -log(1-p)
	// mean is 1/p
	
	inline boost::uint32_t rand_geometric_q(double q) {
		return 1+static_cast<boost::uint32_t>(rand_exp(q));
	}

	// random zeta with slope z
	boost::uint32_t rand_zeta(double z);
	
	// random gamma with mean ab and var abb
	inline double rand_gamma(double a, double b) {
		if(a < 1)
			return rand_gamma_low(a,b);
		else
			return rand_gamma_high(a,b);
	}
	// assumes a >= 1
	double rand_gamma_high(double a, double b);
	// assumes a < 1
	inline double rand_gamma_low(double a, double b) {
		// TODO: check to see if this can be _oc
		return rand_gamma_high(1.0+a,b)*pow(rand_real(), 1.0/a);
	}
	
	// random normal with mean and sigma
	inline double rand_normal(double m, double s) {
		return rand_normal(s)+m;
	}
	// random normal with sigma
	double rand_normal(double sigma);
	
	mutt();
	
private:
	generator gen;
	std::vector<boost::uint32_t> _seed;
	
	static const boost::uint64_t ek[256];
	static const double ew[256], ef[256];
	
	static inline double zH(double x, double z1, double z2) {
		return pow(x+1.0, z1)*z2;
	}
	static inline double zHi(double x, double z1, double z2) {
		return pow(z1*x,z2)-1.0;
	}	
};

// George Marsaglia's quick but good generator
// Short Period: 2^32-1

inline boost::uint32_t create_random_seed() {
	std::size_t v = static_cast<std::size_t>(getpid());
	v += (v << 15) + (v >> 3); // Spread 5-decimal PID over 32-bit number
	boost::hash_combine(v, time(NULL));
	// finish with one round of shr3
	return static_cast<boost::uint32_t>((v^=(v<<17), v^=(v>>13), v^=(v<<5)));
}

inline mutt::mutt() {
	seed(create_random_seed());
}

inline double mutt::rand_exp_zig() {
	static const double r = 7.69711747013104972;
	boost::uint64_t a = rand_uint64();
	int b = a & 255;
	if( a <= ek[b])
		return a*ew[b];
	do {
		if(b == 0)
			return r+rand_exp_inv();
		double x = a*ew[b];
		// we can cache ef[b-1]-ef[b], but it should be minor
		if(ef[b]+rand_real()*(ef[b-1]-ef[b]) < exp(-x) )
			return x;
		a = rand_uint64();
		b = a & 255;
	} while(a > ek[b]);
	
	return a*ew[b];
}

inline double mutt::rand_exp_inv() {
	return -log(rand_real());
}

} // namespace dawg
#endif
