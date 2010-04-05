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

#include <vector>

#include <boost/config.hpp>
#include <boost/cstdint.hpp>
#include <boost/functional/hash.hpp>

#include <dawg/details/mutt.h>

#ifdef BOOST_WINDOWS
#	include <process.h>
#endif

namespace dawg {

class mutt {
public:
	typedef details::dsfmt_mutt_gen generator;
	
	inline void seed(uint32_t s) { 
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
	
	// returns a random double between [0,1)
	inline double operator()() { return rand_01(); }
	inline double rand_01() { return gen.rand_01(); }
	// between (0,1] and (0,1)
	inline double rand_01oc() { return gen.rand_01oc(); }
	inline double rand_01oo() { return gen.rand_01oo(); }
	// returns random 32-bit number
	boost::uint32_t rand_uint32() { return gen.rand_uint32(); }
	// returns random 32-bit number with max n-1
	boost::uint32_t rand_uint32(boost::uint32_t n) {
		return static_cast<boost::uint32_t>(rand_01()*n);
	}
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
		return rand_gamma_high(1.0+a,b)*pow(rand_01oo(), 1.0/a);
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
	std::vector<uint32_t> _seed;
	
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
	shr3(boost::uint32_t s) : jsr(s) { }
	
	inline boost::uint32_t rand_uint32() {
		return (jsr^=(jsr<<17), jsr^=(jsr>>13), jsr^=(jsr<<5));
	}
	inline void seed(uint32_t s) {
		jsr = s;
	}
private:
	boost::uint32_t jsr;
};

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

} // namespace dawg
#endif

