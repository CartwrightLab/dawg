#pragma once
#ifndef DAWG_DETAILS_MUTT_H
#define DAWG_DETAILS_MUTT_H
/****************************************************************************
 *  Copyright (C) 2009 Reed A. Cartwright, PhD <reed@scit.us>               *
 ****************************************************************************/

#define __STDC_CONSTANT_MACROS 1
#include <boost/cstdint.hpp>
#include <boost/functional/hash.hpp>
#include <cmath>

// include stuff for DSFMT
#if defined(_MSC_VER) || defined(__BORLANDC__)
#	define DSFMT_UINT32_DEFINED 1
using boost::uint32_t;
using boost::uint64_t;
#endif
#define DSFMT_MEXP 19937
#define DSFMT_DO_NOT_USE_OLD_NAMES
extern "C" {
#include <dawg/details/dSFMT.h>
}

namespace dawg { namespace details {

struct dsfmt_mutt_gen {
	double rand_01()  { return dsfmt_genrand_close_open(&state);	}
	double rand_01oc() { return dsfmt_genrand_open_close(&state); }
	double rand_01oo() { return dsfmt_genrand_open_open(&state); }
	boost::uint32_t rand_uint32() { return dsfmt_genrand_uint32(&state); }
	boost::uint64_t rand_uint64() { 
		uint64_t a = rand_uint32();
		uint64_t b = rand_uint32();
		return (a << 32) | b;
	}
	void seed(uint32_t x) { dsfmt_init_gen_rand(&state, x); }
	template<int _N>
	void seed(uint32_t (&x)[_N]) {
		dsfmt_init_by_array(&state, &x[0], _N);
	}
	template<typename _It>
	void seed(_It first, _It last) {
		std::size_t sz = last-first;
		boost::uint32_t *p = new boost::uint32_t[sz];
		std::copy(first, last, p);
		dsfmt_init_by_array(&state, p, sz);
		delete[] p;
	}

private:
	dsfmt_t state;
};

// George Marsaglia's quick but good generator
// http://www.jstatsoft.org/v08/i14/paper
// 64-bit period
// 0,0 never seen
struct shr3a_mutt_gen {
	inline boost::uint32_t rand_uint32() {
		uint32_t t = x^(x<<2);
		x=y;
		y=(y^(y>>3))^(t^(t>>7));
		return y;
	}
	inline boost::uint64_t rand_uint64() {
		uint64_t a = rand_uint32();
		uint64_t b = rand_uint32();
		return (a << 32) | b;
	}
	// doubles with 53-bits worth of precision
	inline double rand_01() {
		uint64_t u = rand_uint64() & UINT64_C(4503599627370495);
		return u/4503599627370496.0;
	}
	inline double rand_01oo() {
		uint64_t u = rand_uint64() & UINT64_C(4503599627370495);
		return ((u+1)/4503599627370497.0);
	}
	inline double rand_01oc() {
		uint64_t u = rand_uint64() & UINT64_C(4503599627370495);
		return ((u+1)/4503599627370496.0);
	}
	inline void seed(boost::uint32_t xx) { 
		y = xx;
		// 1 round of 32-bit shr3 to fill in the rest
		x = ((xx^=(xx<<17), xx^=(xx>>13), xx^=(xx<<5)));
	}
	template<int _N>
	inline void seed(boost::uint32_t (&xx)[_N]) {
		seed(&xx[0],&xx[_N]);
	}
	template<typename _It>
	inline void seed(_It first, _It last) {
		if(first == last)
			return; // nothing to do
		y = *first++;
		if(first == last) {
			seed(y);  // single value
			return;
		}
		// use hash to make every element matter
		std::size_t h = static_cast<std::size_t>(*first++);
		boost::hash_range(h, first, last);
		x = static_cast<boost::uint32_t>(h);
	}	
	
	private:
		boost::uint32_t x,y;
};

// 0 never seen
struct shr3b_mutt_gen {
	inline boost::uint32_t rand_uint32() {
		return static_cast<boost::uint32_t>(rand_uint64());
	}
	inline boost::uint64_t rand_uint64() {
		y ^= (y << 13);
		y ^= (y >> 7);
		y ^= (y << 17);
		return y;
	}
	// doubles with 53-bits worth of precision
	inline double rand_01() {
		uint64_t u = rand_uint64() & UINT64_C(4503599627370495);
		return u/4503599627370496.0;
	}
	inline double rand_01oo() {
		uint64_t u = rand_uint64() & UINT64_C(4503599627370495);
		return ((u+1)/4503599627370497.0);
	}
	inline double rand_01oc() {
		uint64_t u = rand_uint64() & UINT64_C(4503599627370495);
		return ((u+1)/4503599627370496.0);
	}
	inline void seed(boost::uint32_t xx) { 
		y = xx;
	}
	template<int _N>
	inline void seed(boost::uint32_t (&xx)[_N]) {
		seed(&xx[0],&xx[_N]);
	}
	template<typename _It>
	inline void seed(_It first, _It last) {
		if(first == last)
			return; // nothing to do
		y = *first++;
		if(first == last)
			return;
		// use hash to make every element matter
		std::size_t h = static_cast<std::size_t>(*first++);
		boost::hash_range(h, first, last);
		y |= (static_cast<boost::uint64_t>(h) << 32);
		
	}
	
	private:
		boost::uint64_t y;
};

// Ziggurat method for rand exponential value
// Memory hog, but really really fast
class rand_exp_ziggurat {
public:

	template<class _R>	
	inline double operator()(_R &rng) const {
		static const double r = 7.69711747013104972;
		boost::uint64_t a = rng.rand_uint64();
		int b = a & 255;
		a >>= 12;
		if( a <= k[b])
			return a*w[b];
		do {
			if(b == 0)
				return r-log(unif01open(rng.rand_uint64()));
			double x = a*w[b];
			// we can cache f[b-1]-f[b], but it should be minor
			if(f[b]+unif01(rng.rand_uint64())*(f[b-1]-f[b]) < exp(-x) )
				return x;
			a = rng.rand_uint64();
			b = a & 255;
			a >>= 12;
		} while(a > k[b]);
		
		return a*w[b];
	}
	template<class _R>
	inline double inv(_R &rng) const {
		return log(unif01open(rng.rand_uint64()));
	}
	rand_exp_ziggurat() {
		static const double v = 0.0039496598225815571993;
		static const double m = 4503599627370496.0; // 2^52
		
		double d=7.69711747013104972, t=d, q;
		
		q = v/exp(-d);
		k[0]=(d/q)*m;
		k[1]=0;

		w[0]=q/m;
		w[255]=d/m;

		f[0]=1;
		f[255]=q=exp(-d);

		for(int i=254; i>=1; i--) {
			d=-log(v/d+q);
			k[i+1] = (d/t)*m;
			t=d;
			f[i]=q=exp(-d);
			w[i]=d/m;
		}
	}
	
protected:
	static double inline unif01open(boost::uint64_t u) {
		u &= UINT64_C(4503599627370495);
		return ((u+1)/4503599627370497.0);
	}

	static double inline unif01(boost::uint64_t u) {
		u &= UINT64_C(4503599627370495);
		return u/4503599627370496.0;
	}

	double w[256], f[256];
	boost::uint64_t k[256];
};

}} /* namespace dawg::details */

#endif /* DAWG_DETAILS_MUTT_H */

