#pragma once
#ifndef DAWG_DETAILS_MUTT_GEN_H
#define DAWG_DETAILS_MUTT_GEN_H
/****************************************************************************
 *  Copyright (C) 2012 Reed A. Cartwright, PhD <reed@scit.us>               *
 ****************************************************************************/

#ifndef __STDC_CONSTANT_MACROS
#	define __STDC_CONSTANT_MACROS 1
#endif
#ifndef __STDC_LIMIT_MACROS
#	define __STDC_LIMIT_MACROS 1
#endif
#include <boost/cstdint.hpp>
#include <utility>

namespace dawg { namespace details {

/*
R. P. Brent's high period xorshift generator
Reference implementation:
http://wwwmaths.anu.edu.au/~brent/ftp/random/xorgens305.tar.gz
*/

struct xorgen4096_64_mutt_gen {
	typedef boost::uint64_t native_t;
	struct _state_t {
		native_t y[64];
		native_t w;
		std::size_t i;
	};
	typedef _state_t state_t;

	inline native_t rand_native() { return rand_uint64(); }

	inline boost::uint32_t rand_uint32() {
		return static_cast<boost::uint32_t>(rand_uint64() >> 32);
	}
	
	inline native_t rand_uint64() {
		native_t t = y[i&63];
		native_t v = y[(i+11)&63];
		t ^= t << 33; t ^= t >> 26;
		v ^= v << 27; v ^= v >> 29; v ^= t;
		y[i] = v;
		i += 1;
		w += UINT64_C(0x61C8864680B583EB);
		return v + (w^(w>>27));
	}
	
	inline double rand_real()   { return to_double52(rand_uint64()); }
	
	inline state_t state() const {
		state_t tmp;
		std::copy(&y[0], &y[64], &tmp.y[0]);
		tmp.w = w;
		tmp.i = i&63;
		return tmp;
	}
	inline void state(const state_t& s) {
		std::copy(&s.y[0], &s.y[64], &y[0]);
		w = s.w;
		i = s.i&63;
	}
	
	inline void seed(boost::uint32_t xx) {
		i = 0;
		w = UINT64_C(0x61C8864680B583EB);
		native_t x = UINT64_C(15191868757011070976);
		// Fill array using xorshift_64's algorithm
		for(int a=0;a<64;++a) {
			x ^= (x << 5); x ^= (x >> 15); x ^= (x << 27);
			w += UINT64_C(0x61C8864680B583EB);
			y[a] = x + (w^(w>>27));
		}
		// Merge seed into array
		if(xx != 0) {
			x = static_cast<native_t>(xx);
			for(int a=0;a<64;++a) {
				x ^= (x << 5); x ^= (x >> 15); x ^= (x << 27);
				w += UINT64_C(0x61C8864680B583EB);
				y[a] ^= x + (w^(w>>27));
			}
		}
		// Discard initial results; four passes though y
		for(int a=0;a<256;++a)
			rand_native();
	}
	
	template<int _N>
	inline void seed(boost::uint32_t (&xx)[_N]) {
		seed(&xx[0],&xx[_N]);
	}
	template<typename _It>
	inline void seed(_It first, _It last) {
		i = 0;
		w = UINT64_C(0x61C8864680B583EB);
		native_t x = UINT64_C(15191868757011070976);
		// Fill array using xorshift_64's algorithm
		for(int a=0;a<64;++a) {
			x ^= (x << 5); x ^= (x >> 15); x ^= (x << 27);
			w += UINT64_C(0x61C8864680B583EB);
			y[a] = x + (w^(w>>27));
		}
		// Merge seeds into array
		for(;first != last;++first) {
			x = static_cast<native_t>(*first);
			if(x == 0)
				continue;
			for(int a=0;a<64;++a) {
				x ^= (x << 5); x ^= (x >> 15); x ^= (x << 27);
				w += UINT64_C(0x61C8864680B583EB);
				y[a] ^= x + (w^(w>>27));
			}
		}
		// Discard initial results; four passes though y
		for(int a=0;a<256;++a)
			rand_native();
	}
	
private:
	native_t y[64];
	native_t w;
	std::size_t i;
};

typedef xorgen4096_64_mutt_gen mutt_gen_default;

}} //namespace dawg::details
#endif
