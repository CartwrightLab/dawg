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

struct xorgen4096_32_mutt_gen {
	typedef boost::uint32_t native_t;
	struct _state_t {
		native_t y[128];
		native_t w;
		std::size_t i;
	};
	typedef _state_t state_t;

	xorgen4096_32_mutt_gen() {
		seed(0);
	}
	
	inline native_t rand_native() { return rand_uint32(); }

	inline boost::uint64_t rand_uint64() {
		native_t a = rand_uint32();
		native_t b = rand_uint32();
		return to_uint64(a,b);
	}
		
	inline native_t rand_uint32() {
		native_t t = y[i&127];
		native_t v = y[(i+33)&127];
		t ^= t << 17; t ^= t >> 12;
		v ^= v << 13; v ^= v >> 15; v ^= t;
		y[i&127] = v;
		i += 1;
		w += 0x61c88647;
		return v + (w^(w>>16));
	}
	
	inline double rand_real()   { return to_double52(rand_uint64()); }
	
	inline state_t state() const {
		state_t tmp;
		std::copy(&y[0], &y[128], &tmp.y[0]);
		tmp.w = w;
		tmp.i = i&127;
		return tmp;
	}
	inline void state(const state_t& s) {
		std::copy(&s.y[0], &s.y[128], &y[0]);
		w = s.w;
		i = s.i&127;
	}
	
	inline void seed(boost::uint32_t xx) {
		i = 0;
		w = 0x61c88647;
		native_t x = 0x6A7BCC42;
		// Fill array using a 32-bit xorshift
		for(int a=0;a<64;++a) {
			x ^= (x >> 7); x ^= (x << 1); x ^= (x >> 9);
			w += 0x61c88647;
			y[a] = x + (w^(w>>16));
		}
		// Merge seed into array
		if(xx != 0) {
			x = static_cast<native_t>(xx);
			for(int a=0;a<64;++a) {
				x ^= (x >> 7); x ^= (x << 1); x ^= (x >> 9);
				w += 0x61c88647;
				y[a] ^= x + (w^(w>>16));
			}
		}
		// Discard initial results; four passes though y
		for(int a=0;a<512;++a)
			rand_native();
	}
	
	template<int _N>
	inline void seed(boost::uint32_t (&xx)[_N]) {
		seed(&xx[0],&xx[_N]);
	}
	template<typename _It>
	inline void seed(_It first, _It last) {
		i = 0;
		w = 0x61c88647;
		native_t x = 0x6A7BCC42;
		// Fill array using xorshift_64's algorithm
		for(int a=0;a<64;++a) {
			x ^= (x >> 7); x ^= (x << 1); x ^= (x >> 9);
			w += 0x61c88647;
			y[a] = x + (w^(w>>16));
		}
		// Merge seeds into array
		for(;first != last;++first) {
			x = static_cast<native_t>(*first);
			if(x == 0)
				continue;
			for(int a=0;a<64;++a) {
				x ^= (x >> 7); x ^= (x << 1); x ^= (x >> 9);
				w += 0x61c88647;
				y[a] ^= x + (w^(w>>16));
			}
		}
		// Discard initial results; four passes though y
		for(int a=0;a<512;++a)
			rand_native();
	}
	
private:
	native_t y[128];
	native_t w;
	std::size_t i;
};

typedef xorgen4096_32_mutt_gen mutt_gen_default;

}} //namespace dawg::details
#endif
