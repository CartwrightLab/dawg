#pragma once
#ifndef DAWG_DETAILS_MUTT_GEN_H
#define DAWG_DETAILS_MUTT_GEN_H
/****************************************************************************
 *  Copyright (C) 2009-2012 Reed A. Cartwright, PhD <reed@scit.us>          *
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
A 32-bit analog of xorshift_64 with nealry the same period

R. P. Brent's high period xorshift generator
Reference implementation:
http://wwwmaths.anu.edu.au/~brent/ftp/random/xorgens305.tar.gz
Uses the 2x32-bit parameters from http://wwwmaths.anu.edu.au/~brent/ftp/random/xortable.txt
*/

struct xorshift_32_mutt_gen {
	typedef boost::uint32_t native_t;
	struct _state_t {
		native_t x,y,w;
	};
	typedef _state_t state_t;

	xorshift_32_mutt_gen() {
		seed(0);
	}
	
	inline native_t rand_native() { return rand_uint32(); }

	inline boost::uint64_t rand_uint64() {
		native_t a = rand_uint32();
		native t b = rand_uint32();
		return to_uint32(a,b);
	}
	inline native_t rand_uint32() {
		x ^= x << 17; x ^= x >> 14;
		y ^= y << 12; y ^= y >> 19;
		y ^= x; x ^= y;
#ifndef DAWG_DISABLE_WEYL_GENERATOR			
		w += 0x61c88647;
		return y+(w^(w>>16));
#else
		return y;
#endif
	}
	// doubles with 52-bits worth of precision
	inline double rand_real()   { return to_double52(rand_uint64()); }

	inline state_t state() const {
		state_t r;
		r.x = x; r.y = y; r.w = w;
		return r;
	}
	inline void state(const state_t z) {
		x = x.x; y = z.y; w = z.w;
	}

	inline void seed(boost::uint32_t xx) {
		x = 0x7F295846; y = 0x6A7BCC42; w = 0x61c88647;
		if(xx != 0) {
			rand_native();
			x ^= static_cast<native_t>(xx);	
		}
		for(int i=0;i<128;++i)
			rand_native();
	}
	template<int _N>
	inline void seed(boost::uint32_t (&xx)[_N]) {
		seed(&xx[0],&xx[_N]);
	}
	template<typename _It>
	inline void seed(_It first, _It last) {
		x = 0x7F295846; y = 0x6A7BCC42; w = 0x61c88647;
		for(;first != last;++first) {
			rand_native();
			x ^= static_cast<native_t>(*first);			
		}
		for(int i=0;i<128;++i)
			rand_native();
	}
	
	private:
		native_t x,y,w;
};

typedef xorshift_32_mutt_gen mutt_gen_default;

}} //namespace dawg::details
#endif
