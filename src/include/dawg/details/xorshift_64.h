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

// Based on George Marsaglia's Xorshift Generator.
// Includes ideas from Richard Brent and Francois Panneton
// http://www.jstatsoft.org/v08/i14/paper
// http://wwwmaths.anu.edu.au/~brent/random.html
// http://www.iro.umontreal.ca/~panneton/These.pdf
//
// Passes the BigCrush tests

struct xorshift_64_mutt_gen {
	typedef boost::uint64_t native_t;
	typedef std::pair<native_t,native_t> state_t;
	xorshift_64_mutt_gen() {
		seed(0);
	}
	
	inline native_t rand_native() { return rand_uint64(); }

	inline boost::uint32_t rand_uint32() {
		return static_cast<boost::uint32_t>(rand_uint64() >> 32);
	}
	inline native_t rand_uint64() {
		y ^= (y << 5); y ^= (y >> 15); y ^= (y << 27);
#ifndef DAWG_DISABLE_WEYL_GENERATOR			
		w += UINT64_C(0x61C8864680B583EB);
		return y+(w^(w>>27));
#else
		return y;
#endif
	}
	// doubles with 52-bits worth of precision
	inline double rand_real()   { return to_double52(rand_uint64()); }

	inline state_t state() const {
		return std::make_pair(y,w);
	}
	inline void state(const state_t x) {
		y = x.first;
		w = x.second;
	}

	inline void seed(boost::uint32_t xx) {
		y = UINT64_C(0x6A7BCC427F295846);
		w = UINT64_C(0x61C8864680B583EB);
		if(xx != 0) {
			rand_native();
			y ^= static_cast<native_t>(xx);	
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
		y = UINT64_C(0x6A7BCC427F295846);
		w = UINT64_C(0x61C8864680B583EB);
		for(;first != last;++first) {
			rand_native();
			y ^= static_cast<native_t>(*first);			
		}
		for(int i=0;i<128;++i)
			rand_native();
	}
	
	private:
		native_t y,w;
};

typedef xorshift_64_mutt_gen mutt_gen_default;

}} //namespace dawg::details
#endif
