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
	xorshift_64_mutt_gen() : y(UINT64_C(15191868757011070976))
#ifndef DAWG_DISABLE_WEYL_GENERATOR		
	, w(UINT64_C(0x61C8864680B583EB))
#endif
	{ }
	
	inline native_t rand_native() { return rand_uint64(); }

	inline boost::uint32_t rand_uint32() {
		return static_cast<boost::uint32_t>(rand_uint64() >> 32);
	}
	inline boost::uint64_t rand_uint64() {
		y ^= (y << 5); y ^= (y >> 15); y ^= (y << 27);
#ifndef DAWG_DISABLE_WEYL_GENERATOR			
		w += UINT64_C(0x61C8864680B583EB);
		return y+w;
#else
		return y;
#endif
	}
	// doubles with 52-bits worth of precision
	double rand_real()   { return to_real52_oo(rand_uint64()); }

	inline void seed(boost::uint32_t xx) {
		y = (xx == 0) ? UINT64_C(15191868757011070976) : xx;
#ifndef DAWG_DISABLE_WEYL_GENERATOR			
		w = UINT64_C(0x61C8864680B583EB);
#endif
		rand_native();
	}
	template<int _N>
	inline void seed(boost::uint32_t (&xx)[_N]) {
		seed(&xx[0],&xx[_N]);
	}
	template<typename _It>
	inline void seed(_It first, _It last) {
#ifndef DAWG_DISABLE_WEYL_GENERATOR			
		w = UINT64_C(0x61C8864680B583EB);
#endif		
		// check to see if there is 0 or 1 element; pass to simple seed function
		y = (first == last) ? 0 : *first++;
		if(first == last) {
			seed(y);
			return;
		}
		// use hash to make every element matter
		std::size_t h = y;
		boost::hash_range(h, first, last);
		y = (static_cast<boost::uint64_t>(h));
		rand_native();
	}
	
	private:
		boost::uint64_t y;
#ifndef DAWG_DISABLE_WEYL_GENERATOR		
		boost::uint64_t w;
#endif
};

typedef xorshift_64_mutt_gen mutt_gen;

}} //namespace dawg::details
#endif

