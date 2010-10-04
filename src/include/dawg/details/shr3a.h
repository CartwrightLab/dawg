#pragma once
#ifndef DAWG_DETAILS_SHR3A_H
#define DAWG_DETAILS_SHR3A_H
/****************************************************************************
 *  Copyright (C) 2009-2010 Reed A. Cartwright, PhD <reed@scit.us>          *
 ****************************************************************************/

#ifndef __STDC_CONSTANT_MACROS
#	define __STDC_CONSTANT_MACROS 1
#endif
#include <boost/cstdint.hpp>

namespace dawg { namespace details {

// George Marsaglia's quick and good generator
// http://www.jstatsoft.org/v08/i14/paper
// 64-bit period
// 0,0 never seen
struct shr3a_mutt_gen {
	inline boost::uint32_t rand_uint32() {
		boost::uint32_t t = x^(x<<2);
		x=y;
		y=(y^(y>>3))^(t^(t>>7));
		return y;
	}
	boost::uint64_t rand_uint64() { return to_uint64(rand_uint32(),rand_uint32()); }
	double rand_real()   { return to_real52_oo(rand_uint64()); }

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

}} //namespace dawg::details
#endif
