#pragma once
#ifndef DAWG_DETAILS_SHR3B_H
#define DAWG_DETAILS_SHR3B_H
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
	double rand_real()   { return to_real52_oo(rand_uint64()); }

	inline void seed(boost::uint32_t xx) { y = xx; }
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

}} //namespace dawg::details
#endif
 