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

extern "C" {
using boost::uint32_t;
using boost::uint64_t;
#define SFMT_UINT32_DEFINED 1
#include <dawg/details/SFMTx.h>
}


namespace dawg { namespace details {

struct sfmt_mutt_gen {
	typedef boost::uint32_t native_t;
	inline native_t rand_native() { return rand_uint32(); }
	boost::uint32_t rand_uint32() { return sfmt_gen_rand32(&state); }
	boost::uint64_t rand_uint64() { return to_uint64(rand_uint32(),rand_uint32()); }
	double rand_real()   { return to_real52_oo(rand_uint64()); }
	
	void seed(uint32_t x) { sfmt_init_gen_rand(&state, x); }
	template<int _N>
	void seed(uint32_t (&x)[_N]) {
		sfmt_init_by_array(&state, &x[0], _N);
	}
	template<typename _It>
	void seed(_It first, _It last) {
		std::size_t sz = last-first;
		boost::uint32_t *p = new boost::uint32_t[sz];
		std::copy(first, last, p);
		sfmt_init_by_array(&state, p, sz);
		delete[] p;
	}

private:
	sfmt_t state;
};

typedef sfmt_mutt_gen mutt_gen;

}} //namespace dawg::details
#endif
