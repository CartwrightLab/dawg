#pragma once
#ifndef DAWG_DETAILS_MUTT_H
#define DAWG_DETAILS_MUTT_H
/****************************************************************************
 *  Copyright (C) 2009-2012 Reed A. Cartwright, PhD <reed@scit.us>          *
 ****************************************************************************/

#ifndef __STDC_CONSTANT_MACROS
#	define __STDC_CONSTANT_MACROS 1
#endif
#ifndef __STDC_LIMIT_MACROS
#	define __STDC_LIMIT_MACROS 1
#endif

#include <dawg/details/config.h>

#include <boost/cstdint.hpp>
#include <boost/functional/hash.hpp>
#include <cmath>
#include <cfloat>

#ifndef RANDOM_GEN_HEADER
#ifndef _NDEBUG
#	define RANDOM_GEN_HEADER <dawg/details/xorshift_64.h>
#elif
#	error RANDOM_GEN_HEADER is not defined.
#endif
#endif

namespace dawg { namespace details {

inline boost::uint64_t to_uint64(boost::uint32_t x, boost::uint32_t y) {
	return y | ((boost::uint64_t)x << 32);
}

/* generate a random number on (0,1)-real-interval with 32-bit precision */
inline double to_double32(boost::uint32_t v) {
	boost::int32_t x = static_cast<boost::int32_t>(v); // + 2147483648;
	return (0.5+0.5/4294967296.0) + (x*(1.0/4294967296.0));
}

/* generate a random number on (0,1) with 52-bit resolution*/
/* my idea, seems nearly unique */
inline double to_double52(boost::uint64_t v) { 
	union {	boost::uint64_t u;	double d; } a;
	a.u = (v >> 12) | UINT64_C(0x3FF0000000000000);
	double q = (1.0-(DBL_EPSILON/2.0));
	return a.d-q;
}
inline double to_double52(boost::uint32_t x, boost::uint32_t y) { 
    return to_double52(to_uint64(x,y));
}

/* generate a random number on [0,1) with 53-bit resolution*/
/* my idea */
inline double to_real53(boost::uint64_t v) {
	union {	boost::uint64_t u;	double d; } a;
	a.u = (v >> 12) | UINT64_C(0x3FF0000000000000);
	double b = (v&2048) ?  (1.0-(DBL_EPSILON/2.0)) : 1.0;
	return a.d-b;
}
inline double to_real53(boost::uint32_t x, boost::uint32_t y) { 
    return to_real53(to_uint64(x,y));
}

}} /* namespace dawg::details */

#include RANDOM_GEN_HEADER

#endif /* DAWG_DETAILS_MUTT_H */
