#pragma once
#ifndef DAWG_DETAILS_MUTT_H
#define DAWG_DETAILS_MUTT_H
/****************************************************************************
 *  Copyright (C) 2009-2010 Reed A. Cartwright, PhD <reed@scit.us>          *
 ****************************************************************************/

#ifndef __STDC_CONSTANT_MACROS
#	define __STDC_CONSTANT_MACROS 1
#endif
#include <boost/cstdint.hpp>
#include <boost/functional/hash.hpp>
#include <cmath>
#include <cfloat>

#ifndef RANDOM_GEN_HEADER
#	define RANDOM_GEN_HEADER <dawg/details/shr3a.h>
#endif
#ifndef RANDOM_GEN_CLASS
#	define RANDOM_GEN_CLASS dawg::details::shr3a_mutt_gen
#endif

namespace dawg { namespace details {

inline boost::uint64_t to_uint64(boost::uint32_t x, boost::uint32_t y) {
	return y | ((boost::uint64_t)x << 32);
}

/* These real versions are derived from Isaku Wada's code. */
/* generate a random number on [0,1)-real-interval */
inline double to_real_co(boost::uint32_t v) {
    return v * (1.0/4294967296.0); /* divided by 2^32 */
}

/* generate a random number on (0,1)-real-interval */
inline double to_real_oo(boost::uint32_t v) {
    return (((double)v) + 0.5)*(1.0/4294967296.0); /* divided by 2^32 */
}

/* generate a random number on (0,1]-real-interval */
inline double to_real_oc(boost::uint32_t v) {
        return 1.0 - (v * (1.0/4294967296.0)); /* divided by 2^32 */
}

/* generate a random number on [0,1) with 52-bit resolution*/
inline double to_real52_co(boost::uint64_t v) { 
	union {	boost::uint64_t u;	double d; } a;
	a.u = (v >> 12) | UINT64_C(0x3FF0000000000000);
	return a.d-1.0;
}

/* generate a random number on [0,1) with 52-bit resolution from two
 * 32 bit integers */
inline double to_real52_co(boost::uint32_t x, boost::uint32_t y) { 
    return to_real52_co(to_uint64(x,y));
}

/* generate a random number on (0,1] with 52-bit resolution*/
inline double to_real52_oc(boost::uint64_t v) { 
	union {	boost::uint64_t u;	double d; } a;
	a.u = (v >> 12) | UINT64_C(0x3FF0000000000000);
	return 2.0-a.d;
}

/* generate a random number on (0,1] with 52-bit resolution from two
 * 32 bit integers */
inline double to_real52_oc(boost::uint32_t x, boost::uint32_t y) { 
    return to_real52_oc(to_uint64(x,y));
}

/* generate a random number on (0,1) with 52-bit resolution*/
/* my idea */
/* uses the high 52-bits, but could use low 52-bits as well */
inline double to_real52_oo(boost::uint64_t v) {
	union {	boost::uint64_t u;	double d; } a;
	a.u = (v >> 12) | UINT64_C(0x3FF0000000000000);
	return a.d-(1.0-(DBL_EPSILON/2.0));
}

/* generate a random number on (0,1) with 52-bit resolution from two
 * 32 bit integers */
inline double to_real52_oo(boost::uint32_t x, boost::uint32_t y) { 
    return to_real52_oo(to_uint64(x,y));
}

/* generate a random number on [0,1) with 53-bit resolution*/
/* my idea */
/* uses the high 53-bits, but could use low 53-bits as well */
inline double to_real53_co(boost::uint64_t v) {
	union {	boost::uint64_t u;	double d; } a;
	double b = (v&2048) ?  (1.0-(DBL_EPSILON/2.0)) : 1.0;
	a.u = (v >> 12) | UINT64_C(0x3FF0000000000000);
	return a.d-b;
}

/* generate a random number on [0,1) with 53-bit resolution from two
 * 32 bit integers */
inline double to_real53_co(boost::uint32_t x, boost::uint32_t y) { 
    return to_real53_co(to_uint64(x,y));
}

}} /* namespace dawg::details */

#include RANDOM_GEN_HEADER

#endif /* DAWG_DETAILS_MUTT_H */
