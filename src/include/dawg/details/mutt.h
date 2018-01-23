#pragma once
#ifndef DAWG_DETAILS_MUTT_H
#define DAWG_DETAILS_MUTT_H
/****************************************************************************
 *  Copyright (C) 2009-2018 Reed A. Cartwright, PhD <reed@scit.us>          *
 ****************************************************************************/

#ifndef DAWG_NO_CONFIG_H
#	include <dawg/details/config.h>
#endif

#include <cstdint>
#include <cmath>
#include <cfloat> // DBL_MANT_DIG, DBL_EPSILON

#ifndef RANDOM_GEN_HEADER
#ifndef _NDEBUG
#	define RANDOM_GEN_HEADER <dawg/details/xorshift_64.h>
#elif
#	error RANDOM_GEN_HEADER is not defined.
#endif
#endif

namespace dawg { namespace details {

inline std::uint64_t to_uint64(std::uint32_t x, std::uint32_t y) {
	return ((std::uint64_t)x << 32) | y;
}

/* adapter to generate a random number on (0,1) with 32-bit resolution */
inline double to_double32(std::uint32_t v) {
	std::int32_t x = static_cast<std::int32_t>(v); // + 2147483648;
	return (0.5+0.5/4294967296.0) + (x*(1.0/4294967296.0));
}

/* adapter to generate a random number on (0,1) with 52-bit resolution*/
inline double to_double52(std::uint64_t v) {
	union {	std::uint64_t u;	double d; } a;
	a.u = (v >> 12) | static_cast<std::uint64_t>(0x3FF0000000000000);
	double q = (1.0-(DBL_EPSILON/2.0));
	return a.d-q;
}
inline double to_double52(std::uint32_t x, std::uint32_t y) {
    return to_double52(to_uint64(x,y));
}

/* adapter to generate a random number on [0,1) with 53-bit resolution*/
inline double to_real53(std::uint64_t v) {
	union {	std::uint64_t u;	double d; } a,b;
	b.u = (v&2048)  ? static_cast<std::uint64_t>(0x3FEFFFFFFFFFFFFF)
	                : static_cast<std::uint64_t>(0x3FF0000000000000);
	a.u = (v >> 12) | static_cast<std::uint64_t>(0x3FF0000000000000);
	return a.d-b.d;
}
inline double to_real53(std::uint32_t x, std::uint32_t y) {
    return to_real53(to_uint64(x,y));
}

/* alternative algorithms; speed (dis)advantage depends on architecture. */
/* use more floating point math. */
inline double to_double52f(std::uint64_t v) {
	std::int64_t x = v >> (64-(DBL_MANT_DIG-1));
	return DBL_EPSILON*x + DBL_EPSILON/2.0;
}
inline double to_double53f(std::uint64_t v) {
	std::int64_t x = v >> (64-(DBL_MANT_DIG));
	return (DBL_EPSILON/2.0)*x;
}


}} /* namespace dawg::details */

#include RANDOM_GEN_HEADER

#endif /* DAWG_DETAILS_MUTT_H */
