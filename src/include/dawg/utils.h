#pragma once
#ifndef DAWG_UTILS_H
#define DAWG_UTILS_H
/****************************************************************************
 *  Copyright (C) 2009 Reed A. Cartwright, PhD <reed@scit.us>               *
 ****************************************************************************/

#include <boost/algorithm/string/predicate.hpp>
#include <boost/cstdint.hpp>
 
namespace dawg {

template<class A, class B, std::size_t _N>
std::size_t key_switch(A &ss, const B (&key)[_N]) {
	using boost::algorithm::istarts_with;
	for(std::size_t i=0;i<_N;++i) {
		if(istarts_with(key[i], ss))
			return i;
	}
	return (std::size_t)-1;
}

template<typename It, typename _V>
inline std::size_t search_binary_cont(It first, It last, const _V &v) {
	std::size_t r = 0;
	for(std::size_t u = (last-first)/2; u > 0; u /= 2) {
		if(v >= first[r+u-1])
			r += u;
	}
	return r;
}

template<typename _V, std::size_t _N>
inline std::size_t search_binary_cont(_V (&a)[_N], const _V &v) {
	return search_binary_cont(&a[0], &a[_N], v);
}

template<typename T>
inline T upper_binary(T u) {
	u--;
	for(int i=1;i<8*sizeof(u);i*=2)
		u |= u >> i;
	u++;
	return u;
}

} /* namespace dawg */
 
#endif /* DAWG_UTILS_H */
 
