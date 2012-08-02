#pragma once
#ifndef DAWG_UTILS_H
#define DAWG_UTILS_H
/****************************************************************************
 *  Copyright (C) 2009 Reed A. Cartwright, PhD <reed@scit.us>               *
 ****************************************************************************/

#include <boost/algorithm/string/predicate.hpp>
#include <boost/cstdint.hpp>
 
namespace dawg {

template<class A, class B, std::size_t N>
std::size_t key_switch(A &ss, const B (&key)[N]) {
	using boost::algorithm::istarts_with;
	for(std::size_t i=0;i<N;++i) {
		if(istarts_with(key[i], ss))
			return i;
	}
	return (std::size_t)-1;
}

template<typename It, typename V>
inline std::size_t search_binary_cont(It first, It last, const V &v) {
	std::size_t r = 0;
	for(std::size_t u = (last-first)/2; u > 0; u /= 2) {
		if(v >= first[r+u-1])
			r += u;
	}
	return r;
}

template<typename V, std::size_t N>
inline std::size_t search_binary_cont(V (&a)[N], const V &v) {
	return search_binary_cont(&a[0], &a[N], v);
}

template<typename T>
inline T upper_binary(T u) {
	u--;
	for(unsigned int i=1;i<8*sizeof(u);i*=2)
		u |= u >> i;
	u++;
	return u;
}

} /* namespace dawg */
 
#endif /* DAWG_UTILS_H */
 
