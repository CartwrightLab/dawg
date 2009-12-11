#pragma once
#ifndef DAWG_UTILS_H
#define DAWG_UTILS_H
/****************************************************************************
 *  Copyright (C) 2009 Reed A. Cartwright, PhD <reed@scit.us>               *
 ****************************************************************************/

 #include <string>
 
namespace dawg {

template<std::size_t _N>
std::size_t key_switch(const std::string &ss, const std::string (&key)[_N]) {
	for(std::size_t i=0;i<_N;++i) {
		if(key[i].find(ss) == 0)
			return i;
	}
	return (size_t)-1;
}
 
} /* namespace dawg */
 
#endif /* DAWG_UTILS_H */
 