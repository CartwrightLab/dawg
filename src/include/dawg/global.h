#pragma once
#ifndef DAWG_GLOBAL_H
#define DAWG_GLOBAL_H
/****************************************************************************
 *  Copyright (C) 2010 Reed A. Cartwright, PhD <reed@scit.us>               *
 ****************************************************************************/

#include <string>
#include <vector>

#include <dawg/pile.h>

namespace dawg {

struct global_options {
#	define XM(name, type, def) type _V(name) ;
#	include <dawg/details/global.xmh>
#	undef XM

	global_options() :
#	define XM(name, type, def) _V(name) (def),
#	include <dawg/details/global.xmh>
#	undef XM
	_unused()
	{ }
	
	void read_section(const pile::data_type::value_type &sec) {
#	define XM(name, type, def) sec.get(_P(name), _V(name));
#	include <dawg/details/global.xmh>
#	undef XM
	}
	
private:
	char _unused;
};

} // namespace dawg

#endif //DAWG_GLOBAL_H

