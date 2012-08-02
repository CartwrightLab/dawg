#pragma once
#ifndef DAWG_GLOBAL_H
#define DAWG_GLOBAL_H
/****************************************************************************
 *  Copyright (C) 2010 Reed A. Cartwright, PhD <reed@scit.us>               *
 ****************************************************************************/

#include <string>
#include <vector>

#include <dawg/trick.h>

namespace dawg {

struct global_options {
#	define XM(name, type, def, desc) type XV(name) ;
#	include <dawg/details/global.xmh>
#	undef XM

	global_options() :
#	define XM(name, type, def, desc) XV(name) (def),
#	include <dawg/details/global.xmh>
#	undef XM
	_unused()
	{ }
	
	void read_section(const trick::data_type::value_type &sec) {
#	define XM(name, type, def, desc) sec.get(XP(name), XV(name));
#	include <dawg/details/global.xmh>
#	undef XM
	}
	
private:
	char _unused;
};

} // namespace dawg

#endif //DAWG_GLOBAL_H

