/****************************************************************************
 *  Copyright (C) 2009 Reed A. Cartwright, PhD <reed@scit.us>               *
 ****************************************************************************/
#include <dawg/details/foreach.h>
#include <dawg/ma.h>
#include <dawg/log.h>
#include <dawg/wood.h>

using namespace dawg;
using namespace std;

// Read a section, setting and values that correspond to 
void dawg::ma::read_section(const pile::map_type::value_type &sec) {
	name = sec.first;
#	define XM(name, type, def) sec.second.get(_P(name), _V(name));
#	include <dawg/details/dawgma.xmh>
#	undef XM
}

