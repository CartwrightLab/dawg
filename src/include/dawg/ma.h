#pragma once
#ifndef DAWG_MA_H
#define DAWG_MA_H
/****************************************************************************
 *  Copyright (C) 2009 Reed A. Cartwright, PhD <reed@scit.us>               *
 ****************************************************************************/

#include <string>
#include <vector>

#include <dawg/pile.h>

namespace dawg {

// dawg::ma is a "model argument" structure
struct ma {
#	define XM(name, type, def) type _V(name) ;
#	include <dawg/details/dawgma.xmh>
#	undef XM
	std::string name;

	ma(const std::string &_n = std::string() ) :
#	define XM(name, type, def) _V(name) (def),
#	include <dawg/details/dawgma.xmh>
#	undef XM
	name(_n)
	{ }
	
	static void from_pile(const dawg::pile &pyle, std::vector<dawg::ma> &v);
	
	void read_section(const pile::data_type::value_type &sec);
private:
};

}
#endif //DAWG_MA_H

