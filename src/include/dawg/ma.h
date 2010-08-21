#pragma once
#ifndef DAWG_MA_H
#define DAWG_MA_H
/****************************************************************************
 *  Copyright (C) 2009 Reed A. Cartwright, PhD <reed@scit.us>               *
 ****************************************************************************/

#include <string>
#include <vector>
#include <iostream>

#include <dawg/pile.h>
#include <dawg/utils/vecio.h>

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
	
	static bool from_pile(const dawg::pile &pyle, std::vector<dawg::ma> &v);
	
	void read_section(const pile::data_type::value_type &sec);
private:
};

template<class CharType, class CharTrait>
inline std::basic_ostream<CharType, CharTrait>& 
operator<<(std::basic_ostream<CharType, CharTrait>& o, const ma &a) {
	if(!o.good()) return o;
		
	o << set_open('\x7f') << set_close('\x7f') << set_delimiter(',');

	o << "[[ " << a.name << " ]]" << std::endl;
#	define XM(name, type, def) o << _P(name) " = " << a._V(name) << std::endl;
#	include <dawg/details/dawgma.xmh>
#	undef XM

	return o;
}

} //namespace dawg
#endif //DAWG_MA_H
