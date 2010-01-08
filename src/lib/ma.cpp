/****************************************************************************
 *  Copyright (C) 2009 Reed A. Cartwright, PhD <reed@scit.us>               *
 ****************************************************************************/
#include <dawg/details/foreach.h>
#include <dawg/ma.h>
#include <dawg/log.h>
#include <dawg/wood.h>

#include <map>

using namespace dawg;
using namespace std;

// Read a section, setting and values that correspond to 
void dawg::ma::read_section(const pile::data_type::value_type &sec) {
	name = sec.first;
#	define XM(name, type, def) sec.second.get(_P(name), _V(name));
#	include <dawg/details/dawgma.xmh>
#	undef XM
}

// Reads the pile format into a vector of dawg::mas
// Use inheritance and defaults
void dawg::ma::from_pile(const pile &pyle, vector<ma> &v) {
	// create lookup map so we can know if something has been touched
	typedef map<std::string, const ma*> map_t;
	map_t lookup;
	// we reserve space to preserve pointers
	v.reserve(pyle.data.size());
	// create default object and put it in the map
	const ma def("_default_");
	lookup.insert(make_pair(def.name, &def));
	
	for(pile::data_type::const_iterator secit = pyle.data.begin();
		secit != pyle.data.end(); ++secit) {
		// ensure that the inherited section exists
		map_t::const_iterator iit = 
		// TODO
		// read inherited ma
		v.push_back(*lookup[secit->second.inherits]);
		// read the section
		v.back().read_section(*secit);
		
	}
}
