/****************************************************************************
 *  Copyright (C) 2009 Reed A. Cartwright, PhD <reed@scit.us>               *
 ****************************************************************************/
#include <dawg/utils/foreach.h>
#include <dawg/ma.h>
#include <dawg/log.h>
#include <dawg/wood.h>

#include <map>

using namespace dawg;
using namespace std;

// Read a section, setting and values that correspond to 
void dawg::ma::read_section(const pile::data_type::value_type &sec) {
	name = sec.name;
#	define XM(aname, type, def) sec.get(_P(aname), _V(aname));
#	include <dawg/details/dawgma.xmh>
#	undef XM
}

// Reads the pile format into a vector of dawg::mas
// Use inheritance and defaults
bool dawg::ma::from_pile(const pile &pyle, vector<ma> &v) {
	// create lookup map so we can know if something has been touched
	typedef map<std::string, const ma*> map_t;
	map_t lookup;
	// we reserve space to preserve pointers
	v.reserve(pyle.data.size());
	// create default object and put it in the map
	const ma def("_default_");
	lookup[def.name] = &def;
	
	for(pile::data_type::const_iterator secit = pyle.data.begin();
		secit != pyle.data.end(); ++secit) {
		// lookup parent section
		map_t::const_iterator iit = lookup.find(secit->inherits);
		if(iit == lookup.end())
			return DAWG_ERROR("section '" << secit->inherits <<
				"' not found (inherited by '" << secit->name << "')");
		// lookup section
		pair<map_t::iterator, bool> me = lookup.insert(make_pair(secit->name, (dawg::ma*)NULL));
		if(!me.second)
			return DAWG_ERROR("section '" << secit->name << "' specified more than once.");

		// read inherited ma
		v.push_back(*iit->second);
		// read the section
		v.back().read_section(*secit);
		// add section to map
		me.first->second = &v.back();
	}
	return true;
}
