/****************************************************************************
 *  Copyright (C) 2009 Reed A. Cartwright, PhD <reed@scit.us>               *
 ****************************************************************************/

#include <boost/foreach.hpp>
#include <dawg/matic.h>
#include <dawg/bark.h>

#define foreach BOOST_FOREACH

template<class It>
bool has_intersection(It first1, It last1, It first2, It last2) {
	while(first1 != last1 && first2 != last2) {
		if(*first1 < *first2)
			++first1;
		else if(*first2 < *first1)
			++first2;
		else
			return true; // bark bad tree, then catch and add segment 
	}
	return false;
}

void dawg::matic::add_config_section(const dawg::ma &ma) {
	if(ma.root_segment >= configs.size())
		configs.resize(ma.root_segment+1);
	segment_info &seg = configs[ma.root_segment];

	// construction section_info
	std::auto_ptr<section_info> info(new section_info);
	
	// [snip]	
	
	// test whether descendents already exist in this segment
	foreach(section_info &r, seg) {
		if(has_intersection(
			r.node_names.begin(), r.node_names.end(),
			info->node_names.begin(), info->node_names.end())
			throw std::runtime_error("invalid trees; descendent already exists"0);
	}
	
	
	
	// find location to insert
	segment_info::iterator it;
	for(it = seg.begin(); it != seg.end()
	    && info->node_names.count(it->root_name); ++it)
		/*noop*/;
	seg.insert(it, info);
		
}
