#pragma once
#ifndef DAWG_MATIC_H
#define DAWG_MATIC_H
/****************************************************************************
 *  Copyright (C) 2009 Reed A. Cartwright, PhD <reed@scit.us>               *
 ****************************************************************************/

#include <boost/ptr_container/ptr_vector.hpp>

#include <dawg/ma.h>

#include <vector>
#include <set>

namespace dawg {
namespace details {
struct matic_section_info {
	std::set<std::string> node_names;
	std::string root_name;
};
}


// Core simulation algorithm class
class matic {
public:
	
	// Configure Simulation
	void add_config_section(const dawg::ma &ma);
	inline void configure(const dawg::ma &ma);
	inline void clear_configuration() {
		configs.clear();
	}

	
	template<class It>
	inline void configure(It first, It last) {
		clear_configuration();
		for(;first != last; ++first)
			add_config_section(*first);
	}

	// Run the simulation
	void walk();
	
protected:
	typedef dawg::details::matic_section_info section_info;
	typedef boost::ptr_vector<section_info> segment_info;
	typedef std::vector<section_info_vector> segment_info_vector;
	
	segment_info_vector configs;
};

}
#endif //DAWG_MATIC_H
