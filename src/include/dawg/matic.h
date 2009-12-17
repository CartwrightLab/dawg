#pragma once
#ifndef DAWG_MATIC_H
#define DAWG_MATIC_H
/****************************************************************************
 *  Copyright (C) 2009 Reed A. Cartwright, PhD <reed@scit.us>               *
 ****************************************************************************/

#include <boost/ptr_container/ptr_vector.hpp>

#include <dawg/ma.h>
#include <dawg/subst.h>
#include <dawg/rate.h>
#include <dawg/indel.h>

#include <vector>
#include <set>

namespace dawg {
namespace details {
struct matic_section_info {
	std::set<std::string> node_names;
	std::string root_name;
	wood usertree;
	
	subst_model sub_mod;
	rate_model  rat_mod;
	indel_mix_model ins_mod;
	indel_mix_model del_mod;
	
};
}

// Core simulation algorithm class
class matic {
public:
	
	// Configure Simulation
	inline bool configure(const dawg::ma &ma) {
		clear_configuration();
		return add_config_section(ma);
	}
	inline void clear_configuration() {
		configs.clear();
	}
	
	template<class It>
	inline bool configure(It first, It last) {
		clear_configuration();
		for(;first != last; ++first)
			if(!add_config_section(*first))
				return DAWG_ERROR("Configuration section '" << first->name << "' failed to process.");
		return true;
	}

	// Run the simulation
	void walk();
	
protected:
	typedef dawg::details::matic_section_info section_info;
	typedef boost::ptr_vector<section_info> segment_info;
	typedef std::vector<segment_info> segment_info_vector;
	
	segment_info_vector configs;
	
	bool add_config_section(const dawg::ma &ma);	
};

}
#endif //DAWG_MATIC_H
