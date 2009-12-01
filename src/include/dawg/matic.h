#pragma once
#ifndef DAWG_MATIC_H
#define DAWG_MATIC_H
/****************************************************************************
 *  Copyright (C) 2009 Reed A. Cartwright, PhD <reed@scit.us>               *
 ****************************************************************************/

#include <dawg/ma.h>

#include <vector>

namespace dawg {

// Core simulation algorithm class
class matic {
public:
	// Configure Simulation
	void add_config_section(const dawg::ma &ma);
	void clear_configuration();
	inline void configure(const dawg::ma &ma);
	
	template<class It>
	inline void configure(It first, It last) {
		clear_configuration();
		for(;first != last; ++first)
			add_config_section(*first);
	}

	// Run the simulation
	void walk();
	
};

}
#endif //DAWG_MATIC_H

