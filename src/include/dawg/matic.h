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
	typedef std::vector<dawg::ma> descriptions;
	
	// Configure Simulation
	bool configure(const descriptions& d);
	
	// Run the simulation
	void walk();
	
};

}
#endif //DAWG_MATIC_H

