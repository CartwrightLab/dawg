#pragma once
#ifndef DAWG_WOOD_H
#define DAWG_WOOD_H

/****************************************************************************
 *  Copyright (C) 2009 Reed A. Cartwright, PhD <reed@scit.us>               *
 ****************************************************************************/

#include <boost/scoped_ptr.hpp>

namespace dawg {

// Represents a phylogenetic tree structure
class wood {
public:
	struct node {
		boost::scoped_ptr<node> left, right;
		node *anc;
		std::string label;
		double length;		
	};
	
	boost::scoped_ptr<node> root;
};

}
#endif //DAWG_WOOD_H

