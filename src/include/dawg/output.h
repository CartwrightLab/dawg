#pragma once
#ifndef DAWG_OUTPUT_H
#define DAWG_OUTPUT_H
/****************************************************************************
 *  Copyright (C) 2010 Reed A. Cartwright, PhD <reed@scit.us>               *
 ****************************************************************************/

#include <ostream>
#include <fstream>

namespace dawg {

class output {
public:
	output() : p_out(NULL) { }
	
	bool open(const char *file_name);
	
protected:
	std::ostream *p_out;
	std::ofstream fout;
};

} // namespace dawg 

#endif // DAWG_OUTPUT_H

