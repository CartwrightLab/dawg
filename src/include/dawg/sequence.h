#pragma once
#ifndef DAWG_SEQUENCE_H
#define DAWG_SEQUENCE_H
/****************************************************************************
 *  Copyright (C) 2010 Reed A. Cartwright, PhD <reed@scit.us>          *
 ****************************************************************************/
#include <dawg/residue.h>

namespace dawg {

namespace details {
struct aligned_sequence {
	std::string label;
	std::string seq;
};

} // details

typedef std::vector<residue> sequence;

struct alignment : public std::vector<details::aligned_sequence> {
	std::string::size_type max_label_width, max_label_width_14;
	int seq_type;
};

} // dawg

#endif // DAWG_SEQUENCE_H
