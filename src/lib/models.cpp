/****************************************************************************
 *  Copyright (C) 2012-2013 Reed A. Cartwright, PhD <reed@scit.us>          *
 ****************************************************************************/

#include <dawg/subst.h>
#include <cassert>
#include <cstdio>
#include <algorithm>

bool dawg::subst_model::create_alias_tables() {
	stat_dist_table_.create(&freqs[0], &freqs[64]);
	for(std::size_t k = 0; k < 64; ++k)
		mutation_table_[k].create(&table[k][0], &table[k][64]);
	return true;
}
