/****************************************************************************
 *  Copyright (C) 2012 Reed A. Cartwright, PhD <reed@scit.us>               *
 ****************************************************************************/

#include <dawg/subst.h>
 
 // vose method for creating alias tables
 // http://www.keithschwarz.com/darts-dice-coins/
 void alias_table_64(const double *p, boost::uint64_t *q, boost::uint32_t *a) {
	double np[64];
	std::vector<int> sm, lg;
	double d = 0.0;
	for(int i=0;i<64;++i)
		d += p[i];
	d = 64.0/d;
	for(int i=0;i<64;++i) {
		q[i] = std::numeric_limits<boost::uint64_t>::max();
		a[i] = i;
		np[i] = p[i]*d;
		if(np[i] < 1.0) {
			sm.push_back(i);
		} else {
			lg.push_back(i);
		}
	}
	while(!sm.empty() && !lg.empty()) {
		int m = sm.back(); sm.pop_back();
		int g = lg.back(); lg.pop_back();
		q[m] = (static_cast<boost::uint64_t>(np[m]*9007199254740991.0) << 11) | 0x7FF;
		a[m] = g;
		np[g] = (np[g]+np[m])-1.0;
		if(np[g] < 1.0)
			sm.push_back(g);
		else
			lg.push_back(g);
	}
 }
 
 bool dawg::subst_model::create_alias_tables() {
	alias_table_64(&freqs[0], &stat_dist_p[0], &stat_dist_a[0]);
	for(int i=0;i<64;++i)
		alias_table_64(&table[i][0], &mutation_p[i][0], &mutation_a[i][0]);
	return true;
 }