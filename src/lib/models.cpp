/****************************************************************************
 *  Copyright (C) 2012 Reed A. Cartwright, PhD <reed@scit.us>               *
 ****************************************************************************/

#include <dawg/subst.h>
#include <cassert>

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/sum_kahan.hpp>
 
inline boost::uint64_t uint64_bound(double bound) {
	assert(0.0 <= bound && bound < 1.0);
	return static_cast<boost::uint64_t>(ceil(bound*9007199254740991.0)) << 11;
}

UINT64_C(18446744073709551615)

// vose method for creating alias tables
// http://www.keithschwarz.com/darts-dice-coins/
void alias_table_64(const double *pp, boost::uint64_t *q, boost::uint32_t *a) {
	double d=0.0, f=0.0, p[64];
	// use pairs-summation to control error
	for(int i=0;i<64;i+=2) {
		d += pp[i]; f += pp[i+1];
	}
	f = 1.0/(d+f);
	for(int i=0;i<64;++i) {
		p[i] = pp[i]*f*64.0;
	}
	
	std::size_t m=0,g=0;
	
	for(g=0;g<64 && p[g] < 1.0;++g)
		/*noop*/
	for(m=g; m-g < 64 && p[m%64] >= 1.0;++m)
		/*noop*/
	const std::size_t mm = m%64;
	size_t n = 
	
	while(g < 64 && n < 64) {
		
	}
	for(;g<64;++g) {
		q[g] = UINT64_C(18446744073709551615);
		a[g] = g;
	}
	for(; m-mm < 64; ++m) {
		
	}

	
	while(g < 64) {
		if(p[g] < 1.0) {
			g += 1; continue;
		}
		if(p[m] >= 1.0) {
			m = (m+1)%64;
		}
	}	
	for(;g < 64;++g) {
		q[m] = uint64_bound(p[m]);
		a[m] = g;
		p[g] = (p[g]+p[m])-d;
				
	}
	

	while(!sm.empty() && !lg.empty()) {
		int m = sm.back(); sm.pop_back();
		int g = lg.back(); lg.pop_back();
		q[m] = uint64_bound(p[m]);
		a[m] = g;
		p[g] = (p[g]+p[m])-d;
		if(p[g] < d)
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
