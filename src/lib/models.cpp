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
#if 0	
	// find the starting pair of sites.
	int g,m,n,m0;
	for(g=0; g < 64 && p[g] < 1.0; ++g)
		/*noop*/;
	for(n=0; n < 64 && p[(g+n)%64] >= 1.0; ++n)
		/*noop*/;
	// save the initial m, check to see if m-vector is empty
	m0 = g+n;
	n = m0+(n & 64);
	while(g < 64 && n < m0+64) {
		m = n%64;
		q[m] = uint64_bound(p[m]);
		a[m] = g;
		p[g] = (p[g]+p[m])-1.0;
		// find next sites as needed
		for(;g < 64 && p[g] < 1.0; ++g)
			/*noop*/;
		for(++n;n < m0+64 && p[n%64] >= 1.0; ++n)
			/*noop*/;
	}
	for(;g < 64; ++g) {
		if(p[g] >= 1.0) {
			q[g] = UINT64_C(18446744073709551615);
			a[g] = g;
		}
	}
	for(;n < m0+64; ++n) {
		g = n%64;
		if(p[g] < 1.0) {
			q[g] = UINT64_C(18446744073709551615);
			a[g] = g;
		}
	}
#else
	std::vector<int> sm,lg;
	for(int i=63;i>=0;--i) {
		if(p[i] < 1.0) {
			sm.push_back(i);
		} else {
			lg.push_back(i);
		}
	}
	while(!sm.empty() && !lg.empty()) {
		int m = sm.back(); sm.pop_back();
		int g = lg.back();
		q[m] = uint64_bound(p[m]);
		a[m] = g;
		p[g] = (p[g]+p[m])-1.0;
		if(p[g] < 1.0) {
			sm.push_back(g);
			lg.pop_back();
		}
	}
	for(int i=0;i<sm.size();++i) {
		int m = sm[i];
		q[m] = std::numeric_limits<boost::uint64_t>::max();
		a[m] = m;
	}
	for(int i=0;i<lg.size();++i) {
		int g = lg[i];
		q[g] = std::numeric_limits<boost::uint64_t>::max();
		a[g] = g;
	}
#endif	
}
 
bool dawg::subst_model::create_alias_tables() {
	alias_table_64(&freqs[0], &stat_dist_p[0], &stat_dist_a[0]);
	for(int i=0;i<64;++i) {
		printf("%u %llu %u\n", i, stat_dist_p[i], stat_dist_a[i]);
	}
	for(int i=0;i<64;++i)
		alias_table_64(&table[i][0], &mutation_p[i][0], &mutation_a[i][0]);
	return true;
}
