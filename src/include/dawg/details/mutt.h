#pragma once
#ifndef DAWG_DETAILS_MUTT_H
#define DAWG_DETAILS_MUTT_H
/****************************************************************************
 *  Copyright (C) 2009 Reed A. Cartwright, PhD <reed@scit.us>               *
 ****************************************************************************/

// include stuff for DSFMT
#define DSFMT_MEXP 19937
#define DSFMT_DO_NOT_USE_OLD_NAMES
extern "C" {
#include <dawg/details/dSFMT.h>
}

namespace dawg { namespace details {

struct dsfmt_mutt_gen {
	double rand_01()  { return dsfmt_genrand_close_open(&state);	}
	double rand_01oc() { return dsfmt_genrand_open_close(&state); }
	double rand_01oo() { return dsfmt_genrand_open_open(&state); }
	uint32_t rand_uint32() { return dsfmt_genrand_uint32(&state); }
	void seed(uint32_t x) { dsfmt_init_gen_rand(&state, x); }
	template<int _N>
	void seed(uint32_t (&x)[_N]) {
		dsfmt_init_by_array(&state, &x[0], _N);
	}
	template<typename _It>
	void seed(_It first, _It last) {
		dsfmt_init_by_array(&state, &*first, last-first);
	}

private:
	dsfmt_t state;
};


}} /* namespace dawg::details */

#endif /* DAWG_DETAILS_MUTT_H */

